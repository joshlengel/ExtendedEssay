import spiceypy as sp
import numpy as np
from spiceypy.spiceypy import furnsh, str2et, unload

def init():
    furnsh("data/solar_system.tm")

def terminate():
    unload("data/solar_system.tm")

def get_time(str):
    return str2et(str)

range = [0,1,2]

class Entity:
    def __init__(self, x, v, mu):
        self.x = x
        self.v = v
        self.a = np.zeros(3)
        self.mu = mu
    
    def step(self, dt):
        self.x = self.x + self.v * dt
        self.v = self.v + self.a * dt
        self.a.fill(0)

class Body(Entity):
    def __init__(self, name, time):
        Entity.__init__(self, np.zeros(3), np.zeros(3), 0)
        self.name = name
        self.time = time
        self._get_state()
        self.is_static = False
    
    def _get_state(self):
        state, _ = sp.spkezr(self.name, self.time, "J2000", "NONE", "SSB")
        _, val = sp.bodvrd(self.name, "GM", 1)
        np.put(self.x, range, state[:3] * 1000)
        np.put(self.v, range, state[3:] * 1000)
        self.mu = val * 1000000000

    def step(self, dt):
        if self.is_static:
            return

        self.time = self.time + dt
        self._get_state()


class Propagator:
    def __init__(self, *bodies):
        self.bodies = bodies
    
    def step(self, dt):
        for b1 in self.bodies:
            if b1 is Body:
                continue
            
            for b2 in self.bodies:
                if b1 is b2:
                    continue

                r = b2.x - b1.x
                R = np.linalg.norm(r)
                a = r * b2.mu / (R ** 3)
                b1.a += a
    
        for b in self.bodies:
            b.step(dt)