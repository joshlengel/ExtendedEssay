import math
import numpy as np

def C2(z):
    if z > 0:
        return (1 - math.cos(math.sqrt(z))) / z
    elif z < 0:
        return (math.cosh(math.sqrt(-z)) - 1) / -z
    else:
        return 0.5

def C3(z):
    sqrt = math.sqrt(abs(z))
    sqrt3 = z * sqrt
    if z > 0:
        return (sqrt - math.sin(sqrt)) / sqrt3
    elif z < 0:
        return (math.sinh(sqrt) - sqrt) / sqrt3
    else:
        return 1/6

class Lambert:
    def __init__(self, r1, r2, tof, mu, itr=100, tol=1e-6):
        # R1 = np.linalg.norm(r1)
        # R2 = np.linalg.norm(r2)
        # gamma = np.dot(r1, r2) / (R1 * R2)
        # A = math.sqrt(R1 * R2 * (1 + gamma))

        # self.possible = False
        # if A == 0:
        #     return
        
        # phi = 0
        # phi_l = -4 * math.pi * math.pi
        # phi_u =  4 * math.pi * math.pi
        
        # for _ in range(itr):
        #     c2 = C2(phi)
        #     c3 = C3(phi)

        #     B = R1 + R2 + A * (phi * c3 - 1) / math.sqrt(c2)

        #     if A > 0 and B < 0:
        #         B = -B
        #         phi_l += math.pi
            
        #     chi = math.sqrt(B / c2)
        #     tof_ = (chi * chi * chi * c3 + A * math.sqrt(B)) / math.sqrt(mu)

        #     if abs(tof_ - tof) < tol:
        #         self.possible = True
        #         break

        #     if tof_ < tof: phi_l = phi
        #     else: phi_u = phi

        #     phi = 0.5 * (phi_l + phi_u)
        
        # if not self.possible:
        #     return
        
        # F = 1 - B / R1
        # G = A * math.sqrt(B / mu)
        # G_dot = 1 - B / R2

        # self.v1 = (r2 - r1 * F) / G
        # self.v2 = (r2 * G_dot - r1) / G

        self.possible = True

        #r1 = np.array([r1[0], r1[1], 0])
        #r2 = np.array([r2[0], r2[1], 0])
        c = r2 - r1
        C = np.linalg.norm(c)
        R1 = np.linalg.norm(r1)
        R2 = np.linalg.norm(r2)
        s = 0.5 * (R1 + R2 + C)
        ir1 = r1 / R1
        ir2 = r2 / R2
        ih = np.cross(ir1, ir2)
        ih = ih / np.linalg.norm(ih)
        ll = 1 - C / s
        l = math.sqrt(ll)
        lll = ll * l

        it1, it2 = np.array(3), np.array(3)

        if ih[2] < 0:
            l = -l
            it1 = np.cross(ir1, ih)
            it2 = np.cross(ir2, ih)
        else:
            it1 = np.cross(ih, ir1)
            it2 = np.cross(ih, ir2)
        
        it1 = it1 / np.linalg.norm(it1)
        it2 = it2 / np.linalg.norm(it2)
        
        T = math.sqrt(2 * mu / (s * s * s)) * tof

        M_max = math.floor(T / math.pi)
        T00 = math.acos(l) + l * math.sqrt(1 - ll)
        T0 = (T00 + M_max * math.pi)
        T1 = 2.0 / 3.0 * (1 - lll)

        DT, DDT, DDDT = 0.0, 0.0, 0.0

        if M_max > 0:
            if T < T0:
                it = 0
                err = 1.0
                T_min = T0
                x_old, x_new = 0.0, 0.0

                while True:
                    DT, DDT, DDDT = self.__dT_dx(x_old, T_min, ll, lll)

                    if (DT != 0.0):
                        x_new = x_old - DT * DDT / (DDT * DDT - 0.5 * DT * DDDT)

                    err = abs(x_old - x_new)

                    if err < 1e-13 or it > 12:
                        break

                    T_min = self.__x2_tof(x_new, M_max, l, ll)
                    x_old = x_new
                    it += 1

                if T_min > T:
                    M_max -= 1
        
        x = 0.0

        if T >= T00:
            x = -(T - T00) / (T - T00 + 4)
        elif T <= T1:
            x = 2.5 * T1 * (T1 - T) / ((1 - ll * lll) * T) + 1
        else:
            x = math.pow(T / T00, 0.69314718055994529 / math.log(T1 / T00)) - 1

        x = self.__householder(T, x, 1e-5, 15, l, ll, lll)

        gamma = math.sqrt(0.5 * mu * s)
        rho = (R1 - R2) / C
        sigma = math.sqrt(1 - rho * rho)

        y = math.sqrt(1 - ll + ll * x * x)
        arg1 = l * y - x
        arg2 = l * y + x
        arg3 = y + l * x
        Vr1 = gamma * (arg1 - rho * arg2) / R1
        Vr2 = -gamma * (arg1 + rho * arg2) / R2
        Vt = gamma * sigma * arg3
        Vt1 = Vt / R1
        Vt2 = Vt / R2
        self.v1 = Vr1 * ir1 + Vt1 * it1
        self.v2 = Vr2 * ir2 + Vt2 * it2

    def __householder(self, T, x0, eps, it_max, l, ll, lll):
        it = 0
        err = 1
        x_new = 0
        tof = 0
        delta = 0
        DT, DDT, DDDT = 0.0, 0.0, 0.0

        while err > eps and it < it_max:
            tof = self.__x2_tof(x0, 0, l, ll)
            DT, DDT, DDDT = self.__dT_dx(x0, tof, ll, lll)
            delta = tof - T
            DT2 = DT * DT
            x_new = x0 - delta * (DT2 - 0.5 * delta * DDT) / (DT * (DT2 - delta * DDT) + DDDT * delta * delta / 6)
            err = abs(x0 - x_new)
            x0 = x_new
            it += 1
        
        return x0

    def __dT_dx(self, x, T, ll, lll):
        umx2 = 1 - x * x
        y = math.sqrt(1 - ll * umx2)
        yy = y * y
        yyy = yy * y
        
        DT = 1 / umx2 * (3 * T * x - 2 + 2 * lll * x / y)
        DDT = 1 / umx2 * (3 * T + 5 * x * DT + 2 * (2 - ll) * lll / yyy)
        DDDT = 1 / umx2 * (7 * x * DDT + 8 * DT - 6 * (1 - ll) * ll * lll * x / (yyy * yy))

        return DT, DDT, DDDT

    def __x2_tof(self, x, N, l, ll):
        battin = 0.01
        lagrange = 0.2
        dist = abs(x - 1)
        if dist < lagrange and dist > battin:
            return self.__x2_tof2(x, N, l, ll)
        
        K = ll
        E = x * x - 1
        rho = abs(E)
        z = math.sqrt(1 + K * E)

        if dist < battin:
            eta = z - l * x
            S1 = 0.5 * (1 - l - x * eta)
            Q = 4.0 / 3.0 * self.__hypergeometricF(S1, 1e-11)
            return 0.5 * (eta * eta * eta * Q + 4 * l * eta) + N * math.pi / pow(rho, 1.5)
        else:
            y = math.sqrt(rho)
            g = x * z - l * E
            d = 0

            if E < 0:
                d = N * math.pi + math.acos(g)
            else:
                d = math.log(y * (z - l * x) + g)
            
            return (x - l * z - d / y) / E

    def __x2_tof2(self, x, N, l, ll):
        a = 1 / (1 - x * x)

        if a > 0:
            alpha = 2 * math.acos(x)
            beta = 2 * math.asin(math.sqrt(ll / a))

            if l < 0:
                beta = -beta
            
            return 0.5 * a * math.sqrt(a) * ((alpha - math.sin(alpha)) - (beta - math.sin(beta)) + 2 * math.pi * N)
        else:
            alpha = 2 * math.acosh(x)
            beta = 2 * math.asinh(math.sqrt(-ll / a))

            if l < 0:
                beta = -beta

            return 0.5 * -a * math.sqrt(-a) * ((beta - math.sinh(beta)) - (alpha - math.sinh(alpha)))

    def __hypergeometricF(self, z, tol):
        Sj, Cj = 1.0, 1.0
        err = 1
        Cj1, Sj1 = 0.0, 0.0
        j = 0

        while err > tol:
            Cj1 = Cj * z * (3 + j) * (1 + j) / ((2.5 + j) * (1 + j))
            Sj1 = Sj + Cj1
            err = abs(Cj1)
            Sj = Sj1
            Cj = Cj1
            j += 1

        return Sj