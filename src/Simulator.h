#pragma once

#include"Vec.h"

#include<vector>

class Body;
class SolarSystem;

class RK4Simulator
{
public:
    Body *rocket;
    SolarSystem &solar_system;

    RK4Simulator(double t, Body *rocket, SolarSystem &solar_system);

    void Step(double gamma, double &dt);

private:
    double t;
    Vec3 prev_position, prev_f;
};