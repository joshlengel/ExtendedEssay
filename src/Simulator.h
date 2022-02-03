#pragma once

#include"Vec.h"

#include<vector>

class Body;
class SolarSystem;

class EulerSimulator
{
public:
    Body *rocket;
    SolarSystem &solar_system;

    EulerSimulator(double t, Body *rocket, SolarSystem &solar_system);

    void Step(double gamma, double &dt);

private:
    double t;
    Vec3 prev_position, prev_f;
};