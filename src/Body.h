#pragma once

#include"Vec.h"

#include<nlohmann/json.hpp>

#include<string>

void GetEphemeris(const std::string &name, double et, Vec3 &position, Vec3 &velocity);
double GetMu(const std::string &name);

double GetET(const std::string &name);

class Body
{
public:
    std::string name;
    Vec3 position, velocity, acceleration;
    double mu;

    Body(const std::string &name, const Vec3 &position, const Vec3 &velocity, double mu);
    Body(const std::string &name, double et);

    void Output(nlohmann::json &json) const;
};

class SolarSystem
{
public:
    double et;

    Body sun;
    Body earth;
    Body mars;
    
    SolarSystem(double et);
    
    void Step(double dt);
    void Output(nlohmann::json &json) const;
};