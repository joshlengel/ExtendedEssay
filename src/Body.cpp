#include"Body.h"

#include<SpiceUsr.h>

void GetEphemeris(const std::string &name, double et, Vec3 &position, Vec3 &velocity)
{
    double state[6];
    double lt;
    spkezr_c(name.c_str(), et, "J2000", "NONE", "SSB", state, &lt);
    position = Vec3(state[0], state[1], state[2]) * 1e3;
    velocity = Vec3(state[3], state[4], state[5]) * 1e3;
}
double GetMu(const std::string &name) { int dim; double mu; bodvrd_c(name.c_str(), "GM", 1, &dim, &mu); return mu * 1e9; }

double GetET(const std::string &name) { double et; utc2et_c(name.c_str(), &et); return et; }

Body::Body(const std::string &name, const Vec3 &position, const Vec3 &velocity, double mu):
    name(name),
    position(position),
    velocity(velocity),
    acceleration(),
    mu(mu)
{}

Body::Body(const std::string &name, double et):
    name(name),
    acceleration(),
    mu(GetMu(name))
{
    GetEphemeris(name.c_str(), et, position, velocity);
}

void Body::Output(nlohmann::json &json) const
{
    nlohmann::json &dst = json["bodies"][name]["ephemerides"];
    dst.push_back(position.x); dst.push_back(position.y); dst.push_back(position.z);
}

SolarSystem::SolarSystem(double et):
    et(et),
    sun("Sun", et),
    earth("Earth", et),
    mars("Mars Barycenter", et)
{}

void SolarSystem::Step(double dt)
{
    et += dt;
    GetEphemeris(sun.name, et, sun.position, sun.velocity);
    GetEphemeris(earth.name, et, earth.position, earth.velocity);
    GetEphemeris(mars.name, et, mars.position, mars.velocity);
}

void SolarSystem::Output(nlohmann::json &json) const
{
    sun.Output(json);
    earth.Output(json);
    mars.Output(json);
}