#include"Simulator.h"
#include"Body.h"

#include<algorithm>
#include<functional>

RK4Simulator::RK4Simulator(double t, Body *rocket, SolarSystem &solar_system):
    rocket(rocket),
    solar_system(solar_system),
    t(t)
{}

void RK4(const std::function<Vec3(double t, const Vec3&)> &f, double &t, Vec3 &x, Vec3 &v, double h)
{
    double h2 = h * 0.5;

    Vec3 k1 = v;
    Vec3 l1 = f(t, x);
    Vec3 k2 = v + l1 * h2;
    Vec3 l2 = f(t + h2, x + k1 * h2);
    Vec3 k3 = v + l2 * h2;
    Vec3 l3 = f(t + h2, x + k2 * h2);
    Vec3 k4 = v + l3 * h;
    Vec3 l4 = f(t + h, x + k4 * h);

    t = t + h;
    x = x + (k1 + k2 * 2 + k3 * 2 + k4) * (h / 6.0);
    v = v + (l1 + l2 * 2 + l3 * 2 + l4) * (h / 6.0);
}

void RK4Simulator::Step(double gamma, double &dt)
{
    auto GetGravity = [this](const Body &body, double t, const Vec3 &x)
    {
        Vec3 position, velocity;
        GetEphemeris(body.name, t, position, velocity);
        Vec3 r = position - x;
        return (r / Vec3::Length(r)) * (body.mu / Vec3::Dot(r, r));
    };

    auto f = [&, this](double t, const Vec3 &x)
    {
        Vec3 a;
        a += GetGravity(solar_system.sun, t, x);
        a += GetGravity(solar_system.earth, t, x);
        a += GetGravity(solar_system.mars, t, x);
        return a;
    };

    rocket->acceleration = f(t, rocket->position);

    dt = std::min(gamma / Vec3::Length(rocket->acceleration), 20.0);
    solar_system.Step(dt);

    // We have the SODE x'' = sum GM/|x_i - x|^3 * (x_i - x)
    // We can introduce the variable v = x' to obtain two FODE:
    // v' = sum GM/|x_i - x|^3 * (x_i - x) = f1(x)
    // x' = v

    // Runge-Kutta-4 iteration
    RK4(f, t, rocket->position, rocket->velocity, dt);
}