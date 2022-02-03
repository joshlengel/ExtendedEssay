#include"Body.h"
#include"Simulator.h"
#include"Lambert.h"
#include"Porkchop.h"

#include<SpiceUsr.h>
#include<nlohmann/json.hpp>

#include<iostream>
#include<fstream>

#define JSON_OUTPUT_DIR "../out/"

static const char *const KERNEL_PATH = "../assets/solar_system.tm";

static const double G = 6.674e-11;

static double T_INJECT;
static double T_EXIT;

int main()
{
    nlohmann::json dst;

    // Setup CSPICE data
    furnsh_c(KERNEL_PATH);

    T_INJECT = GetET("28 Sep 2022 00:00:00");
    T_EXIT = GetET("01 Jun 2023 00:00:00");

    // Setup solar system
    SolarSystem solar_system(T_INJECT);

    double r_leo = 1e6 + 6e6;
    double r_lmo = 1e6 + 5e6;

    // Porkchop plot
    double DT_START = GetET("2022 Apr 1 00:00:00");
    double DT_STOP = GetET("2022 Dec 1 00:00:00");

    double TOF_START = 120 * 24 * 3600;
    double TOF_STOP = 500 * 24 * 3600;

    PorkchopPlot porkchop_plot(DT_START, DT_STOP, TOF_START, TOF_STOP, 50, r_leo, r_lmo);
    porkchop_plot.Write(dst);

    // Setup mission
    Vec3 mars_exit_position, mars_exit_velocity;
    GetEphemeris("Mars Barycenter", T_EXIT, mars_exit_position, mars_exit_velocity);

    Vec3 v1, v2;
    LambertUV(solar_system.earth.position - solar_system.sun.position, mars_exit_position - solar_system.sun.position, T_EXIT - T_INJECT, solar_system.sun.mu, v1, v2);
    Vec3 v_inject = v1 + solar_system.sun.velocity - solar_system.earth.velocity;
    Vec3 v_exit = mars_exit_velocity - (v2 + solar_system.sun.velocity);
    double C3_e = Vec3::Dot(v_inject, v_inject);
    double C3_m = Vec3::Dot(v_exit, v_exit);
    double dv1 = sqrt(C3_e + 2 * solar_system.earth.mu / r_leo) - sqrt(solar_system.earth.mu / r_leo);
    double dv2 = sqrt(C3_m + 2 * solar_system.mars.mu / r_lmo) - sqrt(solar_system.mars.mu / r_lmo);

    Vec3 n_earth = -Vec3::Normalize(Vec3(0.0, 1.0, -v_inject.y / v_inject.z));
    double phi = asin(1.0 / (1.0 + r_leo * C3_e / solar_system.earth.mu));
    Vec3 v_n = Vec3::Normalize(Vec3::Rotate(v_inject, n_earth, phi));
    Vec3 v_i = v_n * sqrt(C3_e + 2 * solar_system.earth.mu / r_leo) + solar_system.earth.velocity;
    Vec3 x_n = Vec3::Normalize(Vec3::Cross(n_earth, v_n));
    Vec3 x_i = x_n * r_leo + solar_system.earth.position;

    std::cout << "v_n: (" << v_n.x << ", " << v_n.y << ", " << v_n.z << ")" << std::endl;
    std::cout << "x_n: (" << x_n.x << ", " << x_n.y << ", " << x_n.z << ")" << std::endl;

    double rocket_mass = 100000.0;
    Body rocket("Rocket", x_i, v_i, rocket_mass * G);

    // Setup output
    double t = T_INJECT;
    EulerSimulator sim(t, &rocket, solar_system);

    double gamma = 2e-2;
    dst["metadata"]["gamma"] = gamma;
    dst["metadata"]["C3_e"] = C3_e;
    dst["metadata"]["C3_m"] = C3_m;
    dst["metadata"]["dv1"] = dv1;
    dst["metadata"]["dv2"] = dv2;
    dst["metadata"]["dv"] = dv1 + dv2;
    dst["metadata"]["phi"] = phi;
    dst["metadata"]["n_earth"]["x"] = n_earth.x;
    dst["metadata"]["n_earth"]["y"] = n_earth.y;
    dst["metadata"]["n_earth"]["z"] = n_earth.z;

    bool exited = false;

    while (t < T_EXIT + 100 * 24 * 3600)
    {
        solar_system.Output(dst);
        rocket.Output(dst);
        dst["time"].push_back(t - T_INJECT);

        double dt;
        sim.Step(gamma, dt);
        t += dt;

        Vec3 mars_disp = rocket.position - solar_system.mars.position;
        double mars_dist_sqr = Vec3::Dot(mars_disp, mars_disp);
        double max_dist = 3e8;

        if (!exited && mars_dist_sqr < max_dist * max_dist)
        {
            double mars_dist = sqrt(mars_dist_sqr);
            Vec3 rocket_vel = rocket.velocity - solar_system.mars.velocity;
            double cosa = Vec3::Dot(mars_disp, rocket_vel) / (mars_dist * Vec3::Length(rocket_vel));
            if (fabs(cosa) < 0.005f)
            {
                exited = true;
                std::cout << "Exit orbit. Distance to mars: " << sqrt(mars_dist_sqr) << std::endl;

                // Do exit orbit
                Vec3 nm = Vec3::Normalize(Vec3::Cross(mars_disp, rocket_vel));
                double v_or = sqrt(solar_system.mars.mu / mars_dist);
                rocket_vel = Vec3::Normalize(Vec3::Cross(nm, mars_disp)) * v_or;
                rocket.velocity = rocket_vel + solar_system.mars.velocity;

                dst["metadata"]["n_mars"]["x"] = nm.x;
                dst["metadata"]["n_mars"]["y"] = nm.y;
                dst["metadata"]["n_mars"]["z"] = nm.z;
            }
        }
    }

    std::cout << "Exit orbit. Distance to mars: " << Vec3::Length(rocket.position - solar_system.mars.position) << std::endl;

    // Write dst to file
    std::ofstream out(JSON_OUTPUT_DIR "results.json");
    out << dst;
}