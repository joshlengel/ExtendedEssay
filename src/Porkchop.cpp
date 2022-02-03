#include"Porkchop.h"
#include"Lambert.h"
#include"Vec.h"

#include<SpiceUsr.h>

PorkchopPlot::PorkchopPlot(double dt1, double dt2, double tof1, double tof2, size_t num_values, double rleo, double rlmo):
    m_dt1(dt1), m_dt2(dt2),
    m_tof1(tof1), m_tof2(tof2),
    m_num_values(num_values)
{
    m_values.reserve(num_values * num_values);

    int dim;
    double gm;

    bodvrd_c("Sun", "GM", 1, &dim, &gm);
    double mu_sun = gm * 1e9;

    bodvrd_c("Earth", "GM", 1, &dim, &gm);
    double mu_earth = gm * 1e9;

    bodvrd_c("Mars Barycenter", "GM", 1, &dim, &gm);
    double mu_mars = gm * 1e9;

    for (size_t i = 0; i < num_values; ++i)
    for (size_t j = 0; j < num_values; ++j)
    {
        double t1 = (dt2 - dt1) * static_cast<float>(i) / num_values + dt1;
        double tof = (tof2 - tof1) * static_cast<float>(j) / num_values + tof1;
        double t2 = t1 + tof;

        double lt;

        double earth_state[6];
        spkezr_c("Earth", t1, "J2000", "NONE", "Sun", earth_state, &lt);
        Vec3 earth_pos = Vec3(earth_state[0], earth_state[1], earth_state[2]) * 1000;
        Vec3 earth_vel = Vec3(earth_state[3], earth_state[4], earth_state[5]) * 1000;

        double mars_state[6];
        spkezr_c("Mars Barycenter", t2, "J2000", "NONE", "Sun", mars_state, &lt);
        Vec3 mars_pos = Vec3(mars_state[0], mars_state[1], mars_state[2]) * 1000;
        Vec3 mars_vel = Vec3(mars_state[3], mars_state[4], mars_state[5]) * 1000;

        Vec3 v1, v2;
        LambertDIzzo(earth_pos, mars_pos, tof, mu_sun, v1, v2);

        Vec3 rel_1 = v1 - earth_vel;
        Vec3 rel_2 = mars_vel - v2;
        double C3_e = Vec3::Dot(rel_1, rel_1);
        double C3_m = Vec3::Dot(rel_2, rel_2);
        double dv_1 = sqrt(C3_e + 2 * mu_earth / rleo) - sqrt(mu_earth / rleo);
        double dv_2 = sqrt(C3_m + 2 * mu_mars / rlmo) - sqrt(mu_mars / rlmo);
        m_values.push_back(dv_1 + dv_2);
    }
}

void PorkchopPlot::Write(nlohmann::json &dst)
{
    dst["metadata"]["porkchop"]["num_values"] = m_num_values;
    dst["metadata"]["porkchop"]["dt1"] = m_dt1;
    dst["metadata"]["porkchop"]["dt2"] = m_dt2;
    dst["metadata"]["porkchop"]["tof1"] = m_tof1;
    dst["metadata"]["porkchop"]["tof2"] = m_tof2;

    for (double value : m_values)
        dst["porkchop"]["values"].push_back(value);
}