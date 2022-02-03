#pragma once

#include<nlohmann/json.hpp>

#include<vector>

class PorkchopPlot
{
public:
    PorkchopPlot(double dt1, double dt2, double tof1, double tof2, size_t num_values, double rleo, double rlmo);

    void Write(nlohmann::json &dst);

private:
    double m_dt1, m_dt2, m_tof1, m_tof2;
    size_t m_num_values;
    std::vector<double> m_values;
};