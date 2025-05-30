#include <random>
#include <chrono>
#include <functional>
#include <sstream>
#include <iomanip>
#include "InitialConditionSetup.hpp"
#include "ParticlesSetup.hpp"
#include "UnitsTable.hpp"
#include "ParticlesTable.hpp"
#include "SamplingFunctionsSet.hpp"


ParticlesTable setup_initial_condition(const ParticlesSetup& setup, UnitsTable unit){
    // Varification 
    SamplingFunctionsSet samplers = setup.get_sampler();

    // Initializing container
    ParticlesTable pt(unit, setup.N);
    pt.SimulationTag = setup.SimulationTag;

    // Sampling particles properties
    // Initialize h
    float happrox = setup.softfactorx * setup.simulation_scale * pow(pt.N, -1.0/3.0);

    float Mtot = 0;
    for (int i = 0; i < pt.N; ++i) {
        pt.particle_index[i] = static_cast<uint32_t>(i + 1);

        auto coor = samplers.coorsampler();

        pt.x[i] = coor[0];
        pt.y[i] = coor[1];
        pt.z[i] = coor[2];

        pt.vx[i] = coor[3];
        pt.vy[i] = coor[4];
        pt.vz[i] = coor[5];

        float mi = samplers.msampler();
        pt.m[i] = mi;
        pt.h[i] = happrox;
        Mtot += mi;

    }
    pt.Mtot = Mtot;

    // Initialize dt
    float max_vel = 0.0f;
    for (int i = 0; i < pt.N; ++i){
        float vi = std::sqrt(pt.vx[i]*pt.vx[i] + pt.vy[i]*pt.vy[i] + pt.vz[i]*pt.vz[i]);
        if (vi > max_vel) max_vel = vi;
    }
    float h_mean = mean(pt.h);
    float dt_vel = setup.tsfactor * h_mean / std::max(max_vel, 1e-8f);
    for (int i = 0; i < pt.N; ++i) {
        pt.dt[i] = dt_vel;
    }

    return pt;
}

std::string format_index(int index, int width) {
    std::ostringstream oss;
    oss << std::setw(width) << std::setfill('0') << index;
    return oss.str();
}