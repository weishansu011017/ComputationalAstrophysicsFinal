#include <random>
#include <chrono>
#include <functional>
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
        Mtot += mi;
    }
    pt.Mtot = Mtot;
    return pt;
}

