#pragma once
#include <functional>
#include "ParticlesTable.hpp"
#include "SimulationSetup.hpp"

struct Integrator {
    std::function<void()> calculate_a;
    std::function<void(float)> kick;
    std::function<void(float)> drift;
};

/*
    Integrator wrap_Integrator(ParticlesTable* pt, SimulationSetup* simsetup);

Wrap integrator into an `Integrator` structure

### Input 
    - ParticlesTable* pt: 
    - SimulationSetup* simsetup
*/
Integrator wrap_Integrator(ParticlesTable* pt, SimulationSetup* simsetup);
