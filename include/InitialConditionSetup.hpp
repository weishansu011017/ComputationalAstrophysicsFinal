#pragma once
#include "ParticlesSetup.hpp"
#include "UnitsTable.hpp"
#include "ParticlesTable.hpp"

/*
    ParticlesTable setup_initial_condition(ParticlesSetup setup, UnitsTable unit)
Initilizing the initial condition of particles from `ParticlesSetup`

## Input 
    - ParticlesSetup setup: The parameter of setup 
    - UnitsTable unit: The code units container.

## Output
    - `ParticlesTable`: The table of particles.
*/
ParticlesTable setup_initial_condition(const ParticlesSetup& setup, UnitsTable unit);


