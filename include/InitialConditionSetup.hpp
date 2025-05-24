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

/*
    std::string format_index(int index, int width = 5)

Format an integer into a zero-padded string, e.g., 42 → "00042".

## Input 
    - int index : The integer to be formatted
    - int width : Total width of the resulting string (default = 5)

## Output
    - std::string : Zero-padded string representation of the input integer
*/
std::string format_index(int index, int width = 5);