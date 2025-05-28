#include "ParticlesTable.hpp"
#include "UnitsTable.hpp"
#include "ParticlesSetup.hpp"
#include "InitialConditionSetup.hpp"
#include "SimulationSetup.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <filesystem>


/*
    TEST: Testing parameter file readed

    Test content: 
    1. Initializing a `ParticlesSetup`
    2. Setup the `ParticlesTable` by uniform distribution
    3. Writing the `ParticlesTable` to HDF5
*/
int main(){

    // Randomly setup code unit (Doesn't matter)
    std::string tag = "setup";
    UnitsTable unit;
    ParticlesSetupUniform setup(tag);
    ParticlesTable pt = setup_initial_condition(setup, unit);
    

    // Make params file
    SimulationSetup::generate_parameters_file(setup, 0.01);   
    SimulationSetup simsetup = SimulationSetup(std::string("setup.in"));
    simsetup.input_file = "awfklhfshkasflkhaefhkljfae.h5";
    simsetup.paramspath = "wdawawdadw.in";
    simsetup.make_parameters_file();

    return 0;
}