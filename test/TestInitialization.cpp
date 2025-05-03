#include "ParticlesTable.hpp"
#include "UnitsTable.hpp"
#include "ParticlesSetup.hpp"
#include "InitialConditionSetup.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <filesystem>


/*
    TEST: Testing particle initilization

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


    // I/O the file
    pt.extract_particles_table(pt.SimulationTag + ".h5");

    // // Verifing the data
    // float msum1 = 0, msum2 = 0;
    // for (int i = 0; i < pt.N; ++i) {
    //     msum1 += pt.m[i];
    //     msum2 += rpt.m[i];
    // }
    // assert(std::abs(msum1 - msum2) < 1e-5f);
    // std::cout << "HDF5 read/write test passed! Total mass: " << msum1 << std::endl;
    // std::filesystem::remove("TEST.hdf5");
    return 0;
}