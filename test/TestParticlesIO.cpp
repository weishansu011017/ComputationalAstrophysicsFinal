#include "ParticlesTable.hpp"
#include "UnitsTable.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <filesystem>


/*
    TEST: Testing the HDF5 I/O funcutionality

    Test content: 
    1. Initializing an empty `ParticlesTable` class
    2. Extracting the `ParticlesTable` into HDF5
    3. Reading the `ParticlesTable` from HDF5
*/
int main(){

    // Randomly setup code unit (Doesn't matter)
    UnitsTable unit;
    ParticlesTable pt(unit, 10000);

    // Assign mass to 10000
    for (int i = 0; i < pt.N; ++i) {
        pt.m[i] = 1.0f;
    }

    // I/O the file
    pt.extract_particles_table("TEST.hdf5");
    ParticlesTable rpt = ParticlesTable::read_particles_table("TEST.hdf5");

    // Verifing the data
    float msum1 = 0, msum2 = 0;
    for (int i = 0; i < pt.N; ++i) {
        msum1 += pt.m[i];
        msum2 += rpt.m[i];
    }
    assert(std::abs(msum1 - msum2) < 1e-5f);
    std::cout << "HDF5 read/write test passed! Total mass: " << msum1 << std::endl;
    std::filesystem::remove("TEST.hdf5");
    return 0;
}