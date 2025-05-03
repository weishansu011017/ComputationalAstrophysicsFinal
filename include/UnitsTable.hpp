#pragma once
#include <hdf5.h>
#include <cmath>
#include "ParticlesSetup.hpp"
#include "PhysicalConstants.hpp"
#include <iostream>

/*
Table of unit for setup
*/
class UnitsTable {
public:
    float udist;
    float utime;
    float umass;

    // Constructor
    UnitsTable(float _udist = 1.0f, float _umass = 1.0f)
        : udist(_udist), umass(_umass)
    {
        if (_udist <= 0.0f || _umass <= 0.0f)
            throw std::invalid_argument(
                "Code units must be positive. Received udist = " + std::to_string(_udist)
                + ", umass = " + std::to_string(_umass));

        utime = std::sqrt(1.0f / (umass * CGSConstants::G / udist / udist / udist ));
    }

    // Destructor
    ~UnitsTable() = default;

    /*
        void read_unit_HDF5(hid_t file_id);
    Read the code unit setup in HDF5

    ## Input
        - hid_t file_id: The HDF5 file handle created by H5Fcreate or H5Fopen
    */
    void read_unit_HDF5(hid_t file_id);

    /*
        UnitsTable setup_units(const ParticlesSetup& setup);
    Construct `UnitsTable` from `ParticlesSetup`

    ## Input 
        - ParticlesSetup setup: The parameter of setup 

    ## Output
        - `UnitsTable`: The code units container.
    */
    static UnitsTable setup_units(const ParticlesSetup& setup){
        UnitsTable units = UnitsTable(setup.udist, setup.umass);
        return units;
    };
private:
protected:    
};
