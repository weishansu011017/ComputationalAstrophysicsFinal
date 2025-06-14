#pragma once
#include <cmath>
#include <iostream>
#include <string>
#include <cstdlib>
#include "PhysicalConstants.hpp"

class ParticlesSetup;
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
        if (_udist <= 0.0f || _umass <= 0.0f){
            std::cerr << "Code units must be positive. Received udist = " << _udist << ", umass = " << std::to_string(_umass) << std::endl;
            std::exit(1);
        }
        double dudist = _udist;
        double dumass = _umass;
        double dutime = utime;
        dutime = std::sqrt((dudist * dudist * dudist) / (dumass * CGSConstants::G));
        utime = dutime;
    }   

    // Destructor
    ~UnitsTable() = default;

    /*
        void read_unit_HDF5(hid_t file_id);
    Read the code unit setup in HDF5

    ## Input
        - hid_t file_id: The HDF5 file handle created by H5Fcreate or H5Fopen
    */
    void read_unit_HDF5(int64_t file_id);

    /*
        UnitsTable setup_units(const ParticlesSetup& setup);
    Construct `UnitsTable` from `ParticlesSetup`

    ## Input 
        - ParticlesSetup setup: The parameter of setup 

    ## Output
        - `UnitsTable`: The code units container.
    */
    static UnitsTable setup_units(const ParticlesSetup& setup);
private:
protected:    
};
