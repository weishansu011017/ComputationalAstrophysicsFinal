#pragma once
#include <hdf5.h>

/*
Table of unit for setup
*/
class UnitsTable {
public:
    float udist;
    float utime;
    float umass;

    // Constructor
    UnitsTable() : udist(1.0f), utime(1.0f), umass(1.0f) {}  // 可選預設值

    // Destructor
    ~UnitsTable() = default;

    /*
        void read_unit_HDF5(hid_t file_id);
    Read the code unit setup in HDF5

    ## Input
        - hid_t file_id: The HDF5 file handle created by H5Fcreate or H5Fopen
    */
    void read_unit_HDF5(hid_t file_id);
private:
protected:    
};
