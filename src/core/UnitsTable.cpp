#include <hdf5.h>
#include <string>
#include "UnitsTable.hpp"
#include "ParticlesSetup.hpp"

UnitsTable UnitsTable::setup_units(const ParticlesSetup& setup){
        UnitsTable units = UnitsTable(setup.udist, setup.umass);
        return units;
    };

void UnitsTable::read_unit_HDF5(hid_t file_id){
    auto read_scalar = [&](const char* group, const char* name, hid_t h5type, void* dst) {
        std::string path = std::string(group) + "/" + name;
        hid_t dset = H5Dopen2(file_id, path.c_str(), H5P_DEFAULT);
        H5Dread(dset, h5type, H5S_ALL, H5S_ALL, H5P_DEFAULT, dst);
        H5Dclose(dset);
    };

    // Read UnitsTable
    read_scalar("/params", "udist", H5T_NATIVE_FLOAT, &udist);
    read_scalar("/params", "utime", H5T_NATIVE_FLOAT, &utime);
    read_scalar("/params", "umass", H5T_NATIVE_FLOAT, &umass); 
}