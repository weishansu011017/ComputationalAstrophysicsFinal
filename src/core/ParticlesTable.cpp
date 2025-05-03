#include <stdlib.h>
#include <hdf5.h>
#include <iostream>
#include "UnitsTable.hpp"
#include "ParticlesTable.hpp"

void ParticlesTable::_write_base_HDF5(hid_t file_id) const {
    // ============== /params/ ==============
    hid_t g_params = H5Gcreate2(file_id, "/params", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // Lambda function for writing parameter
    auto write_scalar = [&](const char* name, auto val, hid_t h5type) {
        hid_t dataspace_id = H5Screate(H5S_SCALAR);
        hid_t dset_id = H5Dcreate2(g_params, name, h5type, dataspace_id,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset_id, h5type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &val);
        H5Dclose(dset_id);
        H5Sclose(dataspace_id);
    };
    auto write_string = [&](const char* name, const std::string& str) {
        hid_t strtype = H5Tcopy(H5T_C_S1);
        H5Tset_size(strtype, str.size() + 1);  // +1 for null terminator
        H5Tset_strpad(strtype, H5T_STR_NULLTERM);
    
        hid_t dataspace_id = H5Screate(H5S_SCALAR);
        hid_t dset_id = H5Dcreate2(g_params, name, strtype, dataspace_id,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
        H5Dwrite(dset_id, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, str.c_str());
    
        H5Dclose(dset_id);
        H5Sclose(dataspace_id);
        H5Tclose(strtype);
    };
    // Write UnitsTable
    write_scalar("udist", unittable.udist, H5T_NATIVE_FLOAT);
    write_scalar("utime", unittable.utime, H5T_NATIVE_FLOAT);
    write_scalar("umass", unittable.umass, H5T_NATIVE_FLOAT);
    // Write Other parameter
    write_scalar("N", N, H5T_NATIVE_INT);
    write_scalar("t", t, H5T_NATIVE_FLOAT);
    write_scalar("Mtot", Mtot, H5T_NATIVE_FLOAT);
    write_string("SimulationTag", SimulationTag);

    // Close Group
    H5Gclose(g_params);

    // ============== /ParticlesTable/ ==============
    hid_t g_table = H5Gcreate2(file_id, "/Table", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    auto write_float_vector = [&](const char* name, const std::vector<float>& vec) {
        hsize_t dims[1] = { vec.size() };
        hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
        hid_t dset_id = H5Dcreate2(g_table, name, H5T_NATIVE_FLOAT, dataspace_id,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec.data());
        H5Dclose(dset_id);
        H5Sclose(dataspace_id);
    };
    auto write_uint32_vector = [&](const char* name, const std::vector<uint32_t>& vec) {
        hsize_t dims[1] = { vec.size() };
        hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
        hid_t dset_id = H5Dcreate2(g_table, name, H5T_NATIVE_UINT32, dataspace_id,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec.data());
        H5Dclose(dset_id);
        H5Sclose(dataspace_id);
    };
    // Write particle attribute
    write_uint32_vector("particle_index", particle_index);
    write_float_vector("x", x);
    write_float_vector("y", y);
    write_float_vector("z", z);
    write_float_vector("vx", vx);
    write_float_vector("vy", vy);
    write_float_vector("vz", vz);
    write_float_vector("h", h);
    write_float_vector("m", m);
    write_float_vector("dt", dt);
    write_float_vector("a", a);

    // Close Group
    H5Gclose(g_table);
}

void ParticlesTable::_read_base_HDF5(hid_t file_id){
    // ============== /params/ ==============
    auto read_scalar = [&](const char* group, const char* name, hid_t h5type, void* dst) {
        std::string path = std::string(group) + "/" + name;
        hid_t dset = H5Dopen2(file_id, path.c_str(), H5P_DEFAULT);
        H5Dread(dset, h5type, H5S_ALL, H5S_ALL, H5P_DEFAULT, dst);
        H5Dclose(dset);
    };
    auto read_string = [&](const char* group, const char* name, std::string& dst) {
        std::string path = std::string(group) + "/" + name;
        hid_t dset = H5Dopen2(file_id, path.c_str(), H5P_DEFAULT);
        hid_t type = H5Dget_type(dset);
    
        if (H5Tis_variable_str(type)) {
            char* rdata = nullptr;
            H5Dread(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rdata);
            dst = std::string(rdata);
            free(rdata); 
        } else {
            size_t size = H5Tget_size(type);
            std::vector<char> buf(size + 1, '\0'); 
            H5Dread(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
            dst = std::string(buf.data());
        }
    
        H5Tclose(type);
        H5Dclose(dset);
    };

    // [Optional] Read Other parameters
    read_scalar("/params","t", H5T_NATIVE_FLOAT, &t);
    read_scalar("/params","Mtot", H5T_NATIVE_FLOAT, &Mtot);
    read_string("/params", "SimulationTag", SimulationTag);

    // ============== /ParticlesTable/ ==============
    auto read_float_vector = [&](const char* group, const char* name, std::vector<float>& vec) {
        std::string path = std::string(group) + "/" + name;
        hid_t dset = H5Dopen2(file_id, path.c_str(), H5P_DEFAULT);
        hid_t space = H5Dget_space(dset);
        hsize_t dims[1]; H5Sget_simple_extent_dims(space, dims, NULL);
        H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec.data());
        H5Sclose(space);
        H5Dclose(dset);
    };

    auto read_uint32_vector = [&](const char* group, const char* name, std::vector<uint32_t>& vec) {
        std::string path = std::string(group) + "/" + name;
        hid_t dset = H5Dopen2(file_id, path.c_str(), H5P_DEFAULT);
        hid_t space = H5Dget_space(dset);
        hsize_t dims[1]; H5Sget_simple_extent_dims(space, dims, NULL);
        H5Dread(dset, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec.data());
        H5Sclose(space);
        H5Dclose(dset);
    };

    // Read table data
    read_uint32_vector("/Table", "particle_index", particle_index);
    read_float_vector("/Table", "x", x);
    read_float_vector("/Table", "y", y);
    read_float_vector("/Table", "z", z);
    read_float_vector("/Table", "vx", vx);
    read_float_vector("/Table", "vy", vy);
    read_float_vector("/Table", "vz", vz);
    read_float_vector("/Table","h", h);
    read_float_vector("/Table","m", m);
    read_float_vector("/Table","dt", dt);
    read_float_vector("/Table","a", a);
}


void ParticlesTable::extract_particles_table(const std::string& filename) const {
    // ============== Create Empty HDF5 File ==============
    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        std::cerr << "Failed to create HDF5 file.\n";
        return;
    }

    // ============== Write base properties ==============
    _write_base_HDF5(file_id);

    // ============== [Optional] Subclass: add additional datasets here ==============

    // ==============================================
    // Close File
    H5Fclose(file_id);
}


ParticlesTable ParticlesTable::read_particles_table(const std::string& filename){
    // ============== Open HDF5 File ==============
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        throw std::runtime_error("Failed to open HDF5 file for reading.");
    }

    // ============== Read base properties ==============
    // Read N for initilization table
    int N;
    hid_t dset = H5Dopen2(file_id, "/params/N", H5P_DEFAULT);
    if (dset < 0) {
        H5Fclose(file_id);
        throw std::runtime_error("Missing /params/N in HDF5 file");
    }
    H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &N);
    H5Dclose(dset);
    // Read code unit
    UnitsTable unit;
    unit.read_unit_HDF5(file_id);
    // Initialized container
    ParticlesTable pt(unit, N);
    // Read base table
    pt._read_base_HDF5(file_id);

    // ============== [Optional] Subclass: add additional datasets here ==============

    // ==============================================
    // Close File
    H5Fclose(file_id);

    return pt;
}