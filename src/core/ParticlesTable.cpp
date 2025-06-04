#include <stdlib.h>
#include <hdf5.h>
#include <iostream>
#include <omp.h>
#include <cstdlib>
#include "UnitsTable.hpp"
#include "ParticlesTable.hpp"
#include "QuadTree.hpp"
#include "OctTree.hpp"

void ParticlesTable::_write_base_HDF5(hid_t file_id, bool debug) const {
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
    write_scalar("dimension", dimension, H5T_NATIVE_INT);
    write_scalar("Mtot", Mtot, H5T_NATIVE_FLOAT);
    write_scalar("Utot", Utot, H5T_NATIVE_FLOAT);
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

    if (debug){
        write_float_vector("_ax", _ax);
        write_float_vector("_ay", _ay);
        write_float_vector("_az", _az);
        write_float_vector("_U", _U);
    }

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
    read_scalar("/params","Utot", H5T_NATIVE_FLOAT, &Utot);
    read_string("/params", "SimulationTag", SimulationTag);
    read_scalar("/params", "dimension", H5T_NATIVE_INT, &dimension);

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
}

void ParticlesTable::extract_particles_table(const std::string& filename, bool debug) const {
    // ============== Create Empty HDF5 File ==============
    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        std::cerr << "Failed to create HDF5 file." << std::endl;
        std::exit(1);
    }

    // ============== Write base properties ==============
    _write_base_HDF5(file_id, debug);

    // ============== [Optional] Subclass: add additional datasets here ==============

    // ==============================================
    // Close File
    H5Fclose(file_id);
}

ParticlesTable ParticlesTable::read_particles_table(const std::string& filename){
    // ============== Open HDF5 File ==============
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        std::cerr << "Failed to open HDF5 file for reading." << std::endl;
        std::exit(1);
    }

    // ============== Read base properties ==============
    // Read N for initilization table
    int N;
    hid_t dset = H5Dopen2(file_id, "/params/N", H5P_DEFAULT);
    if (dset < 0) {
        H5Fclose(file_id);
        std::cerr << "Missing /params/N in HDF5 file";
        std::exit(1);
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

void ParticlesTable::calculate_h(){
    // NEED IMPLEMENT
}

void ParticlesTable::calculate_dt(){
    // NEED IMPLEMENT
}

void ParticlesTable::calculate_Utot(){
    float Utemp = 0.0f;
    for (int i = 0; i < N; ++i){
        Utemp += _U[i];    
    }
    Utot = 0.5f * Utemp;
}


void ParticlesTable::calculate_a_dirnbody(){
    // Direct calculate acc
    #pragma omp parallel for 
    for (int i = 0; i < N; ++i){
        float axtemp = 0.0;
        float aytemp = 0.0;
        float aztemp = 0.0;
        float Utemp = 0.0;
        for (int j = 0; j < N; ++j){
            if (i == j) continue;
            float dx = x[j] - x[i];
            float dy = y[j] - y[i];
            float dz = z[j] - z[i];
            float dr2 = dx * dx + dy * dy + dz * dz + 0.5 * (h[i] * h[i] + h[j] * h[j]);
            float invr = 1.0f / sqrtf(dr2);
            float invr3 = invr * invr * invr;
            float mjinvr3 = m[j] * invr3;
            axtemp += mjinvr3 * dx;
            aytemp += mjinvr3 * dy;
            aztemp += mjinvr3 * dz;
            Utemp -= m[i] * m[j] * invr;
        }
        _ax[i] = axtemp;
        _ay[i] = aytemp;
        _az[i] = aztemp;
        _U[i] = Utemp;
    }
}

void ParticlesTable::calculate_a_dirnbody_2D(){
    // Direct calculate acc
    #pragma omp parallel for 
    for (int i = 0; i < N; ++i){
        float axtemp = 0.0;
        float aytemp = 0.0;
        float Utemp = 0.0;
        for (int j = 0; j < N; ++j){
            if (i == j) continue;
            float dx = x[j] - x[i];
            float dy = y[j] - y[i];
            float dr2 = dx * dx + dy * dy + 0.5 * (h[i] * h[i] + h[j] * h[j]);
            float invr = 1.0f / sqrtf(dr2);
            float invr3 = invr * invr * invr;
            float mjinvr3 = m[j] * invr3;
            axtemp += mjinvr3 * dx;
            aytemp += mjinvr3 * dy;
            Utemp -= m[i] * m[j] * invr;
        }
        _ax[i] = axtemp;
        _ay[i] = aytemp;
        _U[i] = Utemp;
    }
}

void ParticlesTable::calculate_a_BHtree(){
    // Build quadtree
    OctTree tree = buildOctTree();
    
    // initialize _ax, _ay
    std::fill(_ax.begin(), _ax.end(), 0.0f);
    std::fill(_ay.begin(), _ay.end(), 0.0f);
    std::fill(_az.begin(), _az.end(), 0.0f);
    std::fill(_U.begin(), _U.end(), 0.0f);

    // index of root
    const int root = tree.root_idx;

    // Calculate forces for each particle
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        _calculate_a_OctNode(i, tree, root);
    }
}

void ParticlesTable::calculate_a_BHtree_2D(){
    // Build quadtree
    QuadTree tree = buildQuadTree();
    
    // initialize _ax, _ay
    std::fill(_ax.begin(), _ax.end(), 0.0f);
    std::fill(_ay.begin(), _ay.end(), 0.0f);
    std::fill(_U.begin(), _U.end(), 0.0f);

    // index of root
    const int root = tree.root_idx;

    // Calculate forces for each particle
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        _calculate_a_QuadNode(i, tree, root);
    }
}

// Calculate force from an octree node
void ParticlesTable::_calculate_a_OctNode(int idx, const OctTree& tree, int nidx){
    // Assume x,y,m in tree has been well reordered
    const OctNode& node = tree.nodes_list[nidx];
    // If no particles in the node ==> return
    if (node.pcount == 0) return;

    // Current point xi, yi, hi
    float xi = x[idx];
    float yi = y[idx];
    float zi = z[idx];
    float mi = m[idx];
    float hi = h[idx];

    // Leaf node - calculate acc from all particles in this leaf by direct Nbody
    if (node.isLeaf()) {

        int count = node.pcount;
        int start = node.ParticlesLocateidx;
        for (int p = 0; p < count; ++p){                // For the particle in this leaf (count < leafNmax)
            int pidx = tree.order[start + p];           // get the original index from inverse_order
            if (pidx == idx) continue;                  // Comparing the original index to idx. Continue if they are the same (same particles)
            float dx = x[pidx] - xi;                    
            float dy = y[pidx] - yi;
            float dz = z[pidx] - zi;
            float r2 = dx*dx + dy*dy + dz*dz + 0.5 * (hi * hi + h[pidx] * h[pidx]);           // Distance between two particles + softerning
            float invr  = 1.0f / std::sqrt(r2);
            float invr3 = invr * invr * invr;
            float mjinvr3 = m[pidx] * invr3;

            _ax[idx] += mjinvr3 * dx;
            _ay[idx] += mjinvr3 * dy;
            _az[idx] += mjinvr3 * dz;
            _U[idx]  -= mi * m[pidx] * invr;
        }
        return;
    }
    
    // Not Leaf node: recursion or COM approx
    // Calculate distance to node's center of mass
    float dx = node.COMx - xi;
    float dy = node.COMy - yi;
    float dz = node.COMz - zi;
    float r2 = dx*dx + dy*dy + dz*dz;
    
    if (r2 > 1e-8f) {
        float dist  = std::sqrt(r2);

        if (std::max(std::max(node.width(), node.height()), node.depth()) / dist < tree.bhtheta) { 
            r2 += hi*hi;
            float invr  = 1.0f / std::sqrt(r2);
            float invr3 = invr * invr * invr;
            float mInvr3 = node.Mtot * invr3;
            float Mext = node.Mtot;

            _ax[idx] += mInvr3 * dx;
            _ay[idx] += mInvr3 * dy;
            _az[idx] += mInvr3 * dz;
            _U[idx]  -= mi * Mext * invr;
            return;
        }
    }
    // Other cases: Recursion (r <= 1e-8f or Node is too close)
    for (int q = 0; q < 8; ++q) {
        int cidx = node.children[q];
        if (cidx >= 0)
            _calculate_a_OctNode(idx, tree, cidx);
    }
    return;
}

// Calculate force from an quadtree node
void ParticlesTable::_calculate_a_QuadNode(int idx, const QuadTree& tree, int nidx){
    // Assume x,y,m in tree has been well reordered
    const QuadNode& node = tree.nodes_list[nidx];
    // If no particles in the node ==> return
    if (node.pcount == 0) return;

    // Current point xi, yi, hi
    float xi = x[idx];
    float yi = y[idx];
    float mi = m[idx];
    float hi = h[idx];

    // Leaf node - calculate acc from all particles in this leaf by direct Nbody
    if (node.isLeaf()) {

        int count = node.pcount;
        int start = node.ParticlesLocateidx;
        for (int p = 0; p < count; ++p){                // For the particle in this leaf (count < leafNmax)
            int pidx = tree.order[start + p];       // get the original index from inverse_order
            if (pidx == idx) continue;                  // Comparing the original index to idx. Continue if they are the same (same particles)
            float dx = x[pidx] - xi;                    
            float dy = y[pidx] - yi;
            float r2 = dx*dx + dy*dy + 0.5 * (hi * hi + h[pidx] * h[pidx]);;           // Distance between two particles + softerning
            float invr  = 1.0f / std::sqrt(r2);
            float invr3 = invr * invr * invr;
            float mjinvr3 = m[pidx] * invr3;

            _ax[idx] += mjinvr3 * dx;
            _ay[idx] += mjinvr3 * dy;
            _U[idx]  -= mi * m[pidx] * invr;
        }
        return;
    }
    
    // Not Leaf node: recursion or COM approx
    // Calculate distance to node's center of mass
    float dx = node.COMx - xi;
    float dy = node.COMy - yi;
    float r2 = dx*dx + dy*dy;
    
    if (r2 > 1e-8f) {
        float dist  = std::sqrt(r2);

        if (std::max(node.width(), node.height()) / dist < tree.bhtheta) { 
            r2 += hi*hi;
            float invr  = 1.0f / std::sqrt(r2);
            float invr3 = invr * invr * invr;
            float mInvr3 = node.Mtot * invr3;
            float Mext = node.Mtot;

            _ax[idx] += mInvr3 * dx;
            _ay[idx] += mInvr3 * dy;
            _U[idx]  -= mi * Mext * invr;
            return;
        }
    }
    // Other cases: Recursion (r <= 1e-8f or Node is too close)
    for (int q = 0; q < 4; ++q) {
        int cidx = node.children[q];
        if (cidx >= 0)
            _calculate_a_QuadNode(idx, tree, cidx);
    }
    return;
}

void ParticlesTable::kick(float scale) {
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        if (h[i] > 0.0f) { 
            vx[i] += scale * _ax[i] * dt[i];
            vy[i] += scale * _ay[i] * dt[i];
            vz[i] += scale * _az[i] * dt[i];
        }
    }
}

void ParticlesTable::kick_2D(float scale) {
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        if (h[i] > 0.0f) { 
            vx[i] += scale * _ax[i] * dt[i];
            vy[i] += scale * _ay[i] * dt[i];
        }
    }
}

void ParticlesTable::drift(float scale){
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        if (h[i] > 0.0f) { 
            x[i] += scale * vx[i] * dt[i];
            y[i] += scale * vy[i] * dt[i];
            z[i] += scale * vz[i] * dt[i];
        }
    }
}

void ParticlesTable::drift_2D(float scale){
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        if (h[i] > 0.0f) { 
            x[i] += scale * vx[i] * dt[i];
            y[i] += scale * vy[i] * dt[i];
        }
    }
}

void ParticlesTable::particles_validation(){
    // Check for NaN or infinite values
    for (int i = 0; i < N; i++) {
        if (!std::isfinite(x[i]) || !std::isfinite(y[i]) || !std::isfinite(z[i]) ||
            !std::isfinite(vx[i]) || !std::isfinite(vy[i]) || !std::isfinite(vz[i])) {
            std::cerr << "Warning: Particle " << i << " has invalid position/velocity!" << std::endl;
        }
    }
}

// Tree building methods - delegate to tree classes
QuadTree ParticlesTable::buildQuadTree() const {
    QuadTree quadtree(bhTreeTheta, N);
    quadtree.reserve_nodes(N);
    quadtree.build_tree(x, y, m);
    return quadtree;
}
OctTree ParticlesTable::buildOctTree() const {
    OctTree octtree(bhTreeTheta, N);
    octtree.reserve_nodes(N);
    octtree.build_tree(x, y, z, m);
    return octtree;
}

void ParticlesTable::set_bhTreeTheta(float theta){
    bhTreeTheta = theta;
}