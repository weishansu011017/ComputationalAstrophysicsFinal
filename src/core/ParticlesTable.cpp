#include <stdlib.h>
#include <hdf5.h>
#include <iostream>
#include <omp.h>
#include "UnitsTable.hpp"
#include "ParticlesTable.hpp"
#include "QuadTree.hpp"
#include "OctTree.hpp"

// Forward declarations for helper functions (implementation details)
void calculateForceFromQuadNode(ParticlesTable* pt, int particleIndex, QuadTree::Node* node, float theta);
void calculateForceFromOctNode(ParticlesTable* pt, int particleIndex, OctTree::Node* node, float theta);

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

void ParticlesTable::calculate_h(){
    // NEED IMPLEMENT
}

void ParticlesTable::calculate_dt(){
    // NEED IMPLEMENT
}

void ParticlesTable::calculate_a_dirnbody(){
    // Direct calculate acc
    #pragma omp parallel for 
    for (int i = 0; i < N; ++i){
        float axtemp = 0.0;
        float aytemp = 0.0;
        float aztemp = 0.0;
        for (int j = 0; j < N; ++j){
            if (i == j) continue;
            float dx = x[j] - x[i];
            float dy = y[j] - y[i];
            float dz = z[j] - z[i];
            float dr2 = dx * dx + dy * dy + dz * dz + h[i] * h[i];
            float invr = 1.0f / sqrtf(dr2);
            float invr3 = invr * invr * invr;
            float mjinvr3 = m[j] * invr3;
            axtemp += mjinvr3 * dx;
            aytemp += mjinvr3 * dy;
            aztemp += mjinvr3 * dz;
        }
        _ax[i] = axtemp;
        _ay[i] = aytemp;
        _az[i] = aztemp;
    }
}

// Helper function to determine if particles are essentially 2D
bool ParticlesTable::is2D(float tolerance) const {
    float zmin = z[0], zmax = z[0];
    for (int i = 1; i < N; i++) {
        if (z[i] < zmin) zmin = z[i];
        if (z[i] > zmax) zmax = z[i];
    }
    float zrange = zmax - zmin;
    
    // Get x,y range for comparison
    float xmin = x[0], xmax = x[0];
    float ymin = y[0], ymax = y[0];
    for (int i = 1; i < N; i++) {
        if (x[i] < xmin) xmin = x[i];
        if (x[i] > xmax) xmax = x[i];
        if (y[i] < ymin) ymin = y[i];
        if (y[i] > ymax) ymax = y[i];
    }
    float xrange = xmax - xmin;
    float yrange = ymax - ymin;
    float xyrange = std::max(xrange, yrange);
    
    // Consider 2D if z-range is less than tolerance times xy-range
    return (zrange < tolerance * xyrange);
}

void ParticlesTable::calculate_a_BHtree(){
    // Check if essentially 2D
    if (is2D(0.01f)) {
        calculate_a_BHtree_2D();
        return;
    }
    
    // Build octree
    OctTree tree(bhTreeTheta);
    tree.buildFromParticles(*this);
    
    // Calculate forces for each particle
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        _ax[i] = 0.0f;
        _ay[i] = 0.0f;
        _az[i] = 0.0f;
        
        // Tree walk for particle i
        calculateForceFromOctNode(this, i, tree.getRoot(), bhTreeTheta);
    }
}

void ParticlesTable::calculate_a_BHtree_2D(){
    // Build quadtree
    QuadTree tree(bhTreeTheta);
    tree.buildFromParticles(*this);
    
    // Calculate forces for each particle
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        _ax[i] = 0.0f;
        _ay[i] = 0.0f;
        _az[i] = 0.0f;
        
        // Tree walk for particle i
        calculateForceFromQuadNode(this, i, tree.getRoot(), bhTreeTheta);
    }
}

// Calculate force from a quadtree node
void calculateForceFromQuadNode(ParticlesTable* pt, int particleIndex, QuadTree::Node* node, float theta) {
    // Step 1: Safety checks
    if (!node || node->totalMass == 0) return;
    
    // Step 2: Calculate displacement vector
    // Get particle position
    float px = pt->x[particleIndex];
    float py = pt->y[particleIndex];
    
    // Calculate distance to node's center of mass
    float dx = node->centerOfMass[0] - px;  // Vector from particle to COM
    float dy = node->centerOfMass[1] - py;
    float r2 = dx * dx + dy * dy;           // Distance squared
    
    if (r2 < 1e-8f) return;  // Skip if too close
    
    // Calculate cell size
    float cellSize = std::max(node->bounds[1] - node->bounds[0], 
                              node->bounds[3] - node->bounds[2]);
    
    // Opening angle criterion
    float distance = std::sqrt(r2);
    
    // Step 3: Three-way decision tree
    if (node->isLeaf()) {
        // Leaf node - calculate forces from all particles in this leaf
        for (int j : node->particleIndices) {
            if (j == particleIndex) continue;  // Skip self
            
            // Newton's law of gravitation: F = G*m1*m2/r²
            // In code units where G=1: F = m1*m2/r²
            float dxj = pt->x[j] - px;
            float dyj = pt->y[j] - py;
            float r2j = dxj * dxj + dyj * dyj + pt->h[particleIndex] * pt->h[particleIndex];    // Add softening
            float invr = 1.0f / std::sqrt(r2j);
            float invr3 = invr * invr * invr;
            float mjinvr3 = pt->m[j] * invr3;

            // Acceleration: a = F/m = G*m_other/r² * direction
            pt->_ax[particleIndex] += mjinvr3 * dxj;
            pt->_ay[particleIndex] += mjinvr3 * dyj;
        }
    } else if (cellSize / distance < theta) {
        // Node is far enough - use center of mass approximation
        r2 += pt->h[particleIndex] * pt->h[particleIndex];  // Add softening
        float invr = 1.0f / std::sqrt(r2);
        float invr3 = invr * invr * invr;
        float minvr3 = node->totalMass * invr3;
        
        // Treat entire node as single particle at center of mass
        pt->_ax[particleIndex] += minvr3 * dx;
        pt->_ay[particleIndex] += minvr3 * dy;
    } else {
        // Node is too close - recurse into children
        for (int i = 0; i < 4; i++) {
            calculateForceFromQuadNode(pt, particleIndex, node->children[i], theta);
        }
    }
}

// Calculate force from an octree node
void calculateForceFromOctNode(ParticlesTable* pt, int particleIndex, OctTree::Node* node, float theta) {
    // Same logic as QuadTree but with 3D coordinates
    if (!node || node->totalMass == 0) return;
    
    // Get particle position
    float px = pt->x[particleIndex];
    float py = pt->y[particleIndex];
    float pz = pt->z[particleIndex];
    
    // Calculate distance to node's center of mass
    float dx = node->centerOfMass[0] - px;
    float dy = node->centerOfMass[1] - py;
    float dz = node->centerOfMass[2] - pz;
    float r2 = dx * dx + dy * dy + dz * dz;
    
    if (r2 < 1e-8f) return;  // Skip if too close
    
    // Calculate cell size
    float cellSize = std::max({node->bounds[1] - node->bounds[0],
                               node->bounds[3] - node->bounds[2],
                               node->bounds[5] - node->bounds[4]});
    
    // Opening angle criterion
    float distance = std::sqrt(r2);
    
    if (node->isLeaf()) {
        // Leaf node - calculate forces from all particles in this leaf
        for (int j : node->particleIndices) {
            if (j == particleIndex) continue;  // Skip self
            
            float dxj = pt->x[j] - px;
            float dyj = pt->y[j] - py;
            float dzj = pt->z[j] - pz;
            float r2j = dxj * dxj + dyj * dyj + dzj * dzj + pt->h[particleIndex] * pt->h[particleIndex];
            float invr = 1.0f / std::sqrt(r2j);
            float invr3 = invr * invr * invr;
            float mjinvr3 = pt->m[j] * invr3;
            
            pt->_ax[particleIndex] += mjinvr3 * dxj;
            pt->_ay[particleIndex] += mjinvr3 * dyj;
            pt->_az[particleIndex] += mjinvr3 * dzj;
        }
    } else if (cellSize / distance < theta) {
        // Node is far enough - use center of mass approximation
        r2 += pt->h[particleIndex] * pt->h[particleIndex];  // Add softening
        float invr = 1.0f / std::sqrt(r2);
        float invr3 = invr * invr * invr;
        float minvr3 = node->totalMass * invr3;
        
        pt->_ax[particleIndex] += minvr3 * dx;
        pt->_ay[particleIndex] += minvr3 * dy;
        pt->_az[particleIndex] += minvr3 * dz;
    } else {
        // Node is too close - recurse into children
        for (int i = 0; i < 8; i++) {
            calculateForceFromOctNode(pt, particleIndex, node->children[i], theta);
        }
    }
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

// Optional: Add method to set theta parameter
void ParticlesTable::setBHTreeTheta(float theta) {
    bhTreeTheta = theta;
}

// Optional: Add method to validate particles
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
    return QuadTree::buildQuadTree(*this);
}

OctTree ParticlesTable::buildOctTree() const {
    return OctTree::buildOctTree(*this);
}