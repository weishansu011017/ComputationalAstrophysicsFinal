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

void ParticlesTable::calculate_a_dirnbody(){
    // NEED IMPLEMENT
}

void ParticlesTable::calculate_a_BHtree(){
    // NEED IMPLEMENT
}

void ParticlesTable::kick(float scale){
    // NEED IMPLEMENT
}

void ParticlesTable::drift(float scale){
    // NEED IMPLEMENT
}

// Constructor implementations
ParticlesTable::QuadTreeNode::QuadTreeNode(float xmin, float xmax, float ymin, float ymax) 
    : totalMass(0.0f) {
    bounds[0] = xmin; bounds[1] = xmax; bounds[2] = ymin; bounds[3] = ymax;
    centerOfMass[0] = centerOfMass[1] = 0.0f;
    std::fill(children, children + 4, nullptr);
}

ParticlesTable::QuadTreeNode::~QuadTreeNode() {
    for (int i = 0; i < 4; i++) {
        delete children[i];
    }
}

ParticlesTable::OctTreeNode::OctTreeNode(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax)
    : totalMass(0.0f) {
    bounds[0] = xmin; bounds[1] = xmax; bounds[2] = ymin; 
    bounds[3] = ymax; bounds[4] = zmin; bounds[5] = zmax;
    centerOfMass[0] = centerOfMass[1] = centerOfMass[2] = 0.0f;
    std::fill(children, children + 8, nullptr);
}

ParticlesTable::OctTreeNode::~OctTreeNode() {
    for (int i = 0; i < 8; i++) {
        delete children[i];
    }
}

// Quadrant/Octant determination
int ParticlesTable::QuadTreeNode::getQuadrant(float x, float y) const {
    float midX = (bounds[0] + bounds[1]) / 2;
    float midY = (bounds[2] + bounds[3]) / 2;
    
    int quad = 0;
    if (x >= midX) quad |= 1;  // East
    if (y >= midY) quad |= 2;  // North
    return quad; // 0=SW, 1=SE, 2=NW, 3=NE
}

int ParticlesTable::OctTreeNode::getOctant(float x, float y, float z) const {
    float midX = (bounds[0] + bounds[1]) / 2;
    float midY = (bounds[2] + bounds[3]) / 2; 
    float midZ = (bounds[4] + bounds[5]) / 2;
    
    int oct = 0;
    if (x >= midX) oct |= 1;  // East
    if (y >= midY) oct |= 2;  // North  
    if (z >= midZ) oct |= 4;  // Up
    return oct;
}

// Main tree building methods
void ParticlesTable::buildQuadTree() {
    // Clear existing tree
    if (quadTreeRoot) {
        clearQuadTree(quadTreeRoot);
    }
    
    // Find bounding box of all particles
    float xmin = x[0], xmax = x[0], ymin = y[0], ymax = y[0];
    for (int i = 1; i < N; i++) {
        if (x[i] < xmin) xmin = x[i];
        if (x[i] > xmax) xmax = x[i];
        if (y[i] < ymin) ymin = y[i]; 
        if (y[i] > ymax) ymax = y[i];
    }
    
    // Add small padding
    float padding = 0.1f;
    xmin -= padding; xmax += padding;
    ymin -= padding; ymax += padding;
    
    // Create root node
    quadTreeRoot = new QuadTreeNode(xmin, xmax, ymin, ymax);
    
    // Insert all particles
    for (int i = 0; i < N; i++) {
        quadTreeRoot->insert(i, *this);
    }
    
    // Calculate centers of mass
    quadTreeRoot->calculateCenterOfMass(*this);
    useQuadTree = true;
}

void ParticlesTable::buildOctTree() {
    // Clear existing tree
    if (octTreeRoot) {
        clearOctTree(octTreeRoot);
    }
    
    // Find bounding box of all particles
    float xmin = x[0], xmax = x[0], ymin = y[0], ymax = y[0], zmin = z[0], zmax = z[0];
    for (int i = 1; i < N; i++) {
        if (x[i] < xmin) xmin = x[i];
        if (x[i] > xmax) xmax = x[i];
        if (y[i] < ymin) ymin = y[i];
        if (y[i] > ymax) ymax = y[i];
        if (z[i] < zmin) zmin = z[i];
        if (z[i] > zmax) zmax = z[i];
    }
    
    // Add small padding
    float padding = 0.1f;
    xmin -= padding; xmax += padding;
    ymin -= padding; ymax += padding;
    zmin -= padding; zmax += padding;
    
    // Create root node
    octTreeRoot = new OctTreeNode(xmin, xmax, ymin, ymax, zmin, zmax);
    
    // Insert all particles
    for (int i = 0; i < N; i++) {
        octTreeRoot->insert(i, *this);
    }
    
    // Calculate centers of mass
    octTreeRoot->calculateCenterOfMass(*this);
    useQuadTree = false;
}

// Insert methods (simplified - just add to appropriate leaf)
void ParticlesTable::QuadTreeNode::insert(int particleIndex, ParticlesTable& particles) {
    if (isLeaf()) {
        particleIndices.push_back(particleIndex);
        
        // If this leaf now has more than 1 particle, subdivide
        if (particleIndices.size() > 1) {
            // Create children
            float midX = (bounds[0] + bounds[1]) / 2;
            float midY = (bounds[2] + bounds[3]) / 2;
            
            children[0] = new QuadTreeNode(bounds[0], midX, bounds[2], midY); // SW
            children[1] = new QuadTreeNode(midX, bounds[1], bounds[2], midY); // SE
            children[2] = new QuadTreeNode(bounds[0], midX, midY, bounds[3]); // NW  
            children[3] = new QuadTreeNode(midX, bounds[1], midY, bounds[3]); // NE
            
            // Redistribute particles to children
            for (int idx : particleIndices) {
                int quad = getQuadrant(particles.x[idx], particles.y[idx]);
                children[quad]->particleIndices.push_back(idx);
            }
            particleIndices.clear();
        }
    } else {
        // Internal node - pass to appropriate child
        int quad = getQuadrant(particles.x[particleIndex], particles.y[particleIndex]);
        children[quad]->insert(particleIndex, particles);
    }
}

void ParticlesTable::OctTreeNode::insert(int particleIndex, ParticlesTable& particles) {
    if (isLeaf()) {
        particleIndices.push_back(particleIndex);
        
        // If this leaf now has more than 1 particle, subdivide
        if (particleIndices.size() > 1) {
            // Create children
            float midX = (bounds[0] + bounds[1]) / 2;
            float midY = (bounds[2] + bounds[3]) / 2;
            float midZ = (bounds[4] + bounds[5]) / 2;
            
            for (int i = 0; i < 8; i++) {
                float xmin = (i & 1) ? midX : bounds[0];
                float xmax = (i & 1) ? bounds[1] : midX;
                float ymin = (i & 2) ? midY : bounds[2]; 
                float ymax = (i & 2) ? bounds[3] : midY;
                float zmin = (i & 4) ? midZ : bounds[4];
                float zmax = (i & 4) ? bounds[5] : midZ;
                
                children[i] = new OctTreeNode(xmin, xmax, ymin, ymax, zmin, zmax);
            }
            
            // Redistribute particles to children
            for (int idx : particleIndices) {
                int oct = getOctant(particles.x[idx], particles.y[idx], particles.z[idx]);
                children[oct]->particleIndices.push_back(idx);
            }
            particleIndices.clear();
        }
    } else {
        // Internal node - pass to appropriate child
        int oct = getOctant(particles.x[particleIndex], particles.y[particleIndex], particles.z[particleIndex]);
        children[oct]->insert(particleIndex, particles);
    }
}

// Center of mass calculation
void ParticlesTable::QuadTreeNode::calculateCenterOfMass(ParticlesTable& particles) {
    if (isLeaf()) {
        totalMass = 0.0f;
        centerOfMass[0] = centerOfMass[1] = 0.0f;
        
        for (int idx : particleIndices) {
            float mass = particles.m[idx];
            centerOfMass[0] += particles.x[idx] * mass;
            centerOfMass[1] += particles.y[idx] * mass;
            totalMass += mass;
        }
        
        if (totalMass > 0) {
            centerOfMass[0] /= totalMass;
            centerOfMass[1] /= totalMass;
        }
    } else {
        totalMass = 0.0f;
        centerOfMass[0] = centerOfMass[1] = 0.0f;
        
        for (int i = 0; i < 4; i++) {
            if (children[i]) {
                children[i]->calculateCenterOfMass(particles);
                float mass = children[i]->totalMass;
                centerOfMass[0] += children[i]->centerOfMass[0] * mass;
                centerOfMass[1] += children[i]->centerOfMass[1] * mass;
                totalMass += mass;
            }
        }
        
        if (totalMass > 0) {
            centerOfMass[0] /= totalMass;
            centerOfMass[1] /= totalMass;
        }
    }
}

void ParticlesTable::OctTreeNode::calculateCenterOfMass(ParticlesTable& particles) {
    if (isLeaf()) {
        totalMass = 0.0f;
        centerOfMass[0] = centerOfMass[1] = centerOfMass[2] = 0.0f;
        
        for (int idx : particleIndices) {
            float mass = particles.m[idx];
            centerOfMass[0] += particles.x[idx] * mass;
            centerOfMass[1] += particles.y[idx] * mass;
            centerOfMass[2] += particles.z[idx] * mass;
            totalMass += mass;
        }
        
        if (totalMass > 0) {
            centerOfMass[0] /= totalMass;
            centerOfMass[1] /= totalMass;
            centerOfMass[2] /= totalMass;
        }
    } else {
        totalMass = 0.0f;
        centerOfMass[0] = centerOfMass[1] = centerOfMass[2] = 0.0f;
        
        for (int i = 0; i < 8; i++) {
            if (children[i]) {
                children[i]->calculateCenterOfMass(particles);
                float mass = children[i]->totalMass;
                centerOfMass[0] += children[i]->centerOfMass[0] * mass;
                centerOfMass[1] += children[i]->centerOfMass[1] * mass;
                centerOfMass[2] += children[i]->centerOfMass[2] * mass;
                totalMass += mass;
            }
        }
        
        if (totalMass > 0) {
            centerOfMass[0] /= totalMass;
            centerOfMass[1] /= totalMass;
            centerOfMass[2] /= totalMass;
        }
    }
}

void ParticlesTable::setBHTreeTheta(float theta) {
    bhTreeTheta = theta;
}

// Tree saving methods
void ParticlesTable::saveQuadTreeToFile(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open " << filename << " for writing\n";
        return;
    }
    
    if (quadTreeRoot) {
        file << "QUADTREE\n";
        file << "# Format: depth xmin xmax ymin ymax centerX centerY totalMass numParticles [particleIndices...]\n";
        saveQuadTreeNode(file, quadTreeRoot, 0);
    } else {
        std::cerr << "No quadtree built. Call buildQuadTree() first.\n";
    }
    
    file.close();
}

void ParticlesTable::saveOctTreeToFile(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open " << filename << " for writing\n";
        return;
    }
    
    if (octTreeRoot) {
        file << "OCTREE\n";
        file << "# Format: depth xmin xmax ymin ymax zmin zmax centerX centerY centerZ totalMass numParticles [particleIndices...]\n";
        saveOctTreeNode(file, octTreeRoot, 0);
    } else {
        std::cerr << "No octree built. Call buildOctTree() first.\n";
    }
    
    file.close();
}

void ParticlesTable::saveQuadTreeNode(std::ofstream& file, QuadTreeNode* node, int depth) const {
    if (!node) return;
    
    file << depth << " "
         << node->bounds[0] << " " << node->bounds[1] << " "  // xmin xmax
         << node->bounds[2] << " " << node->bounds[3] << " "  // ymin ymax
         << node->centerOfMass[0] << " " << node->centerOfMass[1] << " "  // center
         << node->totalMass << " "
         << node->particleIndices.size();
    
    for (int idx : node->particleIndices) {
        file << " " << idx;
    }
    file << "\n";
    
    // Recursively save children
    for (int i = 0; i < 4; i++) {
        if (node->children[i]) {
            saveQuadTreeNode(file, node->children[i], depth + 1);
        }
    }
}

void ParticlesTable::saveOctTreeNode(std::ofstream& file, OctTreeNode* node, int depth) const {
    if (!node) return;
    
    file << depth << " "
         << node->bounds[0] << " " << node->bounds[1] << " "  // xmin xmax
         << node->bounds[2] << " " << node->bounds[3] << " "  // ymin ymax
         << node->bounds[4] << " " << node->bounds[5] << " "  // zmin zmax
         << node->centerOfMass[0] << " " << node->centerOfMass[1] << " " << node->centerOfMass[2] << " "  // center
         << node->totalMass << " "
         << node->particleIndices.size();
    
    for (int idx : node->particleIndices) {
        file << " " << idx;
    }
    file << "\n";
    
    // Recursively save children
    for (int i = 0; i < 8; i++) {
        if (node->children[i]) {
            saveOctTreeNode(file, node->children[i], depth + 1);
        }
    }
}

// Cleanup methods
void ParticlesTable::clearQuadTree(QuadTreeNode* node) {
    delete node;  // Destructor handles children recursively
}

void ParticlesTable::clearOctTree(OctTreeNode* node) {
    delete node;  // Destructor handles children recursively
}