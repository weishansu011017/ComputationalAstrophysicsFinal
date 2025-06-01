#include "OctTree.hpp"
#include "ParticlesTable.hpp"
#include <iostream>
#include <algorithm>
#include <cstdlib>

// Node implementation
OctTree::Node::Node(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax)
    : totalMass(0.0f) {
    bounds[0] = xmin; bounds[1] = xmax; bounds[2] = ymin; 
    bounds[3] = ymax; bounds[4] = zmin; bounds[5] = zmax;
    centerOfMass[0] = centerOfMass[1] = centerOfMass[2] = 0.0f;
    std::fill(children, children + 8, nullptr);
}

OctTree::Node::~Node() {
    for (int i = 0; i < 8; i++) {
        delete children[i];
    }
}

int OctTree::Node::getOctant(float x, float y, float z) const {
    float midX = (bounds[0] + bounds[1]) / 2;
    float midY = (bounds[2] + bounds[3]) / 2; 
    float midZ = (bounds[4] + bounds[5]) / 2;
    
    int oct = 0;
    if (x >= midX) oct |= 1;  // East
    if (y >= midY) oct |= 2;  // North  
    if (z >= midZ) oct |= 4;  // Up
    return oct;
}

void OctTree::Node::insert(int particleIndex, const ParticlesTable& particles) {
    if (isLeaf()) {
        particleIndices.push_back(particleIndex);
        
        // If this leaf now has more than 1 particle, subdivide
        if (particleIndices.size() > 1) {
            // Check if we can still subdivide (avoid infinite recursion for identical positions)
            float width = bounds[1] - bounds[0];
            float height = bounds[3] - bounds[2];
            float depth = bounds[5] - bounds[4];
            const float MIN_SIZE = 1e-8f;  // Minimum cell size to prevent infinite subdivision
            
            if (width > MIN_SIZE && height > MIN_SIZE && depth > MIN_SIZE) {
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
                    
                    children[i] = new Node(xmin, xmax, ymin, ymax, zmin, zmax);
                }
                
                // Redistribute particles to children using insert (for recursive subdivision)
                std::vector<int> tempIndices = particleIndices;
                particleIndices.clear();
                
                for (int idx : tempIndices) {
                    int oct = getOctant(particles.x[idx], particles.y[idx], particles.z[idx]);
                    children[oct]->insert(idx, particles);
                }
            }
            // If cell is too small, keep particles in this leaf
        }
    } else {
        // Internal node - pass to appropriate child
        int oct = getOctant(particles.x[particleIndex], particles.y[particleIndex], particles.z[particleIndex]);
        children[oct]->insert(particleIndex, particles);
    }
}

void OctTree::Node::calculateCenterOfMass(const ParticlesTable& particles) {
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

// OctTree implementation
OctTree::OctTree() : root(nullptr), theta(0.5f) {}

OctTree::OctTree(float openingAngle) : root(nullptr), theta(openingAngle) {}

OctTree::~OctTree() {
    clear();
}

void OctTree::buildFromParticles(const ParticlesTable& particles) {
    // Clear existing tree
    clear();
    
    if (particles.N == 0) return;
    
    // Find bounding box of all particles
    float xmin = particles.x[0], xmax = particles.x[0];
    float ymin = particles.y[0], ymax = particles.y[0];
    float zmin = particles.z[0], zmax = particles.z[0];
    
    for (int i = 1; i < particles.N; i++) {
        if (particles.x[i] < xmin) xmin = particles.x[i];
        if (particles.x[i] > xmax) xmax = particles.x[i];
        if (particles.y[i] < ymin) ymin = particles.y[i];
        if (particles.y[i] > ymax) ymax = particles.y[i];
        if (particles.z[i] < zmin) zmin = particles.z[i];
        if (particles.z[i] > zmax) zmax = particles.z[i];
    }
    
    // Add small padding
    float padding = 0.1f;
    xmin -= padding; xmax += padding;
    ymin -= padding; ymax += padding;
    zmin -= padding; zmax += padding;
    
    // Create root node
    root = new Node(xmin, xmax, ymin, ymax, zmin, zmax);
    
    // Insert all particles
    for (int i = 0; i < particles.N; i++) {
        root->insert(i, particles);
    }
    
    // Calculate centers of mass
    root->calculateCenterOfMass(particles);
}

OctTree OctTree::buildOctTree(const ParticlesTable& particles) {
    OctTree tree;
    tree.buildFromParticles(particles);
    return tree;
}

void OctTree::clear() {
    if (root) {
        clearNode(root);
        root = nullptr;
    }
}

void OctTree::clearNode(Node* node) {
    delete node;  // Destructor handles children recursively
}

void OctTree::saveToFile(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open " << filename << " for writing\n";
        std::exit(1);
    }
    
    file << "OCTREE\n";
    file << "# Format: depth xmin xmax ymin ymax zmin zmax centerX centerY centerZ totalMass numParticles [particleIndices...]\n";
    
    if (root) {
        saveNode(file, root, 0);
    }
    
    file.close();
}

void OctTree::saveNode(std::ofstream& file, Node* node, int depth) const {
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
            saveNode(file, node->children[i], depth + 1);
        }
    }
}