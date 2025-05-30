#include "QuadTree.hpp"
#include "ParticlesTable.hpp"
#include <iostream>
#include <algorithm>

// Node implementation
QuadTree::Node::Node(float xmin, float xmax, float ymin, float ymax) 
    : totalMass(0.0f) {
    bounds[0] = xmin; bounds[1] = xmax; bounds[2] = ymin; bounds[3] = ymax;
    centerOfMass[0] = centerOfMass[1] = 0.0f;
    std::fill(children, children + 4, nullptr);
}

QuadTree::Node::~Node() {
    for (int i = 0; i < 4; i++) {
        delete children[i];
    }
}

int QuadTree::Node::getQuadrant(float x, float y) const {
    float midX = (bounds[0] + bounds[1]) / 2;
    float midY = (bounds[2] + bounds[3]) / 2;
    
    int quad = 0;
    if (x >= midX) quad |= 1;  // East
    if (y >= midY) quad |= 2;  // North
    return quad; // 0=SW, 1=SE, 2=NW, 3=NE
}

void QuadTree::Node::insert(int particleIndex, const ParticlesTable& particles) {
    if (isLeaf()) {
        particleIndices.push_back(particleIndex);
        
        // If this leaf now has more than 1 particle, subdivide
        if (particleIndices.size() > 1) {
            // Check if we can still subdivide (avoid infinite recursion for identical positions)
            float width = bounds[1] - bounds[0];
            float height = bounds[3] - bounds[2];
            const float MIN_SIZE = 1e-8f;  // Minimum cell size to prevent infinite subdivision
            
            if (width > MIN_SIZE && height > MIN_SIZE) {
                // Create children
                float midX = (bounds[0] + bounds[1]) / 2;
                float midY = (bounds[2] + bounds[3]) / 2;
                
                children[0] = new Node(bounds[0], midX, bounds[2], midY); // SW
                children[1] = new Node(midX, bounds[1], bounds[2], midY); // SE
                children[2] = new Node(bounds[0], midX, midY, bounds[3]); // NW  
                children[3] = new Node(midX, bounds[1], midY, bounds[3]); // NE
                
                // Redistribute particles to children using insert (for recursive subdivision)
                std::vector<int> tempIndices = particleIndices;
                particleIndices.clear();
                
                for (int idx : tempIndices) {
                    int quad = getQuadrant(particles.x[idx], particles.y[idx]);
                    children[quad]->insert(idx, particles);
                }
            }
            // If cell is too small, keep particles in this leaf
        }
    } else {
        // Internal node - pass to appropriate child
        int quad = getQuadrant(particles.x[particleIndex], particles.y[particleIndex]);
        children[quad]->insert(particleIndex, particles);
    }
}

void QuadTree::Node::calculateCenterOfMass(const ParticlesTable& particles) {
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

// QuadTree implementation
QuadTree::QuadTree() : root(nullptr), theta(0.5f) {}

QuadTree::QuadTree(float openingAngle) : root(nullptr), theta(openingAngle) {}

QuadTree::~QuadTree() {
    clear();
}

void QuadTree::buildFromParticles(const ParticlesTable& particles) {
    // Clear existing tree
    clear();
    
    if (particles.N == 0) return;
    
    // Find bounding box of all particles
    float xmin = particles.x[0], xmax = particles.x[0];
    float ymin = particles.y[0], ymax = particles.y[0];
    
    for (int i = 1; i < particles.N; i++) {
        if (particles.x[i] < xmin) xmin = particles.x[i];
        if (particles.x[i] > xmax) xmax = particles.x[i];
        if (particles.y[i] < ymin) ymin = particles.y[i]; 
        if (particles.y[i] > ymax) ymax = particles.y[i];
    }
    
    // Add small padding
    float padding = 0.1f;
    xmin -= padding; xmax += padding;
    ymin -= padding; ymax += padding;
    
    // Create root node
    root = new Node(xmin, xmax, ymin, ymax);
    
    // Insert all particles
    for (int i = 0; i < particles.N; i++) {
        root->insert(i, particles);
    }
    
    // Calculate centers of mass
    root->calculateCenterOfMass(particles);
}

QuadTree QuadTree::buildQuadTree(const ParticlesTable& particles) {
    QuadTree tree;
    tree.buildFromParticles(particles);
    return tree;
}

void QuadTree::clear() {
    if (root) {
        clearNode(root);
        root = nullptr;
    }
}

void QuadTree::clearNode(Node* node) {
    delete node;  // Destructor handles children recursively
}

void QuadTree::saveToFile(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open " << filename << " for writing\n";
        return;
    }
    
    file << "QUADTREE\n";
    file << "# Format: depth xmin xmax ymin ymax centerX centerY totalMass numParticles [particleIndices...]\n";
    
    if (root) {
        saveNode(file, root, 0);
    }
    
    file.close();
}

void QuadTree::saveNode(std::ofstream& file, Node* node, int depth) const {
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
            saveNode(file, node->children[i], depth + 1);
        }
    }
}