#pragma once

#include <vector>
#include <string>
#include <fstream>

// Forward declaration to avoid circular dependency
class ParticlesTable;

/*
    class OctTree
    
3D Barnes-Hut tree implementation for efficient N-body calculations.
Organizes particles in a hierarchical spatial structure.
*/
class OctTree {
public:
    /*
        struct Node
    Represents a node in the octree (either internal or leaf)
    */
    struct Node {
        float bounds[6];           // xmin, xmax, ymin, ymax, zmin, zmax
        float centerOfMass[3];     // x, y, z center of mass
        float totalMass;           // Total mass in this node
        Node* children[8];         // 8 octants (nullptr for empty)
        std::vector<int> particleIndices;  // Particle indices in this leaf
        
        Node(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax);
        ~Node();
        bool isLeaf() const { return children[0] == nullptr; }
        int getOctant(float x, float y, float z) const;
        void insert(int particleIndex, const ParticlesTable& particles);
        void calculateCenterOfMass(const ParticlesTable& particles);
    };

private:
    Node* root;
    float theta;  // Opening angle parameter

public:
    // Constructor/Destructor
    OctTree();
    explicit OctTree(float openingAngle);
    ~OctTree();

    // Core functionality
    /*
        void buildFromParticles(const ParticlesTable& particles)
    Build octree from particle positions (instance method approach)
    */
    void buildFromParticles(const ParticlesTable& particles);
    
    /*
        static OctTree buildOctTree(const ParticlesTable& particles)
    Static factory method to create and build tree in one call
    */
    static OctTree buildOctTree(const ParticlesTable& particles);

    /*
        void clear()
    Clear the tree and free all memory
    */
    void clear();

    /*
        void setTheta(float openingAngle)
    Set Barnes-Hut opening angle parameter
    */
    void setTheta(float openingAngle) { theta = openingAngle; }
    float getTheta() const { return theta; }

    /*
        void saveToFile(const std::string& filename) const
    Save tree structure to text file for visualization
    */
    void saveToFile(const std::string& filename) const;

    /*
        Node* getRoot() const
    Get root node for tree traversal
    */
    Node* getRoot() const { return root; }

private:
    void clearNode(Node* node);
    void saveNode(std::ofstream& file, Node* node, int depth) const;
};