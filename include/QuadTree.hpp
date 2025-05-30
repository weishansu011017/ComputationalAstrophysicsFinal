#pragma once

#include <vector>
#include <string>
#include <fstream>

// Forward declaration to avoid circular dependency
class ParticlesTable;

/*
    class QuadTree
    
2D Barnes-Hut tree implementation for efficient N-body calculations.
Organizes particles in a hierarchical spatial structure.
*/
class QuadTree {
public:
    /*
        struct Node
    Represents a node in the quadtree (either internal or leaf)
    */
    struct Node {
        float bounds[4];           // xmin, xmax, ymin, ymax
        float centerOfMass[2];     // x, y center of mass
        float totalMass;           // Total mass in this node
        Node* children[4];         // NW, NE, SW, SE (nullptr for empty)
        std::vector<int> particleIndices;  // Particle indices in this leaf
        
        Node(float xmin, float xmax, float ymin, float ymax);
        ~Node();
        bool isLeaf() const { return children[0] == nullptr; }
        int getQuadrant(float x, float y) const;
        void insert(int particleIndex, const ParticlesTable& particles);
        void calculateCenterOfMass(const ParticlesTable& particles);
    };

private:
    Node* root;
    float theta;  // Opening angle parameter

public:
    // Constructor/Destructor
    QuadTree();
    explicit QuadTree(float openingAngle);
    ~QuadTree();

    // Core functionality
    /*
        void buildFromParticles(const ParticlesTable& particles)
    Build quadtree from particle positions (static factory method approach)
    */
    void buildFromParticles(const ParticlesTable& particles);
    
    /*
        static QuadTree buildQuadTree(const ParticlesTable& particles)
    Static factory method to create and build tree in one call
    */
    static QuadTree buildQuadTree(const ParticlesTable& particles);

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