#pragma once

#include "UnitsTable.hpp"
#include <vector>
#include <string>
#include <fstream>

/*
    class ParticlesTable

The main container for particles data.
*/
class ParticlesTable{
public:
    // Variable declaration
    // Simulation condition
    UnitsTable unittable;                                       // Code unit (to cgs).
    
    // Table of particles(External)
    std::vector<uint32_t> particle_index;                       // Particles index
    std::vector<float> x;                                       // x-position
    std::vector<float> y;                                       // y-position
    std::vector<float> z;                                       // z-position
    std::vector<float> vx;                                      // x-velocity
    std::vector<float> vy;                                      // y-velocity
    std::vector<float> vz;                                      // z-velocity
    std::vector<float> m;                                       // Mass of particle
    std::vector<float> h;                                       // Smoothing length
    std::vector<float> dt;                                      // time step of particles

    // Table of particles(Internal)
    std::vector<float> _ax;                                     // x-Particle acceleration
    std::vector<float> _ay;                                     // y-Particle acceleration 
    std::vector<float> _az;                                     // z-Particle acceleration 

    // Other scalar parameter(External)
    int N;
    float t = 0.0;
    float Mtot;
    std::string SimulationTag;

    // Table of particles(Internal)
    
    // Other scalar parameter(Internal)

    // Method declaration
    // Constructor
    ParticlesTable(const UnitsTable& unit , int N_)
        : unittable(unit), N(N_),
          quadTreeRoot(nullptr), octTreeRoot(nullptr),
          bhTreeTheta(0.5f), useQuadTree(true)
    {
        _resize_vectors(N);
        std::fill(x.begin(), x.end(), 0.0f);
        std::fill(y.begin(), y.end(), 0.0f);
        std::fill(z.begin(), z.end(), 0.0f);
        std::fill(vx.begin(), vx.end(), 0.0f);
        std::fill(vy.begin(), vy.end(), 0.0f);
        std::fill(vz.begin(), vz.end(), 0.0f);
        std::fill(m.begin(),  m.end(),  0.0f);
        std::fill(h.begin(),  h.end(),  0.0f);
        std::fill(dt.begin(), dt.end(), 0.0f);
        std::fill(_ax.begin(), _ax.end(), 0.0f);
        std::fill(_ay.begin(), _ay.end(), 0.0f);
        std::fill(_az.begin(), _az.end(), 0.0f);
    }
    // Destructor
    ~ParticlesTable() {
        if (quadTreeRoot) clearQuadTree(quadTreeRoot);
        if (octTreeRoot) clearOctTree(octTreeRoot);
    }

    // Other method

    /*
        void extract_particles_table(const std::string& filename) const;

    Extract all the particle into a HDF5 file

    ## Input
        - string filename: Name of output
    */
    void extract_particles_table(const std::string& filename) const;

    /*
        ParticlesTable read_particles_table(const std::string& filename);

    Read the particle from a HDF5 file

    ## Input
        - string filename: Name of input
    
    ## Output
        - ParticlesTable pt: The `ParticlesTable` object.
    */
    static ParticlesTable read_particles_table(const std::string& filename);
    

    /*
        void calculate_a_dirnbody()

    Calculate the acceleration by direct N-body method
    */
    void calculate_a_dirnbody();

    /*
        void calculate_a_BHtree()
    Calculate the acceleration by Barnes–Hut tree
    */
    void calculate_a_BHtree();

    // === Barnes-Hut Tree Methods ===
    
    /*
        void buildQuadTree()
    Build 2D quadtree from current particle positions
    Use this for 2D simulations or 2D projections
    */
    void buildQuadTree();
    
    /*
        void buildOctTree()
    Build 3D octree from current particle positions  
    Use this for full 3D simulations
    */
    void buildOctTree();
    
    /*
        void setBHTreeTheta(float theta)
    Set opening angle parameter for Barnes-Hut approximation
    ## Input
        - float theta: Opening angle (typical values: 0.5-1.0)
                      Smaller = more accurate but slower
                      Larger = faster but less accurate
    */
    void setBHTreeTheta(float theta);
    
    /*
        void saveQuadTreeToFile(const std::string& filename) const
    Save quadtree structure to text file for visualization
    ## Input
        - string filename: Output filename
    */
    void saveQuadTreeToFile(const std::string& filename) const;
    
    /*
        void saveOctTreeToFile(const std::string& filename) const
    Save octree structure to text file for visualization
    ## Input
        - string filename: Output filename
    */
    void saveOctTreeToFile(const std::string& filename) const;

    // === Integration Methods ===
    
    /*
        void kick(float scale = 1.0)
    Do a kick operation to all of the active particles
    ## Input
        - float scale: scale of dt (stepping ... + scale * dt)
    */
    void kick(float scale = 1.0);

    /*
        void drift(float scale = 1.0)
    Do a drift operation to all of the active particles
    ## Input
        - float scale: scale of dt (stepping ... + scale * dt)
    */
    void drift(float scale = 1.0);

protected:

    // Allocating vector
    virtual void _resize_vectors(std::size_t N) {
        particle_index.resize(N);
        x.resize(N); y.resize(N); z.resize(N);
        vx.resize(N); vy.resize(N); vz.resize(N);
        m.resize(N); h.resize(N); dt.resize(N);
        _ax.resize(N),_ay.resize(N),_az.resize(N);
    }

    // === Base I/O Methods ===

    /*
        void _write_base_HDF5(hid_t file_id) const;
    (Internal Method) Write the base properties of particles to HDF5
    ## Input
        - hid_t file_id: The HDF5 file handle created by H5Fcreate or H5Fopen
    */
    void _write_base_HDF5(hid_t file_id) const;

    /*
        void _read_base_HDF5(hid_t file_id);
    (Internal Method) Read the properties of particles from HDF5
    ## Input
        - hid_t file_id: The HDF5 file handle created by H5Fcreate or H5Fopen
    */
    void _read_base_HDF5(hid_t file_id);

    // === Tree Node Structures ===
    
    /*
        struct QuadTreeNode
    2D Quadtree node for Barnes-Hut algorithm
    Each node represents a rectangular region and can have up to 4 children
    */
    struct QuadTreeNode {
        float bounds[4];           // xmin, xmax, ymin, ymax
        float centerOfMass[2];     // x, y center of mass
        float totalMass;           // Total mass in this node
        QuadTreeNode* children[4]; // NW, NE, SW, SE (nullptr for empty)
        std::vector<int> particleIndices;  // Particle indices in this leaf
        
        QuadTreeNode(float xmin, float xmax, float ymin, float ymax);
        ~QuadTreeNode();
        bool isLeaf() const { return children[0] == nullptr; }
        int getQuadrant(float x, float y) const;
        void insert(int particleIndex, ParticlesTable& particles);
        void calculateCenterOfMass(ParticlesTable& particles);
    };
    
    /*
        struct OctTreeNode
    3D Octree node for Barnes-Hut algorithm  
    Each node represents a cubic region and can have up to 8 children
    */
    struct OctTreeNode {
        float bounds[6];           // xmin, xmax, ymin, ymax, zmin, zmax
        float centerOfMass[3];     // x, y, z center of mass
        float totalMass;           // Total mass in this node
        OctTreeNode* children[8];  // 8 octants (nullptr for empty)
        std::vector<int> particleIndices;  // Particle indices in this leaf
        
        OctTreeNode(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax);
        ~OctTreeNode();
        bool isLeaf() const { return children[0] == nullptr; }
        int getOctant(float x, float y, float z) const;
        void insert(int particleIndex, ParticlesTable& particles);
        void calculateCenterOfMass(ParticlesTable& particles);
    };
    
    // === Tree Data Members ===
    QuadTreeNode* quadTreeRoot;    // Root of 2D quadtree  
    OctTreeNode* octTreeRoot;      // Root of 3D octree
    float bhTreeTheta;             // Opening angle parameter (default 0.5)
    bool useQuadTree;              // true for 2D, false for 3D
    
    // === Tree Helper Methods ===
    void clearQuadTree(QuadTreeNode* node);
    void clearOctTree(OctTreeNode* node); 
    void saveQuadTreeNode(std::ofstream& file, QuadTreeNode* node, int depth) const;
    void saveOctTreeNode(std::ofstream& file, OctTreeNode* node, int depth) const;
};