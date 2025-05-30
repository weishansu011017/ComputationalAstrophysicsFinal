#pragma once

#include "UnitsTable.hpp"
#include <vector>
#include <string>
#include <unordered_map>

// Forward declarations
class QuadTree;
class OctTree;

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
        : unittable(unit), N(N_)
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
    ~ParticlesTable() = default;

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
        void calculate_h();

    Calculate the smoothing radius of each particles
    */
    void calculate_h();

    /*
        void calculate_a_dirnbody()

    Calculate dt of each particles.
    */
    void calculate_dt();

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

    /*
        void calculate_a_BHtree_2D()
    Calculate the acceleration by Barnes–Hut tree in 2D (using QuadTree)
    */
    void calculate_a_BHtree_2D();

    /*
        void kick();
    Do a kick operation to all of the active particles
    ## Input
        - float scale: scale of dt (stepping ... + scale * dt)
    */
    void kick(float scale = 1.0);

    /*
        void drift();
    Do a drift operation to all of the active particles
    ## Input
        - float scale: scale of dt (stepping ... + scale * dt)
    */
    void drift(float scale = 1.0);

    /*
        void particles_validation();
        
    Performs sanity checks and aggressive validation on particle data.

    */        
    void particles_validation();
  
    /*
        QuadTree buildQuadTree() const
    Build a quadtree from current particle positions
    Returns a QuadTree object
    */
    QuadTree buildQuadTree() const;

    /*
        OctTree buildOctTree() const
    Build an octree from current particle positions
    Returns an OctTree object
    */
    OctTree buildOctTree() const;

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
};