#pragma once
typedef struct CUstream_st* cudaStream_t;
#include <vector>
#include <string>
#include "QuadTree.hpp"
#include "OctTree.hpp"

class UnitsTable;

/*
    class ParticlesTable

The main container for particles data.
*/
class ParticlesTable{
public:
    // Variable declaration
    // Simulation condition
    const UnitsTable& unittable;                                      // Code unit (to cgs).
    
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
    std::vector<float> _U;                                      // Potential energy

    // Other scalar parameter(External)
    int N;
    int dimension;
    float t = 0.0;
    float Mtot = 0.0;
    float Utot = 0.0;
    float bhTreeTheta = 0.0f;
    std::string SimulationTag;

    // Method declaration
    // Constructor
    ParticlesTable(const UnitsTable& unit , int N_);
    
    // Destructor
    ~ParticlesTable() = default;

    // Other method

    /*
        void extract_particles_table(const std::string& filename) const;

    Extract all the particle into a HDF5 file

    ## Input
        - string filename: Name of output
    */
    void extract_particles_table(const std::string& filename, bool debug = false) const;

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
        void calculate_dt()

    Calculate dt of each particles.
    */
    void calculate_dt();

    /*
        void calculate_Utot()

    Calculate the total potential energy Utot and store into `Utot`
    */
    void calculate_Utot();

    /*
        void calculate_a_dirnbody()

    Calculate the acceleration by direct N-body method
    */
    void calculate_a_dirnbody();

    /*
        void calculate_a_dirnbody_2D()

    Calculate the acceleration by direct N-body method in 2D
    */
    void calculate_a_dirnbody_2D();

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
        void kick_2D();
    Do a kick operation to all of the active particles
    ## Input
        - float scale: scale of dt (stepping ... + scale * dt)
    */
    void kick_2D(float scale = 1.0);

    /*
        void drift();
    Do a drift operation to all of the active particles
    ## Input
        - float scale: scale of dt (stepping ... + scale * dt)
    */
    void drift_2D(float scale = 1.0);

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

    /*
        void set_bhTreeTheta(float theta);

    Set the internel variable `bhTreeTheta` to the given value

    ## Input
        - float theta: The new bhTreeTheta
    */
    void set_bhTreeTheta(float theta);


    //=== CUDA GPU-related===//
    // Define pointer of array
    float *d_x=nullptr, *d_y=nullptr, *d_z=nullptr;
    float *d_vx=nullptr, *d_vy=nullptr, *d_vz=nullptr;
    float *d_m=nullptr,  *d_h=nullptr, *d_dt=nullptr;
    float *d_ax=nullptr,*d_ay=nullptr,*d_az=nullptr,*d_U=nullptr;

    int block = 256;
    cudaStream_t stream = 0;         
    bool gpu_init = false;
    QuadNode* d_nodes_2D = nullptr;
    OctNode* d_nodes_3D = nullptr;
    int* d_order = nullptr;

    /*
        void device_init();

    Initialize GPU device memory and allocate necessary buffers.
    This function must be called before any GPU operations.
    */    
    void device_init();

    /*
        void device_finalize();

    Free all GPU memory previously allocated by `device_init`.
    */
    void device_finalize();

    /*
        void upload_all();

    Upload all particle data from host to device memory.
    */
    void upload_all();

    /*
        void download_state();

    Download all relevant particle state from device to host.
    */
    void download_state();


    /*
        void calculate_a_dirnbody_gpu();

    Calculate the acceleration by direct N-body method using CUDA on GPU.
    */
    void calculate_a_dirnbody_gpu();

    /*
        void calculate_a_dirnbody_2D_gpu();

    Calculate the acceleration by direct N-body method using CUDA on GPU in 2D
    */
    void calculate_a_dirnbody_2D_gpu();

    /*
        void calculate_a_BHtree_gpu();

    Calculate the acceleration by Barnes–Hut tree
    */
    void calculate_a_BHtree_gpu();

    /*
        void calculate_a_BHtree_2D_gpu();

    Calculate the acceleration by Barnes–Hut tree in 2D
    */
    void calculate_a_BHtree_2D_gpu();

    /*
        void kick_gpu(float scale);

    Perform a kick operation on GPU (update velocity using acceleration).

    ## Input
        - float scale: scale of dt (applied as v += scale * dt * a)
    */
    void kick_gpu(float scale);

    /*
        void drift_gpu(float scale);

    Perform a drift operation on GPU (update position using velocity).

    ## Input
        - float scale: scale of dt (applied as x += scale * dt * v)
    */
    void drift_gpu(float scale);

    /*
        void kick_2D_gpu(float scale);

    Perform a kick operation on GPU (update velocity using acceleration) in 2D

    ## Input
        - float scale: scale of dt (applied as v += scale * dt * a)
    */
    void kick_2D_gpu(float scale);

    /*
        void drift_2D_gpu(float scale);

    Perform a drift operation on GPU (update position using velocity) in 2D.

    ## Input
        - float scale: scale of dt (applied as x += scale * dt * v)
    */
    void drift_2D_gpu(float scale);



protected:
    // Allocating vector
    void _resize_vectors(std::size_t N);

    // === Base I/O Methods ===

    /*
        void _write_base_HDF5(hid_t file_id, bool debug) const;
    (Internal Method) Write the base properties of particles to HDF5
    ## Input
        - hid_t file_id: The HDF5 file handle created by H5Fcreate or H5Fopen
    */
    void _write_base_HDF5(long long file_id, bool debug) const;

    /*
        void _read_base_HDF5(hid_t file_id);
    (Internal Method) Read the properties of particles from HDF5
    ## Input
        - hid_t file_id: The HDF5 file handle created by H5Fcreate or H5Fopen
    */
    void _read_base_HDF5(long long file_id);
    
    /*
        void _calculate_a_OctNode(int idx, const OctTree& tree, int nidx);

    (Internel method) Calculate acc from a octree node for a sigle particles
    */
    void _calculate_a_OctNode(int idx, const OctTree& tree, int nidx);


    /*
        void _calculate_a_QuadNode(int idx, const QuadTree& tree, int nidx)

    (Internel method) Calculate acc from a quadtree node for a sigle particles in 2D
    */
    void _calculate_a_QuadNode(int idx, const QuadTree& tree, int nidx);
};