#pragma once
#include <vector>
#include <string>
#include <algorithm> 
#include <math.h>

#ifdef __CUDACC__
  #define HD __host__ __device__
#else
  #define HD
#endif

class OctNode {
public:
    // Node range
    float xmin; float xmax;
    float ymin; float ymax;
    float zmin; float zmax;
    
    // COM & Mass
    float Mtot;
    float COMx; float COMy; float COMz;

    // Other quantities
    int ParticlesLocateidx;     // Location of first particles in this node in the reordered particles list
    int pcount;                 // Particles counts (How much particles do this node contains)

    // Child index in the node list
    int children[8];             // Ordering: NE(0) -> NW(1) -> SW(2) -> SE(3); Set children[n] = -1 if no chirdren node in that Quadrant

    // Constructor
    OctNode(float _xmin, float _xmax, float _ymin, float _ymax, float _zmin, float _zmax)
        : xmin(_xmin), xmax(_xmax),
        ymin(_ymin), ymax(_ymax), zmin(_zmin), zmax(_zmax),
        Mtot(0.0f), COMx(0.0f), COMy(0.0f),
        ParticlesLocateidx(-1), pcount(0)
    {
        children[0] = children[1] = children[2] = children[3] = children[4] = children[5] = children[6] = children[7] = -1;
    }

    // Destructor
    ~OctNode() = default;

    //=== Both CPU & GPU ===//
    /*
        bool isLeaf() const

    Check whether this node is a leaf node.

    ## Output
        - `bool`
    */
   HD bool isLeaf() const {
        return (children[0] < 0) && (children[1] < 0) && (children[2] < 0) && (children[3] < 0)  && (children[4] < 0)  && (children[5] < 0)  && (children[6] < 0)  && (children[7] < 0);
    };

    /*
        int getOctrant(float x, float y, float z) const;

    Get the Octrant in this node which contains given position

    ## Input 
        - float x: x
        - float y: y
        - float z: z

    ## Output
        - `int`: The Octrant. Ordering: NEU(0) -> NWU(1) -> SWU(2) -> SEU(3) -> NED(4) -> NWD(5) -> SWD(6) -> SED(7);
    */
   HD int getOctrant(float x, float y, float z) const{
        float xmid = 0.5 * (xmin + xmax);
        float ymid = 0.5 * (ymin + ymax);
        float zmid = 0.5 * (zmin + zmax);
        if (z >= zmid){
            if (y >= ymid){
                if (x >= xmid){return 0;} 
                else {return 1;}
            } else {
                if (x >= xmid){return 3;} 
                else {return 2;}
            }
        } else {
            if (y >= ymid){
                if (x >= xmid){return 4;} 
                else {return 5;}
            } else {
                if (x >= xmid){return 8;} 
                else {return 7;}
            }
        }
    }

    /*
        bool contains(float x, float y, float z) const;

    Verify whether the given point is in this node

    ## Input 
        - float x: x
        - float y: y
        - float z: z

    ## Output
        - `bool`
    */
   HD bool contains(float x, float y, float z) const{
        return (x >= xmin && x < xmax && y >= ymin && y < ymax && z >= zmin && z < zmax);
    }

    /*
        float width() const

    Get the length of node along x-axis(xmax - xmin)
   
    ## Output
        - `float`: The width of node.
    */
   HD float width() const {
        return xmax - xmin;
    }
    /*
        float height() const

    Get the length of node along y-axis(ymax - ymin)
   
    ## Output
        - `float`: The height of node.
    */
   HD float height() const{
        return ymax - ymin;
    }

    /*
        float depth() const

    Get the length of node along z-axis(zmax - zmin)
   
    ## Output
        - `float`: The depth of node.
    */
   HD float depth() const{
        return zmax - zmin;
    }

    /*
        float cellsize() const

    Get the cellsize of node (diagonal distance)
   
    ## Output
        - `float`: cellsize of node.
    */
   HD float cellsize() const{
        float w = width();
        float h = height();
        float d = depth();
        float s = sqrtf(w * w + h * h + d * d);
        return s;
    }
   /*
        bool useMultipoleApproximation(float px, float py, float pz, float theta) const;

    Verify whether the given point (px, py, pz) is far enough to used the multopole approximation. (i.e. cellsize / distance < theta)

    ## Input
        - float px, float py, float pz: The point.

    ## Output 
        - `bool`
   */
   HD bool useMultipoleApproximation(float px, float py, float pz, float theta) const{
        float dx = COMx - px;
        float dy = COMy - py;
        float dz = COMz - pz;
        float r2 = dx*dx + dy*dy + dz*dz;
        if (r2 <= 1e-8f) return false;
        float distance = std::sqrt(r2);
        float cellSize = cellsize();
        return (cellSize / distance) < theta;
    }

    //=== CPU Only ===//
     /*  
        void subdivide(std::vector<OctNode>& nodes_list);

    Divide node into four subnodes

    ## Input
        - std::vector<OctNode>& nodes_list: The vector that contain all of the node in the tree
    */
   void subdivide(std::vector<OctNode>& nodes_list);
};

class OctTree {
public:
    // Information of tree
    std::vector<OctNode> nodes_list;        // Array of current nodes
    int root_idx;                           // Index of root node (usually 0)
    int leafNmax = 8;                       // Maximum number of particles in a sigle leaf
    float bhtheta;                          // Theta for mutipole expansion
    
    // Particles data
    std::vector<float> x; std::vector<float> y; std::vector<float> z;
    std::vector<float> m;                   // Mass of each particles
    std::vector<int> order;                 // Order of particles after constructing tree.


    // Constructor
    OctTree(float theta  = 0.7f, std::size_t est_particles = 0)
        :  bhtheta(theta)
        {
            if (est_particles) reserve_nodes(est_particles);
        }
    
    // Destructor
    ~OctTree() = default;

    /*
        void reserve_nodes(std::size_t N)

    Reserve an empty contineous memory for arrays

    ## Input
        - std::size_t N: Number of particles
    */
    void reserve_nodes(std::size_t N);

    /*
        void update_Mtot_COM();

    Calculate Mtot and COM for each node
    */
    void update_Mtot_COM();

    /*
        void build_tree(const std::vector<float>& xin, const std::vector<float>& yin, const std::vector<float>& zin, const std::vector<float>& mass)

    Build tree from given particle information

    ## Input
        - xin: x-coordinates of all particles (length N)
        - yin: y-coordinates of all particles (length N)
        - zin: z-coordinates of all particles (length N)
        - mass: mass of all particles (length N)
    */
    void build_tree(const std::vector<float>& xin, const std::vector<float>& yin, const std::vector<float>& zin, const std::vector<float>& mass);

};