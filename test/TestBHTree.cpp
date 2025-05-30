#include <iostream>
#include <random>
#include <cassert>
#include <cmath>
#include <filesystem>
#include "ParticlesTable.hpp"
#include "UnitsTable.hpp"
#include "QuadTree.hpp"
#include "OctTree.hpp"

int main(){
    std::cout << "\n=== Barnes-Hut Tree End-to-End Test ===\n\n";

    // 1) Create 100 particles in a small box [-0.1,0.1]^3
    std::cout << "[1] Generating 100 particles in [-0.1,0.1]^3…\n";
    UnitsTable units(1.0f, 1.0f);
    const int N = 100;
    ParticlesTable pt(units, N);
    std::mt19937 rng{42};
    std::uniform_real_distribution<float> pos(-0.1f, 0.1f);
    std::uniform_real_distribution<float> mass_dist(0.8f, 1.2f);

    float totalMass = 0;
    for(int i=0; i<N; ++i){
        pt.x[i] = pos(rng);
        pt.y[i] = pos(rng);
        pt.z[i] = pos(rng);
        pt.m[i] = mass_dist(rng);
        totalMass += pt.m[i];
        // v, h, dt can be left default or zero
    }
    pt.Mtot = totalMass;
    std::cout << "    total mass = " << totalMass << "\n";

    // 2) Write to HDF5 and read back
    const std::string h5file = "bhtree_particles.h5";
    std::cout << "[2] Writing to " << h5file << "…\n";
    pt.extract_particles_table(h5file);

    std::cout << "[3] Reading back…\n";
    ParticlesTable pt2 = ParticlesTable::read_particles_table(h5file);
    float mass2 = 0;
    for(int i=0; i<pt2.N; ++i) mass2 += pt2.m[i];
    assert(std::abs(totalMass - mass2) < 1e-6f);
    std::cout << "    verified mass = " << mass2 << "\n";

    // 3) Build QuadTree in three ways & dump
    std::cout << "[4] Building & dumping QuadTrees…\n";
    {
        // a) instance
        auto qt1 = pt2.buildQuadTree();
        qt1.saveToFile("quadtree_instance.txt");

        // b) static
        auto qt2 = QuadTree::buildQuadTree(pt2);
        qt2.saveToFile("quadtree_static.txt");

        // c) direct
        QuadTree qt3;
        qt3.buildFromParticles(pt2);
        qt3.saveToFile("quadtree_direct.txt");
    }

    // 4) Build OctTree in three ways & dump
    std::cout << "[5] Building & dumping OctTrees…\n";
    {
        // a) instance
        auto ot1 = pt2.buildOctTree();
        ot1.saveToFile("octree_instance.txt");

        // b) static
        auto ot2 = OctTree::buildOctTree(pt2);
        ot2.saveToFile("octree_static.txt");

        // c) direct
        OctTree ot3;
        ot3.buildFromParticles(pt2);
        ot3.saveToFile("octree_direct.txt");
    }

    std::cout << "\nDone! Files created:\n"
              << "  bhtree_particles.h5\n"
              << "  quadtree_{instance,static,direct}.txt\n"
              << "  octree_{instance,static,direct}.txt\n\n";
    return 0;
}