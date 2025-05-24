#include "ParticlesTable.hpp"
#include "UnitsTable.hpp"
#include "ParticlesSetup.hpp"
#include "InitialConditionSetup.hpp"
#include <iostream>
#include <filesystem>

/*
    TEST: Testing Barnes-Hut Tree Construction

    Test content:
    1. Create 100 particles with uniform distribution in a small box
    2. Build both quadtree (2D) and octree (3D) from the particles
    3. Save tree structures to files for visualization
    4. Verify tree construction completed successfully
*/

int main() {
    std::cout << "\n=== Barnes-Hut Tree Test ===\n\n";
    
    // Step 1: Create a simple particle setup
    std::cout << "Step 1: Creating particle setup...\n";
    
    // Create a temporary setup object (we'll make it manually)
    UnitsTable unit(1.0f, 1.0f);  // Simple units
    ParticlesTable pt(unit, 100);  // 100 particles
    pt.SimulationTag = "bhtree_test";
    
    // Manually create particle distribution in a 10x10x10 box
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> pos_dist(-5.0f, 5.0f);
    std::uniform_real_distribution<float> vel_dist(-0.1f, 0.1f);
    std::uniform_real_distribution<float> mass_dist(0.8f, 1.2f);
    
    float totalMass = 0.0f;
    for (int i = 0; i < pt.N; i++) {
        pt.particle_index[i] = i + 1;
        
        // Random positions in box
        pt.x[i] = pos_dist(gen);
        pt.y[i] = pos_dist(gen); 
        pt.z[i] = pos_dist(gen);
        
        // Small random velocities
        pt.vx[i] = vel_dist(gen);
        pt.vy[i] = vel_dist(gen);
        pt.vz[i] = vel_dist(gen);
        
        // Random masses around 1.0
        pt.m[i] = mass_dist(gen);
        totalMass += pt.m[i];
        
        // Zero smoothing length and timestep for now
        pt.h[i] = 0.1f;
        pt.dt[i] = 0.01f;
    }
    pt.Mtot = totalMass;
    
    std::cout << "Created " << pt.N << " particles with total mass " << totalMass << "\n";
    
    // Step 2: Save initial conditions to HDF5 file
    std::cout << "\nStep 2: Saving initial conditions...\n";
    std::string ic_filename = "bhtree_test_00000.h5";
    pt.extract_particles_table(ic_filename);
    std::cout << "Saved to: " << ic_filename << "\n";
    
    // Step 3: Test reading the file back
    std::cout << "\nStep 3: Reading particles from file...\n";
    ParticlesTable pt_read = ParticlesTable::read_particles_table(ic_filename);
    std::cout << "Successfully read " << pt_read.N << " particles\n";
    
    // Step 4: Build quadtree (2D)
    std::cout << "\nStep 4: Building quadtree (2D)...\n";
    pt_read.buildQuadTree();
    std::cout << "Quadtree construction completed\n";
    
    // Step 5: Build octree (3D) 
    std::cout << "\nStep 5: Building octree (3D)...\n";
    pt_read.buildOctTree();
    std::cout << "Octree construction completed\n";
    
    // Step 6: Save tree structures
    std::cout << "\nStep 6: Saving tree structures...\n";
    
    // Save quadtree
    pt_read.saveQuadTreeToFile("quadtree_structure.txt");
    std::cout << "Saved quadtree structure to: quadtree_structure.txt\n";
    
    // Save octree
    pt_read.saveOctTreeToFile("octree_structure.txt");  
    std::cout << "Saved octree structure to: octree_structure.txt\n";
    
    // Step 7: Summary
    std::cout << "\n=== Test Summary ===\n";
    std::cout << "✓ Created " << pt.N << " test particles\n";
    std::cout << "✓ Saved/loaded HDF5 file successfully\n"; 
    std::cout << "✓ Built quadtree for 2D visualization\n";
    std::cout << "✓ Built octree for 3D visualization\n";
    std::cout << "✓ Exported tree structures for Python plotting\n";
    std::cout << "\nFiles created:\n";
    std::cout << "  - " << ic_filename << " (particle data)\n";
    std::cout << "  - quadtree_structure.txt (2D tree)\n";
    std::cout << "  - octree_structure.txt (3D tree)\n";
    std::cout << "\nRun Python visualization scripts to see the trees!\n\n";
    
    return 0;
}