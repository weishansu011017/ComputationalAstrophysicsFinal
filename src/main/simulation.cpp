#include <string>
#include <iostream>
#include <iomanip> 
#include <toml++/toml.h>
#include <omp.h>
#include <functional>
#include <chrono>
#include "ParticlesSetup.hpp"
#include "UnitsTable.hpp"
#include "ParticlesTable.hpp"
#include "PhysicalConstants.hpp"
#include "SimulationSetup.hpp"

int main(int argc, char** argv){
    // Checking arguments
    if (argc < 2) {
        std::cerr << "Usage: ./simulation paramsfile \n";
        return 1;
    }
    std::cout << "\n";
    std::cout << " Barnes-Hut tree-based N body simulations\n";
    std::cout << "     Version 0.0.1\n";

    // Read simulation setup
    std::string paramsfile = argv[1];
    SimulationSetup simsetup = SimulationSetup(paramsfile);

    // Setup parallel stuff (OpenMP and CUDA).
    omp_set_num_threads(simsetup.OMP_NUM_THREAD);
    std::cout << "Assigning " << simsetup.OMP_NUM_THREAD << " threads to OpenMP." << std::endl;

    // Read current dump file
    int current_idt = simsetup.extract_current_index();
    ParticlesTable pt = ParticlesTable::read_particles_table(simsetup.input_file);

    // Configure Barnes-Hut if needed
    if (simsetup.a_mode == 1) {
        // Optional: allow theta to be set in config file
        // pt.setBHTreeTheta(simsetup.bh_theta);  // Add bh_theta to SimulationSetup
        pt.setBHTreeTheta(0.5f);  // Default value
        std::cout << "Using Barnes-Hut tree with theta = " << pt.bhTreeTheta << std::endl;
    }

    // Specify the calculate_a() function
    std::function<void()> ptcalculate_a;
    if (simsetup.a_mode == 0){
        ptcalculate_a = [&pt]() {
            pt.calculate_a_dirnbody();
        };
    } else if (simsetup.a_mode == 1) {
        ptcalculate_a = [&pt]() {
            pt.calculate_a_BHtree();
        };
    } else {
        throw std::runtime_error("The a_mode " + std::to_string(simsetup.a_mode) + 
                                " is not supported! Use 0 for direct, 1 for Barnes-Hut.");
    }

    // Other option
    // Calculate initial energy for conservation check
    float initial_kinetic = 0.0f;
    float initial_potential = 0.0f;
    
    for (int i = 0; i < pt.N; i++) {
        float v2 = pt.vx[i]*pt.vx[i] + pt.vy[i]*pt.vy[i] + pt.vz[i]*pt.vz[i];
        initial_kinetic += 0.5f * pt.m[i] * v2;
    }
    
    // Approximate potential energy (using forces)
    for (int i = 0; i < pt.N; i++) {
        float r = std::sqrt(pt.x[i]*pt.x[i] + pt.y[i]*pt.y[i] + pt.z[i]*pt.z[i]);
        initial_potential -= pt.m[i] * pt.Mtot / (r + pt.h[i]);
    }
    
    float initial_energy = initial_kinetic + initial_potential;
    std::cout << "Initial energy: E = " << initial_energy 
              << " (K = " << initial_kinetic << ", U = " << initial_potential << ")" << std::endl;

    // Main iteration
    int iter = 0;
    std::cout << "============================= Start simulation =============================" << std::endl;
    while (true) {
        if (pt.t >= simsetup.tmax){
            break;
        }
        // Timer start 
        auto start = std::chrono::high_resolution_clock::now();

        // KDK algorithm
        ptcalculate_a();
        pt.kick(0.5);
        pt.drift(1.0);
        ptcalculate_a();
        pt.kick(0.5);

        // Sanity check
        pt.particles_validation();

        // Timer end
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        // Other updating
        float deltat = mean(pt.dt);         // Currently, we don't have hierarchial time stepping (DO NOT USE THIS IN HIERARCHIAL TIMESTEPPING!!!!!!!!)
        // Time add & output dumpfile
        std::cout 
            << "t = " << std::scientific << std::setw(13) << std::setprecision(6) << pt.t
            << " (code unit)  =====>  "
            << std::scientific << std::setw(13) << std::setprecision(6) << pt.t + deltat
            << " (code unit), Walltime/iter = "
            << std::setw(8) << duration << " us"
            << std::endl;

        pt.t += deltat;        
        ++iter;

        // Periodic energy check (every 100 steps)
        if (iter % 100 == 0) {
            float kinetic = 0.0f;
            for (int i = 0; i < pt.N; i++) {
                float v2 = pt.vx[i]*pt.vx[i] + pt.vy[i]*pt.vy[i] + pt.vz[i]*pt.vz[i];
                kinetic += 0.5f * pt.m[i] * v2;
            }
            // Energy drift check
            float energy_error = std::abs((kinetic + initial_potential - initial_energy) / initial_energy);
            if (energy_error > 0.1) {
                std::cout << "Warning: Large energy drift detected: " << energy_error * 100 << "%" << std::endl;
            }
        }

        if (iter % simsetup.num_per_dump == 0){
            ++current_idt;
            std::string timeindex = format_index(current_idt, 5);
            std::string output = pt.SimulationTag + "_" + timeindex + ".h5";
            std::cout << "Dump: "<< output << std::endl;
            pt.extract_particles_table(output);
            simsetup.input_file = output;
            simsetup.make_parameters_file();
        }
    }
    std::cout << "============================== End simulation ==============================" << std::endl;
    std::cout << "Total iterations: " << iter << std::endl;
    return 0;
}