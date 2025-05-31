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
        throw std::runtime_error("The a_mode is not support to this program!");
    }

    // Other option

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
    return 0;
}