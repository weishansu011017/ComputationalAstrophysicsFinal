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
#include "Integrator.hpp"

int main(int argc, char** argv){
    // Checking arguments
    if (argc < 2) {
        std::cerr << "Usage: ./simulation paramsfile \n";
        return 1;
    }
    std::cout << "\n";
    std::cout << " Barnes-Hut tree-based N body simulations\n";
    std::cout << "     Version 0.1.0\n";
    std::cout << "      Main Developer: Wei-Shan Su,\n";
    std::cout << "      Contributors:  Yu-Hsuan Tu, Yu-Jen Lin, (Conceptual discussion, physical modeling suggestions)";

    // Read simulation setup
    std::string paramsfile = argv[1];
    SimulationSetup simsetup = SimulationSetup(paramsfile);

    // Setup parallel stuff (OpenMP and CUDA).
    omp_set_num_threads(simsetup.OMP_NUM_THREAD);
    std::cout << "Assigning " << simsetup.OMP_NUM_THREAD << " threads to OpenMP." << std::endl;

    // Read current dump file
    int current_idt = simsetup.extract_current_index();
    ParticlesTable pt = ParticlesTable::read_particles_table(simsetup.input_file);

    // Setup theta
    pt.set_bhTreeTheta(simsetup.bhTreeTheta);

    // Specify the Integrator
    Integrator integrator = wrap_Integrator(&pt, &simsetup);

    // Other option
    if (simsetup.use_GPU){
        std::cout << "[INFO] GPU mode enabled. Initializing CUDA..." << std::endl;
        pt.device_init();
        pt.upload_all();
    }

    // Main iteration
    static bool has_dumped_initial = false;
    int iter = 0;
    std::cout << "============================= Start simulation =============================" << std::endl;
    while (true) {
        if (pt.t >= simsetup.tmax){
            std::cout << "Reached maximum simulation time (t = " << pt.t 
            << " >= tmax = " << simsetup.tmax << "). Terminating." << std::endl;
            break;
        }
        if (current_idt == 0 && !has_dumped_initial){
            integrator.calculate_a();
            if (simsetup.use_GPU){
                pt.download_state();
            }
            pt.calculate_Utot();
            std::string timeindex = format_index(current_idt, 5);
            std::string output = pt.SimulationTag + "_" + timeindex + ".h5";
            std::cout << "Rewrite and Dump: "<< output << std::endl;
            pt.extract_particles_table(output, simsetup.print_internal);
            simsetup.input_file = output;
            simsetup.make_parameters_file();
            has_dumped_initial = true;
        }
        // Timer start 
        auto start = std::chrono::high_resolution_clock::now();

        // KDK algorithm
        integrator.calculate_a();
        integrator.kick(0.5);
        integrator.drift(1.0);
        integrator.calculate_a();
        integrator.kick(0.5);

        // Sanity check
        pt.particles_validation();

        // Timer end
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

        // Other updating
        float deltat = mean(pt.dt);         // Currently, we don't have hierarchial time stepping (DO NOT USE THIS IN HIERARCHIAL TIMESTEPPING!!!!!!!!)
        // Time add & output dumpfile
        std::cout 
            << "t = " << std::scientific << std::setw(13) << std::setprecision(6) << pt.t
            << " (code unit)  =====>  "
            << std::scientific << std::setw(13) << std::setprecision(6) << pt.t + deltat
            << " (code unit), Walltime/iter = "
            << std::setw(11) << duration << " ns"
            << std::endl;

        pt.t += deltat;        
        ++iter;
        
        if (iter % simsetup.num_per_dump == 0){
            ++current_idt;
            if (simsetup.use_GPU){
                pt.download_state();
            }
            pt.calculate_Utot();
            std::string timeindex = format_index(current_idt, 5);
            std::string output = pt.SimulationTag + "_" + timeindex + ".h5";
            std::cout << "Dump: "<< output << std::endl;
            pt.extract_particles_table(output, simsetup.print_internal);
            simsetup.input_file = output;
            simsetup.make_parameters_file();
        }
    }
    std::cout << "============================== End simulation ==============================" << std::endl;
    if (simsetup.use_GPU){
        pt.device_finalize();
    }
    return 0;
}