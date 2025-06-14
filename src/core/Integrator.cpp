#include <functional>
#include "Integrator.hpp"
#include "ParticlesTable.hpp"
#include "SimulationSetup.hpp"

Integrator wrap_Integrator(ParticlesTable* pt, SimulationSetup* simsetup) {
    Integrator integrator;
    if (simsetup->a_mode == 0) {
        if (pt->dimension == 3) {
            #ifdef USE_CUDA
            if (simsetup->use_GPU){
                integrator.calculate_a = [pt]() { pt->calculate_a_dirnbody_gpu(); };
                integrator.kick        = [pt](float scale) { pt->kick_gpu(scale); };
                integrator.drift       = [pt](float scale) { pt->drift_gpu(scale); };
            } else {
                integrator.calculate_a = [pt]() { pt->calculate_a_dirnbody(); };
                integrator.kick        = [pt](float scale) { pt->kick(scale); };
                integrator.drift       = [pt](float scale) { pt->drift(scale); };
            }
            return integrator;
            #endif
            integrator.calculate_a = [pt]() { pt->calculate_a_dirnbody(); };
            integrator.kick        = [pt](float scale) { pt->kick(scale); };
            integrator.drift       = [pt](float scale) { pt->drift(scale); };
            
        } else if (pt->dimension == 2) {
            #ifdef USE_CUDA
            if (simsetup->use_GPU){
                integrator.calculate_a = [pt]() { pt->calculate_a_dirnbody_2D_gpu(); };
                integrator.kick        = [pt](float scale) { pt->kick_2D_gpu(scale); };
                integrator.drift       = [pt](float scale) { pt->drift_2D_gpu(scale); };
            } else {
                integrator.calculate_a = [pt]() { pt->calculate_a_dirnbody_2D(); };
                integrator.kick        = [pt](float scale) { pt->kick(scale); };
                integrator.drift       = [pt](float scale) { pt->drift(scale); };
            }
            return integrator;
            #endif
            integrator.calculate_a = [pt]() { pt->calculate_a_dirnbody_2D(); };
            integrator.kick        = [pt](float scale) { pt->kick_2D(scale); };
            integrator.drift       = [pt](float scale) { pt->drift_2D(scale); };
        } else {
            std::cerr << "Unsupported dimension (should be 2 or 3)" << std::endl;
            std::exit(1);
        }
    } else if (simsetup->a_mode == 1) {
        if (pt->dimension == 3) {
            #ifdef USE_CUDA
            if (simsetup->use_GPU){
                integrator.calculate_a = [pt]() { pt->calculate_a_BHtree_gpu(); };
                integrator.kick        = [pt](float scale) { pt->kick_gpu(scale); };
                integrator.drift       = [pt](float scale) { pt->drift_gpu(scale); };
            } else {
                integrator.calculate_a = [pt]() { pt->calculate_a_BHtree(); };
                integrator.kick        = [pt](float scale) { pt->kick(scale); };
                integrator.drift       = [pt](float scale) { pt->drift(scale); };
            }
            return integrator;
            #endif
            integrator.calculate_a = [pt]() { pt->calculate_a_BHtree(); };
            integrator.kick        = [pt](float scale) { pt->kick(scale); };
            integrator.drift       = [pt](float scale) { pt->drift(scale); };
        } else if (pt->dimension == 2) {
            #ifdef USE_CUDA
            if (simsetup->use_GPU){
                integrator.calculate_a = [pt]() { pt->calculate_a_BHtree_2D_gpu(); };
                integrator.kick        = [pt](float scale) { pt->kick_2D_gpu(scale); };
                integrator.drift       = [pt](float scale) { pt->drift_2D_gpu(scale); };
            } else {
                integrator.calculate_a = [pt]() { pt->calculate_a_BHtree_2D(); };
                integrator.kick        = [pt](float scale) { pt->kick(scale); };
                integrator.drift       = [pt](float scale) { pt->drift(scale); };
            }
            return integrator;
            #endif
            integrator.calculate_a = [pt]() { pt->calculate_a_BHtree_2D(); };
            integrator.kick        = [pt](float scale) { pt->kick_2D(scale); };
            integrator.drift       = [pt](float scale) { pt->drift_2D(scale); };
        } else {
           std::cerr << "Unsupported dimension (should be 2 or 3)" << std::endl;
            std::exit(1);
        }
    } else {
        std::cerr << "Unsupported a_mode (should be 0 or 1)" << std::endl;
        std::exit(1);
    }

    return integrator;
}