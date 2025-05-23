#pragma once

#include <string>
#include "UnitsTable.hpp"
#include "ParticlesSetup.hpp"

class SimulationSetup{
public:
    // Simulation setup
    std::string input_file;                                                 // File for reading (Update whenever extract new dumpfile)
    float tmax;                                                             // Max simulation time (IN CODE UNIT)
    bool enable_adaptive_dt = true;                                         // Enable adaptive time stepping
    int max_adaptive_dt_substeps = 3;                                       // Max number of substeps per time step

    // CPU setup
    int OMP_NUM_THREADS = 1;                                    


protected:

};