#pragma once

#include <string>
#include <iostream>
#include <cstdlib>
#include "UnitsTable.hpp"
#include "ParticlesSetup.hpp"
#include "InitialConditionSetup.hpp"

class SimulationSetup{
public:
    std::string paramspath;
    // Simulation setup
    std::string input_file;                                                 // File for reading (Update whenever extract new dumpfile)
    float tmax;                                                             // Max simulation time (IN CODE UNIT)
    float bhTreeTheta = 0.2;                                                // Angle of BH Tree (Suggestion: 0.01 < theta < 1.0)
    int dt_substepsmax = 1;                                                 // Max number of substeps per time step (Current No used)
    int num_per_dump = 10;                                                  // Dump output data per given timestep.
    int a_mode = 0;                                                         // Mode for calculate acceleration (0 => direct N-body, 1 => BHTree)
    bool print_internal = 0;                                                // Dump internel vector in `ParticlesTable` (e.g. _ax, _ay, _az...) (0 => Don't print, 1 => print)

    // CPU setup
    int OMP_NUM_THREAD = 1;                                                 // Number of OpenMP threads                                
    
    // GPU setup
    bool use_GPU = false;                                                   // Enable GPU
    int BLOCK_SIZE = 32;                                                    // GPU block size



    // Method declaration
    // Constructor
    SimulationSetup(const std::string& paramspath){
        if (std::filesystem::exists(paramspath)) {
            std::cout << "Reading parameters from " << paramspath << std::endl;
            _read_params_toml(paramspath);
        } else {
            std::cerr << "ParametersfileNotfound: The parameters file "+ paramspath +" is not in the current location. Please check again!" << std::endl;
            std::exit(1);
        }
    }
    // Destructor
    ~SimulationSetup() = default;

    /*
        static void generate_parameters_file()
    Generate a `new` parameters file for simulation

    ## Input
        - ParticlesSetup& setup: `ParticlesSetup` object for generate files
        - float dtref: The reference dt for estimating tmax
    */
   static void generate_parameters_file(const ParticlesSetup& setup, const float dtref);
    
    /*
        void make_parameters_file()
    Generate a parameters file for simulation

    */
   void make_parameters_file();

   /*
        int extract_current_index()
    Extract the current timestep from "paramspath"

    ## Output
        - int : Current index of stepping.
    */
   int extract_current_index();


protected:

    template <typename T>
    static void _write_toml_kvc(std::ofstream& fout, const std::string& key, const T& value, const std::string& comment,
                                int key_width = 14, int eq_pos = 18, int val_pos = 38, int comment_pos = 50); 

    

    /*
        void _read_params_toml(const std::string& paramspath) const;
    (Internal Method) Read a TOML parameters file for simulation parameters

    ## Input
        - std::string paramsfilepath: Name of input
    */
    void _read_params_toml(const std::string& paramsfilepath);
};