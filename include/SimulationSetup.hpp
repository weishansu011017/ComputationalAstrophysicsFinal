#pragma once

#include <string>
#include <iostream>
#include "UnitsTable.hpp"
#include "ParticlesSetup.hpp"
#include "InitialConditionSetup.hpp"

class SimulationSetup{
public:
    std::string paramspath;
    // Simulation setup
    std::string input_file;                                                 // File for reading (Update whenever extract new dumpfile)
    float tmax;                                                             // Max simulation time (IN CODE UNIT)
    int dt_substepsmax = 1;                                                 // Max number of substeps per time step (Current No used)

    // CPU setup
    int OMP_NUM_THREAD = 1;                                                 // Number of OpenMP threads                                
    
    // GPU setup
    int BLOCK_SIZE = 32;                                                    // GPU block size



    // Method declaration
    // Constructor
    SimulationSetup(const std::string& paramspath){
        if (std::filesystem::exists(paramspath)) {
            std::cout << "Reading parameters from " << paramspath << std::endl;
            _read_params_toml(paramspath);
        } else {
            throw std::runtime_error("ParametersfileNotfound: The parameters file "+ paramspath +" is not in the current location. Please check again!");
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

protected:

    template <typename T>
    static void _write_toml_kvc(std::ofstream& fout, const std::string& key, const T& value, const std::string& comment,
                                int key_width = 14, int eq_pos = 18, int val_pos = 38, int comment_pos = 50); 

    

    /*
        void _read_params_toml(const std::string& paramspath) const;
    (Internal Method) Read a TOML parameters file for simulation parameters

    ## Input
        - string paramsfilepath: Name of input
    */
    void _read_params_toml(const std::string& paramsfilepath);
};