#include <string>
#include <iostream>
#include "ParticlesSetup.hpp"
#include "InitialConditionSetup.hpp"
#include "UnitsTable.hpp"
#include "ParticlesTable.hpp"
#include "PhysicalConstants.hpp"

int main(int argc, char** argv){
    // Checking arguments
    if (argc < 3) {
        std::cerr << "Usage: ./setup setup_mode SimulationTag \n";
        return 1;
    }
    std::cout << "\n";
    std::cout << " Initial Condition Sampler\n";
    std::cout << "     Version 0.0.2\n";

    // Preparing input output argument
    std::string ICsetup = argv[1];
    std::string simulation_tag = argv[2];
    std::string outputname = simulation_tag + "_00000.h5";

    // Assign empty place for `ParticlesSetup`
    std::unique_ptr<ParticlesSetup> setupptr;

    if (ICsetup == "uniform") {
        std::cout << " Select Initial Condition mode: `uniform`\n";
        std::cout << "     (Uniform box inside a given cube with given mass sampling range)\n\n";
        setupptr = std::make_unique<ParticlesSetupUniform>(simulation_tag);
    } else if (ICsetup == "isotropic") {
        std::cout << " Select Initial Condition mode: `isotropic`\n";
        std::cout << "     (Isotropic sphere with power law distribution along spacial direction.)\n\n";
        setupptr = std::make_unique<ParticlesSetupIsotropic>(simulation_tag);
    } else {
        throw std::runtime_error("Unsupported setup mode: " + ICsetup);
    }
    // Constructing `UnitsTable` from setup
    UnitsTable units = UnitsTable::setup_units(*setupptr);

    // Sampling particles
    ParticlesTable pt = setup_initial_condition(*setupptr, units);

    // Print output log
    std::cout << "\n ===============================================================\n\n";

    std::cout << " End sampling:\n";
    std::cout << " Number of particles: " << pt.N << "\n\n";
    std::cout << "  Code units in cgs:\n";
    std::cout << std::scientific << std::setprecision(5)
              << " udist = " << pt.unittable.udist << " cm, "
              << "umass = " << pt.unittable.umass << " g, "
              << "utime = " << pt.unittable.utime << " s\n";
    float vcode = pt.unittable.udist / pt.unittable.utime;
    float ccode = CGSConstants::c / vcode;
    std::cout << "\n  with G = 1.0 (code units) and c = "
    << ccode << " (code units)\n";
    std::cout << "\n ===============================================================\n";

    // Extract into h5
    pt.extract_particles_table(outputname);

    return 0;
}