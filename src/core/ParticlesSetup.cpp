#include <fstream>
#include <iostream>
#include <toml++/toml.h>
#include "ParticlesSetup.hpp"


void ParticlesSetupUniform::_make_setupin_toml(const std::string& simulation_tag) const{
    std::string filename = simulation_tag + ".setup";
    std::ofstream fout(filename, std::ios::out | std::ios::trunc);
    if (!fout.is_open()) {
        throw std::runtime_error("Failed to open file for writing: " + filename);
    }
    _write_toml_kvc(fout, "ICSetup", std::string("uniform"), "Initial Considion setup module: `Uniform` (DO NOT CHANGE THIS!!)");
    fout << "[SimulationParameters]\n";
    _write_toml_kvc(fout, "SimulationTag", simulation_tag, "Tag of simulation (Formatting the output into `Tag_00XXX.h5`)");
    _write_toml_kvc(fout, "N", N, "Number of particles");
    _write_toml_kvc(fout, "udist", udist, "Code unit of length in cgs");
    _write_toml_kvc(fout, "umass", umass, "Code unit of mass in cgs");
    fout << "\n[InitialConditions]\n";
    fout << "# Initial Considion setup module: `Uniform` (Uniform box inside a given cube with given mass sampling range)\n";
    _write_toml_kvc(fout, "xmin", xmin, "Minimum x-position sampling for particles IN CODE UNITS.");
    _write_toml_kvc(fout, "xmax", xmax, "Maximum x-position sampling for particles IN CODE UNITS.");
    _write_toml_kvc(fout, "ymin", ymin, "Minimum y-position sampling for particles IN CODE UNITS.");
    _write_toml_kvc(fout, "ymax", ymax, "Maximum y-position sampling for particles IN CODE UNITS.");
    _write_toml_kvc(fout, "zmin", zmin, "Minimum z-position sampling for particles IN CODE UNITS.");
    _write_toml_kvc(fout, "zmax", zmax, "Maximum z-position sampling for particles IN CODE UNITS.");

    _write_toml_kvc(fout, "vxmin", vxmin, "Minimum vx-position sampling for particles IN CODE UNITS.");
    _write_toml_kvc(fout, "vxmax", vxmax, "Maximum vx-position sampling for particles IN CODE UNITS.");
    _write_toml_kvc(fout, "vymin", vymin, "Minimum vy-position sampling for particles IN CODE UNITS.");
    _write_toml_kvc(fout, "vymax", vymax, "Maximum vy-position sampling for particles IN CODE UNITS.");
    _write_toml_kvc(fout, "vzmin", vzmin, "Minimum vz-position sampling for particles IN CODE UNITS.");
    _write_toml_kvc(fout, "vzmax", vzmax, "Maximum vz-position sampling for particles IN CODE UNITS.");

    _write_toml_kvc(fout, "mmin", mmin, "Minimum mass sampling for particles IN CODE UNITS.");
    _write_toml_kvc(fout, "mmax", mmax, "Maximum mass sampling for particles IN CODE UNITS.");

    fout << "\n";
    fout.close();
};

void ParticlesSetupUniform::_read_setupin_toml(const std::string& filename) {
    toml::table config;
    try {
        config = toml::parse_file(filename);
    } catch (const toml::parse_error& err) {
        throw std::runtime_error("TOML parsing error: " + std::string(err.description()));
    }

    // ========= Simulation Parameters =========
    N      = config["SimulationParameters"]["N"].value_or(N);
    udist  = config["SimulationParameters"]["udist"].value_or(udist);
    umass  = config["SimulationParameters"]["umass"].value_or(umass);
    SimulationTag = config["SimulationParameters"]["SimulationTag"].value_or(SimulationTag);

    // ========= Initial Conditions =========
    xmin = config["InitialConditions"]["xmin"].value_or(xmin);
    xmax = config["InitialConditions"]["xmax"].value_or(xmax);
    ymin = config["InitialConditions"]["ymin"].value_or(ymin);
    ymax = config["InitialConditions"]["ymax"].value_or(ymax);
    zmin = config["InitialConditions"]["zmin"].value_or(zmin);
    zmax = config["InitialConditions"]["zmax"].value_or(zmax);

    vxmin = config["InitialConditions"]["vxmin"].value_or(vxmin);
    vxmax = config["InitialConditions"]["vxmax"].value_or(vxmax);
    vymin = config["InitialConditions"]["vymin"].value_or(vymin);
    vymax = config["InitialConditions"]["vymax"].value_or(vymax);
    vzmin = config["InitialConditions"]["vzmin"].value_or(vzmin);
    vzmax = config["InitialConditions"]["vzmax"].value_or(vzmax);

    mmin = config["InitialConditions"]["mmin"].value_or(mmin);
    mmax = config["InitialConditions"]["mmax"].value_or(mmax);
};