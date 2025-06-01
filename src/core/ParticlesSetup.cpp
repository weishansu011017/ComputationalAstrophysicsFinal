#include <fstream>
#include <iostream>
#include <toml++/toml.h>
#include "ParticlesSetup.hpp"


// Global: // Case: ParticlesSetup

template <typename T>
void ParticlesSetup::_write_toml_kvc(std::ofstream& fout, const std::string& key, const T& value, const std::string& comment,
                            int key_width, int eq_pos, int val_pos, int comment_pos) const {
    std::ostringstream val_stream;

    if constexpr (std::is_same_v<T, std::string>)
        val_stream << "\"" << value << "\"";
    else if constexpr (std::is_same_v<T, bool>)
        val_stream << (value ? "true" : "false");
    else if constexpr (std::is_floating_point_v<T>)
        val_stream << std::scientific << std::setprecision(4) << std::showpoint << value;
    else
        val_stream << value;

    std::string val_str = val_stream.str();

    // indent
    fout << std::string(4, ' ');  // 4 space indent

    // key + padding to equal sign
    fout << std::left << std::setw(key_width) << key;

    // equal sign
    fout << " = ";

    // value (right-aligned)
    int val_field_width = val_pos - eq_pos;
    fout << std::right << std::setw(val_field_width) << val_str;

    // comment alignment
    if (!comment.empty()) {
        int curr_pos = eq_pos + val_field_width + 4; // estimated current column
        if (curr_pos < comment_pos)
            fout << std::string(comment_pos - curr_pos, ' ');
        fout << "# " << comment;
    }

    fout << '\n';
}


// Case: ParticlesSetupUniform

void ParticlesSetupUniform::_make_setupin_toml(const std::string& simulation_tag) const{
    std::string filename = simulation_tag + ".setup";
    std::ofstream fout(filename, std::ios::out | std::ios::trunc);
    if (!fout.is_open()) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        std::exit(1);
    }
    _write_toml_kvc(fout, "ICSetup", std::string("uniform"), "Initial Considion setup module: `Uniform` (DO NOT CHANGE THIS!!)");
    fout << "[SimulationParameters]\n";
    _write_toml_kvc(fout, "SimulationTag", simulation_tag, "Tag of simulation (Formatting the output into `Tag_00XXX.h5`)");
    _write_toml_kvc(fout, "N", N, "Number of particles");
    _write_toml_kvc(fout, "udist", udist, "Code unit of length in cgs");
    _write_toml_kvc(fout, "umass", umass, "Code unit of mass in cgs");
    _write_toml_kvc(fout, "softfactorx", softfactorx, "Softening factor for gravity (Suggestion: 0.01 < etax < 0.1)");
    _write_toml_kvc(fout, "tsfactor", tsfactor, "Factor for estimating timestep (Suggestion: 0.1 < etat < 1.0)");
    _write_toml_kvc(fout, "bhTreeTheta", bhTreeTheta, "Angle of BH Tree (Suggestion: 0.3 < theta < 1.0)");
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
        std::cerr << "TOML parsing error: " << std::string(err.description()) << std::endl;
        std::exit(1);
    }
    std::string icsetup = config["ICSetup"].as_string()->get();
    if (icsetup != "uniform") {
        std::cerr << "Wrong setup error: The setup file should be consistent with \"uniform\", but got \"" << icsetup << "\"" << std::endl;
        std::exit(1);
    }

    // ========= Simulation Parameters =========
    N      = config["SimulationParameters"]["N"].value_or(N);
    udist  = config["SimulationParameters"]["udist"].value_or(udist);
    umass  = config["SimulationParameters"]["umass"].value_or(umass);
    SimulationTag = config["SimulationParameters"]["SimulationTag"].value_or(SimulationTag);
    softfactorx  = config["SimulationParameters"]["softfactorx"].value_or(softfactorx);
    tsfactor  = config["SimulationParameters"]["tsfactor"].value_or(tsfactor);
    bhTreeTheta = config["SimulationParameters"]["bhTreeTheta"].value_or(bhTreeTheta);

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

    // ========= Other quantities =========
    // Simulation scale
    float dx = xmax - xmin;
    float dy = ymax - ymin;
    float dz = zmax - zmin;
    simulation_scale = sqrtf(dx * dx + dy * dy + dz * dz);
    dimension = 3;
};

// Case: ParticlesSetupIsotropic

void ParticlesSetupIsotropic::_make_setupin_toml(const std::string& simulation_tag) const{
    std::string filename = simulation_tag + ".setup";
    std::ofstream fout(filename, std::ios::out | std::ios::trunc);
    if (!fout.is_open()) {
        std::cerr <<"Failed to open file for writing: " << filename << std::endl;
        std::exit(1);
    }
    _write_toml_kvc(fout, "ICSetup", std::string("isotropic"), "Initial Considion setup module: `Isotropic` (DO NOT CHANGE THIS!!)");
    fout << "[SimulationParameters]\n";
    _write_toml_kvc(fout, "SimulationTag", simulation_tag, "Tag of simulation (Formatting the output into `Tag_00XXX.h5`)");
    _write_toml_kvc(fout, "N", N, "Number of particles");
    _write_toml_kvc(fout, "udist", udist, "Code unit of length in cgs");
    _write_toml_kvc(fout, "umass", umass, "Code unit of mass in cgs");
    _write_toml_kvc(fout, "softfactorx", softfactorx, "Softening factor for gravity (Suggestion: 0.01 < etax < 0.1)");
    _write_toml_kvc(fout, "tsfactor", tsfactor, "Factor for estimating timestep (Suggestion: 0.1 < etat < 1.0)");
    _write_toml_kvc(fout, "bhTreeTheta", bhTreeTheta, "Angle of BH Tree (Suggestion: 0.3 < theta < 1.0)");
    fout << "\n[InitialConditions]\n";
    fout << "# Initial Considion setup module: `Isotropic` (Isotropic sphere with power law distribution along spacial direction.)\n";
    _write_toml_kvc(fout, "rmax", rmax, "Maximum radius for particles IN CODE UNITS.");
    _write_toml_kvc(fout, "alpha", alpha, "Power law index r^alpha");

    _write_toml_kvc(fout, "rcx", rcx, "x-position of center of particles distribution IN CODE UNITS.");
    _write_toml_kvc(fout, "rcy", rcy, "y-position of center of particles distribution IN CODE UNITS.");
    _write_toml_kvc(fout, "rcz", rcz, "z-position of center of particles distribution IN CODE UNITS.");

    _write_toml_kvc(fout, "mmin", mmin, "Minimum mass sampling for particles IN CODE UNITS.");
    _write_toml_kvc(fout, "mmax", mmax, "Maximum mass sampling for particles IN CODE UNITS.");

    fout << "\n";
    fout.close();
};

void ParticlesSetupIsotropic::_read_setupin_toml(const std::string& filename) {
    toml::table config;
    try {
        config = toml::parse_file(filename);
    } catch (const toml::parse_error& err) {
        std::cerr << "TOML parsing error: " << std::string(err.description()) << std::endl;
        std::exit(1);
    }
    std::string icsetup = config["ICSetup"].as_string()->get();
    if (icsetup != "isotropic") {
        std::cerr << "Wrong setup error: The setup file should be consistent with \"isotropic\", but got \"" << icsetup << "\"" << std::endl;
        std::exit(1);
    }

    // ========= Simulation Parameters =========
    N      = config["SimulationParameters"]["N"].value_or(N);
    udist  = config["SimulationParameters"]["udist"].value_or(udist);
    umass  = config["SimulationParameters"]["umass"].value_or(umass);
    SimulationTag = config["SimulationParameters"]["SimulationTag"].value_or(SimulationTag);
    softfactorx  = config["SimulationParameters"]["softfactorx"].value_or(softfactorx);
    tsfactor  = config["SimulationParameters"]["tsfactor"].value_or(tsfactor);
    bhTreeTheta = config["SimulationParameters"]["bhTreeTheta"].value_or(bhTreeTheta);

    // ========= Initial Conditions =========
    rmax = config["InitialConditions"]["rmax"].value_or(rmax);
    alpha = config["InitialConditions"]["alpha"].value_or(alpha);
    rcx = config["InitialConditions"]["rcx"].value_or(rcx);
    rcy = config["InitialConditions"]["rcy"].value_or(rcy);
    rcz = config["InitialConditions"]["rcz"].value_or(rcz);

    mmin = config["InitialConditions"]["mmin"].value_or(mmin);
    mmax = config["InitialConditions"]["mmax"].value_or(mmax);

    // ========= Other quantities =========
    // Simulation scale
    simulation_scale = 2 * rmax;
    dimension = 3;
};

// Case: ParticlesSetupPlummer

void ParticlesSetupPlummer::_make_setupin_toml(const std::string& simulation_tag) const{
    std::string filename = simulation_tag + ".setup";
    std::ofstream fout(filename, std::ios::out | std::ios::trunc);
    if (!fout.is_open()) {
        std::cerr <<"Failed to open file for writing: " << filename << std::endl;
        std::exit(1);
    }
    _write_toml_kvc(fout, "ICSetup", std::string("plummer"), "Initial Considion setup module: `Plummer` (DO NOT CHANGE THIS!!)");
    fout << "[SimulationParameters]\n";
    _write_toml_kvc(fout, "SimulationTag", simulation_tag, "Tag of simulation (Formatting the output into `Tag_00XXX.h5`)");
    _write_toml_kvc(fout, "dimension", dimension, "Dimension of data.");
    _write_toml_kvc(fout, "N", N, "Number of particles");
    _write_toml_kvc(fout, "udist", udist, "Code unit of length in cgs");
    _write_toml_kvc(fout, "umass", umass, "Code unit of mass in cgs");
    _write_toml_kvc(fout, "softfactorx", softfactorx, "Softening factor for gravity (Suggestion: 0.01 < etax < 0.1)");
    _write_toml_kvc(fout, "tsfactor", tsfactor, "Factor for estimating timestep (Suggestion: 0.1 < etat < 1.0)");
    _write_toml_kvc(fout, "bhTreeTheta", bhTreeTheta, "Angle of BH Tree (Suggestion: 0.3 < theta < 1.0)");
    fout << "\n[InitialConditions]\n";
    fout << "# Initial Condition setup module: `Plummer` (Isotropic sphere with with Plummer density profile along spacial direction.)\n";
    _write_toml_kvc(fout, "Mtot", Mtot, "Total mass IN CODE UNITS.");
    _write_toml_kvc(fout, "rmax", rmax, "Maximum radius for particles IN CODE UNITS.");
    _write_toml_kvc(fout, "scalefactor", scalefactor ,"Factor of scale length a = scale * rmax");

    _write_toml_kvc(fout, "rcx", rcx, "x-position of center of particles distribution IN CODE UNITS.");
    _write_toml_kvc(fout, "rcy", rcy, "y-position of center of particles distribution IN CODE UNITS.");
    _write_toml_kvc(fout, "rcz", rcz, "z-position of center of particles distribution IN CODE UNITS. (If dimension == 2 then this argument is useless)");


    fout << "\n";
    fout.close();
};

void ParticlesSetupPlummer::_read_setupin_toml(const std::string& filename) {
    toml::table config;
    try {
        config = toml::parse_file(filename);
    } catch (const toml::parse_error& err) {
        std::cerr << "TOML parsing error: " << std::string(err.description()) << std::endl;
        std::exit(1);
    }
    std::string icsetup = config["ICSetup"].as_string()->get();
    if (icsetup != "plummer") {
        std::cerr << "Wrong setup error: The setup file should be consistent with \"plummer\", but got \"" << icsetup << "\"" << std::endl;
        std::exit(1);
    }

    // ========= Simulation Parameters =========
    N      = config["SimulationParameters"]["N"].value_or(N);
    udist  = config["SimulationParameters"]["udist"].value_or(udist);
    umass  = config["SimulationParameters"]["umass"].value_or(umass);
    SimulationTag = config["SimulationParameters"]["SimulationTag"].value_or(SimulationTag);
    softfactorx  = config["SimulationParameters"]["softfactorx"].value_or(softfactorx);
    tsfactor  = config["SimulationParameters"]["tsfactor"].value_or(tsfactor);
    bhTreeTheta = config["SimulationParameters"]["bhTreeTheta"].value_or(bhTreeTheta);
    dimension = config["SimulationParameters"]["dimension"].value_or(dimension);

    // ========= Initial Conditions =========
    rmax = config["InitialConditions"]["rmax"].value_or(rmax);
    Mtot = config["InitialConditions"]["Mtot"].value_or(Mtot);
    scalefactor = config["InitialConditions"][scalefactor].value_or(scalefactor);
    rcx = config["InitialConditions"]["rcx"].value_or(rcx);
    rcy = config["InitialConditions"]["rcy"].value_or(rcy);
    rcz = config["InitialConditions"]["rcz"].value_or(rcz);

    // ========= Other quantities =========
    // Simulation scale
    simulation_scale = 2 * rmax;
};