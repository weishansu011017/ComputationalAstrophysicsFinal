#include <string>
#include <iostream>
#include <regex>
#include "UnitsTable.hpp"
#include "SimulationSetup.hpp"
#include "InitialConditionSetup.hpp"

template <typename T>
void SimulationSetup::_write_toml_kvc(std::ofstream& fout, const std::string& key, const T& value, const std::string& comment,
                            int key_width, int eq_pos, int val_pos, int comment_pos) {
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

void SimulationSetup::generate_parameters_file(const ParticlesSetup& setup, const float dtref){
    std::string SimulationTag = setup.SimulationTag;
    std::string input_file = SimulationTag + "_00000.h5";
    std::string paramspath = SimulationTag + ".in";

    std::ofstream fout(paramspath, std::ios::out | std::ios::trunc);
    if (!fout.is_open()) {
        std::cerr << "Failed to open file for writing: " << paramspath << std::endl;
        std::exit(1);
    }
    fout << "[SimulationParameters]\n";
    _write_toml_kvc(fout, "input_file", input_file, "File for reading (Update whenever extract new dumpfile)");
    _write_toml_kvc(fout, "tmax", 4000.0 * dtref, "Max simulation time (IN CODE UNIT)");
    _write_toml_kvc(fout, "bhTreeTheta", 0.2, "Angle of BH Tree (Suggestion: 0.01 < theta < 1.0)");
    _write_toml_kvc(fout, "dt_substepsmax", 1, "Max number of substeps per time step (Current No used)");
    _write_toml_kvc(fout, "num_per_dump", 10, "Dump output data per given timestep.");
    _write_toml_kvc(fout, "a_mode", 0, "Mode for calculate acceleration (0 => direct N-body, 1 => BHTree)");
    _write_toml_kvc(fout, "print_internal", 0, "Dump internel vector in `ParticlesTable` (e.g. _ax, _ay, _az...) (0 => Don't print, 1 => print).");

    fout << "[CPUSetup]\n";
    _write_toml_kvc(fout, "OMP_NUM_THREAD", 1, "Number of OpenMP threads");

    fout << "[GPUSetup]\n";
    _write_toml_kvc(fout, "BLOCK_SIZE", 32, "GPU block size");

    fout << "\n";
    fout.close();
}

void SimulationSetup::make_parameters_file(){

    std::ofstream fout(paramspath, std::ios::out | std::ios::trunc);
    if (!fout.is_open()) {
        std::cerr << "Failed to open file for writing: " << paramspath << std::endl;
        std::exit(1);
    }
    int print_int = (print_internal) ? 1 : 0;
    fout << "[SimulationParameters]\n";
    _write_toml_kvc(fout, "input_file", input_file, "File for reading (Update whenever extract new dumpfile)");
    _write_toml_kvc(fout, "tmax", tmax, "Max simulation time (IN CODE UNIT)");
    _write_toml_kvc(fout, "bhTreeTheta", bhTreeTheta, "Angle of BH Tree (Suggestion: 0.01 < theta < 1.0)");
    _write_toml_kvc(fout, "dt_substepsmax", dt_substepsmax, " Max number of substeps per time step (Current No used)");
    _write_toml_kvc(fout, "num_per_dump", num_per_dump, "Dump output data per given timestep.");
    _write_toml_kvc(fout, "a_mode", a_mode, "Mode for calculate acceleration (0 => direct N-body, 1 => BHTree)");
    _write_toml_kvc(fout, "print_internal", print_int, "Dump internel vector in `ParticlesTable` (e.g. _ax, _ay, _az...) (0 => Don't print, 1 => print).");

    fout << "[CPUSetup]\n";
    _write_toml_kvc(fout, "OMP_NUM_THREAD", OMP_NUM_THREAD, "Number of OpenMP threads");

    fout << "[GPUSetup]\n";
    _write_toml_kvc(fout, "BLOCK_SIZE",BLOCK_SIZE, "GPU block size");

    fout << "\n";
    fout.close();
}

void SimulationSetup::_read_params_toml(const std::string& paramsfilepath){
    toml::table config;
    try {
        config = toml::parse_file(paramsfilepath);
    } catch (const toml::parse_error& err) {
        std::cerr << "TOML parsing error: " << std::string(err.description()) << std::endl;
        std::exit(1);
    }
    paramspath               = paramsfilepath;
    input_file               = config["SimulationParameters"]["input_file"].value_or(input_file);
    tmax                     = config["SimulationParameters"]["tmax"].value_or(tmax);
    bhTreeTheta              = config["SimulationParameters"]["bhTreeTheta"].value_or(bhTreeTheta);
    dt_substepsmax           = config["SimulationParameters"]["dt_substepsmax"].value_or(dt_substepsmax);
    num_per_dump             = config["SimulationParameters"]["num_per_dump"].value_or(num_per_dump);
    a_mode                   = config["SimulationParameters"]["a_mode"].value_or(a_mode);
    int print_int            = config["SimulationParameters"]["print_internal"].value_or(0);
    OMP_NUM_THREAD           = config["CPUSetup"]["OMP_NUM_THREAD"].value_or(OMP_NUM_THREAD);
    BLOCK_SIZE               = config["GPUSetup"]["BLOCK_SIZE"].value_or(BLOCK_SIZE);

    if (print_int == 0){
        print_internal = false;
    } else {
        print_internal = true;
    }

}

int SimulationSetup::extract_current_index() {
    std::regex pattern(R"(_(\d+)\.h5$)");
    std::smatch match;
    if (std::regex_search(input_file, match, pattern)) {
        return std::stoi(match[1]);
    } else {
        std::cerr << "Filename format not matched: " + input_file << std::endl;
        std::exit(1);
    }
}