#pragma once
#include <string.h>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <random>
#include <vector>
#include <math.h>
#include <toml++/toml.h>
#include "SamplingFunctionsSet.hpp"
#include "MathematicsTools.hpp"

/*
Particles setup properties
*/
class ParticlesSetup{
public:
    // Number of particles
    int N = 10000;

    // Code units
    float udist = 1.0;
    float umass = 1.0;
    // Other parameters

    // Destructor
    virtual SamplingFunctionsSet get_sampler() const = 0;
    virtual ~ParticlesSetup() = default;

    // Other variable
    std::string SimulationTag = "";


private:
protected:
    template <typename T>
    void _write_toml_kvc(std::ofstream& fout, const std::string& key, const T& value, const std::string& comment,
                                int key_width = 14, int eq_pos = 18, int val_pos = 38, int comment_pos = 50) const {
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
};

class ParticlesSetupUniform : public ParticlesSetup{
public:
    // Range of x
    float xmin = -1.0f;
    float xmax = 1.0f;
    // Range of y
    float ymin = -1.0f;
    float ymax = 1.0f;
    // Range of z
    float zmin = -1.0f;
    float zmax = 1.0f;
    // Range of vx
    float vxmin = -1.0f;
    float vxmax = 1.0f;
    // Range of vy
    float vymin = -1.0f;
    float vymax = 1.0f;
    // Range of vz
    float vzmin = -1.0f;
    float vzmax = 1.0f;

    // Range of mass
    float mmin = 1e-20f;
    float mmax = 1.0f;
    

    // Constructor
    ParticlesSetupUniform(const std::string& simulation_tag){
        std::string filename = simulation_tag + ".setup";
        if (std::filesystem::exists(filename)) {
            std::cout << "Reading setup from " << filename << std::endl;
            _read_setupin_toml(filename);
        } else {
            std::cerr << "Setup file not found: " << filename << "\n"
          << "A default setup file has been created at the same location.\n"
          << "Please edit the file and re-run the program.\n";
            _make_setupin_toml(simulation_tag);
            std::exit(0); 
        }
    }
    
    // Other Method
    /*
        SamplingFunctionsSet get_sampler();
    Make the initial condition sampler that correspoinding to the mode of setup.
    */
    SamplingFunctionsSet get_sampler() const override {
        SamplingFunctionsSet sampler;
        std::mt19937 rng(std::random_device{}());

        std::uniform_real_distribution<float> distx(xmin, xmax);
        std::uniform_real_distribution<float> disty(ymin, ymax);
        std::uniform_real_distribution<float> distz(zmin, zmax);
        std::uniform_real_distribution<float> distvx(vxmin, vxmax);
        std::uniform_real_distribution<float> distvy(vymin, vymax);
        std::uniform_real_distribution<float> distvz(vzmin, vzmax);
        std::uniform_real_distribution<float> distm(mmin, mmax);
        
        sampler.coorsampler = [=]() mutable -> std::array<float, 6> {
            return {
                distx(rng),
                disty(rng),
                distz(rng),
                distvx(rng),
                distvy(rng),
                distvz(rng)
            };
        };    
        sampler.msampler = [=]() mutable { return distm(rng); };

        return sampler;
    }
    
        
private:
protected:  
    /*
        void _make_setupin_toml() const;
    (Internal Method) Extract a TOML setup file for setup
    ## Input
        - string filename: Name of input
    */
    void _make_setupin_toml(const std::string& filename) const;
    /*
        void _read_setupin_toml() const;
    (Internal Method) Read a TOML setup file for setup

     ## Input
        - string filename: Name of input
    */
   void _read_setupin_toml(const std::string& filename);
};

class ParticlesSetupIsotropic : public ParticlesSetup{
    public:
        // Center of sphere
        float rcx = 0.0f;
        float rcy = 0.0f;
        float rcz = 0.0f;
        // Range of r
        float rmax = 1.0f;
        float alpha = -2.0f; // r-sampling with power law rho(r) ∝ r^alpha

        // Range of vr
        float vrmin = -1.0f;
        float vrmax = 1.0f;
        // Range of vphi
        float vphimin = -1.0f;
        float vphimax = 1.0f;
        // Range of vtheta
        float vthetamin = -1.0f;
        float vthetamax = 1.0f;

    
        // Range of mass
        float mmin = 1e-20f;
        float mmax = 1.0f;
        
    
        // Constructor
        ParticlesSetupIsotropic(const std::string& simulation_tag){
            std::string filename = simulation_tag + ".setup";
            if (std::filesystem::exists(filename)) {
                std::cout << "Reading setup from " << filename << std::endl;
                _read_setupin_toml(filename);
            } else {
                std::cerr << "Setup file not found: " << filename << "\n"
              << "A default setup file has been created at the same location.\n"
              << "Please edit the file and re-run the program.\n";
                _make_setupin_toml(simulation_tag);
                std::exit(0); 
            }
        }
        
        // Other Method
        /*
            SamplingFunctionsSet get_sampler();
        Make the initial condition sampler that correspoinding to the mode of setup.

        # Output
        - `SamplingFunctionsSet`
        */
        SamplingFunctionsSet get_sampler() const override {
            std::mt19937 rng(std::random_device{}());

            SamplingFunctionsSet sampler;
            
            // velocity sampler
            std::uniform_real_distribution<float> distvr(vrmin, vrmax);
            std::uniform_real_distribution<float> distvphi(vphimin, vphimax);
            std::uniform_real_distribution<float> distvtheta(vthetamin, vthetamax);
            std::uniform_real_distribution<float> distm(mmin, mmax);

            // rsampler
            std::uniform_real_distribution<float> urand(0.0f,1.0f);
            auto r_sampler = [=]() mutable {
                float u = urand(rng);
                float a = -1.0f; 
                float exp = a + 3.0f;
                float norm = std::pow(rmax, exp) - std::pow(1e-3, exp);
                return std::pow(u * norm + std::pow(1e-3, exp), 1.0f / exp);
            };
            auto xyz_sampler = [=]() mutable {
                // sample position
                float r = r_sampler();
                float phi = 2.0f * M_PI * urand(rng);
                float theta = std::acos(1.0f - 2.0f * urand(rng));  // sampling on cos(theta)
                // sample velocity
                float vr = distvr(rng);
                float vphi = distvphi(rng);
                float vtheta = distvtheta(rng);

                std::array<float, 6> sph_coor{};
                sph_coor[0] = r;
                sph_coor[1] = phi;
                sph_coor[2] = theta;
                sph_coor[3] = vr;
                sph_coor[4] = vphi;
                sph_coor[5] = vtheta;

                std::array<float, 6> cart_coor = sph2cart(sph_coor);
                cart_coor[0] += rcx;                                // Shift to center
                cart_coor[1] += rcy;                                // Shift to center
                cart_coor[2] += rcz;                                // Shift to center

                return cart_coor;
            };
            sampler.coorsampler = xyz_sampler;
            sampler.msampler = [=]() mutable { return distm(rng); };
            return sampler;
        }
        
            
    private:
    protected:  
        /*
            void _make_setupin_toml() const;
        (Internal Method) Extract a TOML setup file for setup
        ## Input
            - string filename: Name of input
        */
        void _make_setupin_toml(const std::string& filename) const;
        /*
            void _read_setupin_toml() const;
        (Internal Method) Read a TOML setup file for setup
    
         ## Input
            - string filename: Name of input
        */
       void _read_setupin_toml(const std::string& filename);
    };