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
    float udist = 3.08567758128e18;
    float umass = 1.989e33;
    // Other parameters
    float softfactorx = 0.02;
    float tsfactor = 0.2;
    float bhTreeTheta = 0.2;

    // Destructor
    virtual ~ParticlesSetup() = default;

    // Other method (for decline only)
    virtual SamplingFunctionsSet get_sampler() const = 0;

    // Other variable
    int dimension = 3;
    std::string SimulationTag = "";
    float simulation_scale; 

private:
protected:
    

    template <typename T>
    void _write_toml_kvc(std::ofstream& fout, const std::string& key, const T& value, const std::string& comment,
                                int key_width = 14, int eq_pos = 18, int val_pos = 38, int comment_pos = 50) const; 
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

        // Mass sampler
        if (std::abs(mmax - mmin) < 1e-6f) {
            float fixed_mass = mmin;
            sampler.msampler = [=]() { return fixed_mass; };
        } else {
            std::uniform_real_distribution<float> distm(mmin, mmax);
            sampler.msampler = [=]() mutable { return distm(rng); };
        }


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
    float rmax = 5.0f;
    float alpha = -2.0f; // r-sampling with power law rho(r) ∝ r^alpha

    // Range of mass
    float mmin = 1e-3f;
    float mmax = 1e-2f;
    

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

        // rsampler
        std::uniform_real_distribution<float> urand(0.0f,1.0f);
        auto r_sampler = [=]() mutable {
            float u = urand(rng);
            float exp = alpha + 3.0f;
            float norm = std::pow(rmax, exp) - std::pow(1e-3, exp);
            return std::pow(u * norm + std::pow(1e-3, exp), 1.0f / exp);
        };

        // Estimate sigma from mass & radius: sigma^2 ~ G M / R
        float Mtot_estimate = 0.5f * (mmin + mmax) * static_cast<float>(N);
        float sigma2 = (Mtot_estimate) / (rmax * 0.5f);            // assume R_eff = rmax / 2
        float sigma = std::sqrt(sigma2);
        std::normal_distribution<float> gauss_v(0.0f, sigma);

        auto xyz_sampler = [=]() mutable {
            // sample position
            float r = r_sampler();
            float phi = 2.0f * M_PI * urand(rng);
            float theta = std::acos(1.0f - 2.0f * urand(rng));  // sampling on cos(theta)
            // sample velocity
            float vr = gauss_v(rng);
            float vphi = gauss_v(rng);
            float vtheta = gauss_v(rng);

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

        // Mass sampler
        if (std::abs(mmax - mmin) < 1e-6f) {
            float fixed_mass = mmin;
            sampler.msampler = [=]() { return fixed_mass; };
        } else {
            std::uniform_real_distribution<float> distm(mmin, mmax);
            sampler.msampler = [=]() mutable { return distm(rng); };
        }
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


class ParticlesSetupPlummer : public ParticlesSetup{
    public:
        // Center of sphere
        float rcx = 0.0f;
        float rcy = 0.0f;
        float rcz = 0.0f;
        // Range of r
        float rmax = 10.0f;
        float scalefactor = 0.1f; // Factor of scale length a = scalefactor * rmax
    
        // Other quantities
        float Mtot = 1.0e5f;    // Total mass
        
    
        // Constructor
        ParticlesSetupPlummer(const std::string& simulation_tag){
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

            // rsampler
            std::uniform_real_distribution<float> urand(0.0f, 1.0f);
            std::uniform_real_distribution<float> rtrial(0.0f, rmax);
            float a = scalefactor * rmax;
            auto r_sampler = [=]() mutable {
                float r;
                while (true) {
                    float u = std::clamp(urand(rng), 1e-6f, 1.0f - 1e-6f);
                    r = a / std::sqrt(std::pow(u, -2.0f / 3.0f) - 1.0f);
                    if (r <= rmax) return r;
                }
            };
            // vsampler
            std::uniform_real_distribution<float> urv(0.0f, 1.0f);
            auto v_sampler = [=](float r) mutable {
                float v_esc = std::sqrt(2.0f * Mtot / std::sqrt(r*r + a*a));
                
                while (true) {
                    float v_try = v_esc * urv(rng); 
                    float f = std::pow(1.0f - v_try * v_try / (v_esc * v_esc), 3.5f); 
                    float y = urv(rng);
                    if (y < f) return v_try;
                }
            };

            auto xyz_sampler = [=]() mutable {
                // sample position
                float r = r_sampler();
                float phi = 2.0f * M_PI * urand(rng);
                float theta = std::acos(1.0f - 2.0f * urand(rng));  // sampling on cos(theta)
                // sample velocity
                float v = v_sampler(r);

                std::array<float, 6> cart_coor{};
                if (dimension == 3){
                    cart_coor[0] = r * sinf(theta) * cosf(phi);
                    cart_coor[1] = r * sinf(theta) * sinf(phi);
                    cart_coor[2] = r * cosf(theta);
                    // Random sample solid angle
                    float costheta = 1.0f - 2.0f * urand(rng);
                    float sintheta = std::sqrt(1.0f - costheta * costheta);
                    float angle = 2.0f * M_PI * urand(rng);
                    cart_coor[3] = v * sintheta * cosf(angle);
                    cart_coor[4] = v * sintheta * sinf(angle);
                    cart_coor[5] = v * costheta;

                    cart_coor[0] += rcx;                                // Shift to center
                    cart_coor[1] += rcy;                                // Shift to center
                    cart_coor[2] += rcz;                                // Shift to center
                } else if (dimension == 2){
                    cart_coor[0] = r * cosf(phi);
                    cart_coor[1] = r * sinf(phi);
                    cart_coor[2] = 0.0f;
                    // Random sample angle
                    float angle = 2.0f * M_PI * urand(rng);
                    cart_coor[3] = v * cosf(angle);
                    cart_coor[4] = v * sinf(angle);
                    cart_coor[5] = 0.0f;

                    cart_coor[0] += rcx;                                // Shift to center
                    cart_coor[1] += rcy;                                // Shift to center
                }
                return cart_coor;
            };
            sampler.coorsampler = xyz_sampler;

            // Mass sampler
            float single_m = Mtot / N;
            sampler.msampler = [=]() { return single_m; };

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