#pragma once
#include <string.h>

/*
Particles setup properties
*/
class ParticlesSetup{
public:
    // Number of particles
    int N = 10000;
    
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
    
    // Other parameters

    // Destructor
    ~ParticlesSetup() = default;

    // Other variable
    std::string sampling_function = "Uniform";


private:
protected:    
};
