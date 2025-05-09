#include <vector>
#include <math.h>
#include <random>
#include <iostream>
#include <string>
#include "MathematicsTools.hpp"


/*
    TEST: Testing the Coordinate transfering funcutionality

    Test content: 
    - Sph2Cart
    1. Generate v = (r, phi, theta, vr,vphi, vtheta)
    2. Test v - cart2sph(sph2cart(v)) is identical
    3. Repeat 10 times
    - Cart2Sph
    1. Generate v = (x, y, z, vx, vy, vz)
    2. Test v - sph2cart(cart2sph(v)) is identical
    3. Repeat 10 times
*/
int main(int argc, char** argv){
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> distr(0.0f, 10.0f);
    std::uniform_real_distribution<float> distphi(0.0f, 2*M_PI);
    std::uniform_real_distribution<float> disttheta(0.0f, M_PI);
    std::uniform_real_distribution<float> distvr(-10.0f, 10.0f);
    std::uniform_real_distribution<float> distvphi(-10.0f, 10.0f);
    std::uniform_real_distribution<float> distvtheta(-10.0f, 10.0f);
    std::uniform_real_distribution<float> distx(-10.0f, 10.0f);
    std::uniform_real_distribution<float> disty(-10.0f, 10.0f);
    std::uniform_real_distribution<float> distz(-10.0f, 10.0f);
    std::uniform_real_distribution<float> distvx(-10.0f, 10.0f);
    std::uniform_real_distribution<float> distvy(-10.0f, 10.0f);
    std::uniform_real_distribution<float> distvz(-10.0f, 10.0f);

    std::array<float, 6> cart_coor{};

    std::array<float, 6> sph_coor{};

    // Test Sph2Cart
    std::cout << "----------------------------------------\nTest spherical coordinate to cartisian coordinate.\n";
    for (int i = 0; i<11; i++){
        std::cout << "============\n";
        float r = distr(rng);
        float phi = distphi(rng);
        float theta = disttheta(rng);
        float vr = distvr(rng);
        float vphi = distvphi(rng);
        float vtheta = distvtheta(rng);
        sph_coor[0] = r;
        sph_coor[1] = phi;
        sph_coor[2] = theta;
        sph_coor[3] = vr;
        sph_coor[4] = vphi;
        sph_coor[5] = vtheta;

        std::cout << "Input: [ ";
        for (float val : sph_coor){
            std::cout << val << " ";
        }
        std::cout << "]" << std::endl;
        
        std::cout << "Output: [ ";
        for (float val : cart2sph(sph2cart(sph_coor))){
            std::cout << val << " ";
        }
        std::cout << "]" << std::endl;

        std::cout << "Diff: [ ";
        auto sph_out = cart2sph(sph2cart(sph_coor));
        for (int j = 0; j<6; j++){
            std::cout << sph_out[j] - sph_coor[j] << " ";
        };
        std::cout << "]" << std::endl;     
    };
    // Test Cart2Sph
    std::cout << "----------------------------------------\nTest cartisian coordinate to spherical coordinate.\n";
    for (int i = 0; i<11; i++){
        std::cout << "============\n";
        float x = distx(rng);
        float y = disty(rng);
        float z = distz(rng);
        float vx = distvx(rng);
        float vy = distvy(rng);
        float vz = distvz(rng);
        cart_coor[0] = x;
        cart_coor[1] = y;
        cart_coor[2] = z;
        cart_coor[3] = vx;
        cart_coor[4] = vy;
        cart_coor[5] = vz;

        std::cout << "Input: [ ";
        for (float val : cart_coor){
            std::cout << val << " ";
        }
        std::cout << "]" << std::endl;
        
        std::cout << "Output: [ ";
        for (float val : sph2cart(cart2sph(cart_coor))){
            std::cout << val << " ";
        }
        std::cout << "]" << std::endl;
        std::cout << "Diff: [ ";
        auto cart_out = sph2cart(cart2sph(cart_coor));
        for (int j = 0; j<6; j++){
            std::cout << cart_out[j] - cart_coor[j] << " ";
        };
        std::cout << "]" << std::endl;  
    };
    return 0;
}   