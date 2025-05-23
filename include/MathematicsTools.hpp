#pragma once
#include <math.h>
#include <vector>
#include <array>

/*
    std::array<float, 6> sph2cart(std::array<float, 6> sph_coor);

Transfering the general coordinate (r, phi, theta, vr, vphi, vtheta) to (x, y, z, vx, vy,vz)

## Input 
    - std::array<float, 6> sph_coor: The point in spherical coordinate
## Output 
    - std::array<float, 6> cart_coor: The point in cartisian coordinate
*/
std::array<float, 6> sph2cart(std::array<float, 6> sph_coor);

/*
    std::array<float, 6> cart2sph(std::array<float, 6> cart_coor);

Transfering the general coordinate (x, y, z, vx, vy,vz) to (r, phi, theta, vr, vphi, vtheta)

## Input 
    - std::array<float, 6> cart_coor: The point in cartisian coordinate
## Output 
    - std::array<float, 6> sph_coor: The point in spherical coordinate
*/
std::array<float, 6> cart2sph(std::array<float, 6> cart_coor);