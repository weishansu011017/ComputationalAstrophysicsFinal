#pragma once
#include <math.h>
#include <vector>
#include <array>
#include <numeric>

/*
    float mean(const std::vector<T>& vec);

Calculate the arithmetic mean of a 1D container (e.g., std::vector).

## Input
    - const std::vector<T>& vec: The input 1D container containing numeric values.

## Output
    - float: The average (mean) value of the container. Returns 0.0 if the input is empty.
*/
template <typename T>
float mean(const std::vector<T>& vec) {
    if (vec.empty()) return 0.0f;
    return std::accumulate(vec.begin(), vec.end(), 0.0f) / static_cast<float>(vec.size());
}
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