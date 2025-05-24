#include <math.h>
#include <vector>
#include <numeric>
#include "MathematicsTools.hpp"

std::array<float, 6> sph2cart(std::array<float, 6> sph_coor){
    // sph_coor = (r, phi, theta, vr, vphi, vtheta)
    // cart_coor = (x, y, z, vx, vy, vz)
    std::array<float, 6> cart_coor{};

    float r = sph_coor[0];
    float phi = sph_coor[1];
    float theta = sph_coor[2];
    float vr = sph_coor[3];
    float vphi = sph_coor[4];
    float vtheta = sph_coor[5];

    float cosphi = cosf(phi);
    float costheta = cosf(theta);
    float sinphi = sinf(phi);
    float sintheta = sinf(theta);

    // x
    cart_coor[0] = r * sintheta * cosphi;
    // y
    cart_coor[1] = r * sintheta * sinphi;
    // z
    cart_coor[2] = r * costheta;
    // vx
    cart_coor[3] = (sintheta * cosphi * vr) - (sinphi * vphi) + (costheta * cosphi * vtheta);
    // vy
    cart_coor[4] = (sintheta * sinphi * vr) + (cosphi * vphi) + (costheta * sinphi * vtheta);
    // vz
    cart_coor[5] = (costheta * vr) - (sintheta * vtheta);
    return cart_coor;
}

std::array<float, 6> cart2sph(std::array<float, 6> cart_coor){
    // sph_coor = (r, phi, theta, vr, vphi, vtheta)
    // cart_coor = (x, y, z, vx, vy, vz)
    std::array<float, 6> sph_coor{};

    float x = cart_coor[0];
    float y = cart_coor[1];
    float z = cart_coor[2];
    float vx = cart_coor[3];
    float vy = cart_coor[4];
    float vz = cart_coor[5];

    float s = sqrtf(x*x + y*y);
    float r = sqrtf(x*x + y*y + z*z);
    float phi = atan2f(y, x);
    if (phi < 0) phi += 2.0f * M_PI;
    float theta = acosf(z / r);
    
    float cosphi = cosf(phi);
    float costheta = cosf(theta);
    float sinphi = sinf(phi);
    float sintheta = sinf(theta);

    // r
    sph_coor[0] = r;
    // phi
    sph_coor[1] = phi;
    // theta
    sph_coor[2] = theta;
    // vr
    sph_coor[3] = (sintheta * cosphi * vx) + (sintheta * sinphi * vy) + (costheta * vz);
    // vphi
    sph_coor[4] = -(sinphi * vx) + (cosphi * vy);
    // vtheta
    sph_coor[5] = (costheta * cosphi * vx) + (costheta * sinphi * vy) - (sintheta * vz);

    return sph_coor;
}