#pragma once
#include <random>
#include <functional>
#include "ParticlesSetup.hpp"

class SamplingFunctionsSet {
public:
    // Function set for sampling
    std::function<float()> xsampler, ysampler, zsampler;
    std::function<float()> vxsampler, vysampler, vzsampler;
    std::function<float()> msampler;

    // Constructor
    SamplingFunctionsSet(const ParticlesSetup& setup)
    {
        std::mt19937 rng(std::random_device{}());

        if (setup.sampling_function == "Uniform") {
            std::uniform_real_distribution<float> distx(setup.xmin, setup.xmax);
            std::uniform_real_distribution<float> disty(setup.ymin, setup.ymax);
            std::uniform_real_distribution<float> distz(setup.zmin, setup.zmax);
            std::uniform_real_distribution<float> distvx(setup.vxmin, setup.vxmax);
            std::uniform_real_distribution<float> distvy(setup.vymin, setup.vymax);
            std::uniform_real_distribution<float> distvz(setup.vzmin, setup.vzmax);
            std::uniform_real_distribution<float> distm(setup.mmin, setup.mmax);

            xsampler = [=]() mutable { return distx(rng); };
            ysampler = [=]() mutable { return disty(rng); };
            zsampler = [=]() mutable { return distz(rng); };

            vxsampler = [=]() mutable { return distvx(rng); };
            vysampler = [=]() mutable { return distvy(rng); };
            vzsampler = [=]() mutable { return distvz(rng); };

            msampler = [=]() mutable { return distm(rng); };
        } else {
            throw std::runtime_error("Unsupported sampling function: " + setup.sampling_function);
        }
    }

    // Destructor
    ~SamplingFunctionsSet() = default;

protected:
    // 可擴充其他策略、rng 或參數
};