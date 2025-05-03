#pragma once
#include <random>
#include <functional>

class SamplingFunctionsSet {
public:
    // Function set for sampling
    std::function<float()> xsampler, ysampler, zsampler;
    std::function<float()> vxsampler, vysampler, vzsampler;
    std::function<float()> msampler;

    // Destructor
    ~SamplingFunctionsSet() = default;
};