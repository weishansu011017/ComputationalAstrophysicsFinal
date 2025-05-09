#pragma once
#include <random>
#include <functional>

class SamplingFunctionsSet {
public:
    // Function set for sampling
    std::function<std::array<float, 6>()> coorsampler;
    std::function<float()> msampler;

    // Destructor
    ~SamplingFunctionsSet() = default;
};