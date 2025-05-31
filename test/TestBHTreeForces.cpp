#include <iostream>
#include <random>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <vector>
#include <functional>
#include <algorithm>
#include <numeric>
#include <cmath>
#include "ParticlesTable.hpp"
#include "UnitsTable.hpp"

/*
    Complete Updated TEST: Barnes-Hut Tree Force Calculation with Robust Analysis
    
    Features:
    - Multiple timing runs with statistical analysis
    - Improved particle distributions for realistic forces
    - Comprehensive validation and diagnostics
    - Better error detection and reporting
    - Prevents compiler optimization artifacts
*/

struct TimingResult {
    long mean_time_us;
    long std_dev_us;
    long min_time_us;
    long max_time_us;
    int num_runs;
    
    void print(const std::string& label) const {
        std::cout << label << ": " << mean_time_us << " ± " << std_dev_us 
                  << " μs [" << min_time_us << "-" << max_time_us << "] (" << num_runs << " runs)\n";
    }
};

struct PerformanceData {
    int N;
    float theta;
    TimingResult direct_timing;
    TimingResult bh_timing;
    float speedup;
    float rms_error;
    float max_error;
    float avg_force_magnitude;
    std::string method;
    std::string system_type;
};

struct ParticleForceData {
    int N;
    std::vector<float> x, y, z;
    std::vector<float> m;
    std::vector<float> fx_direct, fy_direct, fz_direct;
    std::vector<float> fx_bh, fy_bh, fz_bh;
    std::vector<float> force_errors;
    float theta;
    std::string method;
    std::string system_type;
};

// Prevent compiler optimization
volatile float force_sink = 0.0f;

void preventOptimization(const ParticlesTable& pt) {
    float sum = 0.0f;
    for (int i = 0; i < pt.N; i++) {
        sum += pt._ax[i] + pt._ay[i] + pt._az[i];
    }
    force_sink = sum;
}

TimingResult benchmarkFunction(std::function<void()> func, int num_runs = 5) {
    std::vector<long> times;
    
    // Warm up run
    func();
    
    for (int i = 0; i < num_runs; i++) {
        auto start = std::chrono::high_resolution_clock::now();
        func();
        auto end = std::chrono::high_resolution_clock::now();
        
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        times.push_back(duration.count());
    }
    
    // Calculate statistics
    TimingResult result;
    result.num_runs = num_runs;
    result.mean_time_us = std::accumulate(times.begin(), times.end(), 0L) / num_runs;
    result.min_time_us = *std::min_element(times.begin(), times.end());
    result.max_time_us = *std::max_element(times.begin(), times.end());
    
    // Standard deviation
    double variance = 0.0;
    for (long time : times) {
        variance += (time - result.mean_time_us) * (time - result.mean_time_us);
    }
    result.std_dev_us = static_cast<long>(std::sqrt(variance / num_runs));
    
    return result;
}

void validateParticleSystem(ParticlesTable& pt, const std::string& system_name) {
    std::cout << "\n--- Validating " << system_name << " (N=" << pt.N << ") ---\n";
    
    // Check particle distribution
    float x_min = *std::min_element(pt.x.begin(), pt.x.end());
    float x_max = *std::max_element(pt.x.begin(), pt.x.end());
    float y_min = *std::min_element(pt.y.begin(), pt.y.end());
    float y_max = *std::max_element(pt.y.begin(), pt.y.end());
    float z_min = *std::min_element(pt.z.begin(), pt.z.end());
    float z_max = *std::max_element(pt.z.begin(), pt.z.end());
    
    std::cout << "Spatial extent: x[" << std::fixed << std::setprecision(3) 
              << x_min << ", " << x_max << "], y[" << y_min << ", " << y_max 
              << "], z[" << z_min << ", " << z_max << "]\n";
    
    // Check minimum separation
    float min_separation = 1e10f;
    int close_pairs = 0;
    for (int i = 0; i < std::min(pt.N, 100); i++) {  // Sample to avoid O(N²) for large N
        for (int j = i + 1; j < std::min(pt.N, 100); j++) {
            float dx = pt.x[i] - pt.x[j];
            float dy = pt.y[i] - pt.y[j];
            float dz = pt.z[i] - pt.z[j];
            float r = std::sqrt(dx*dx + dy*dy + dz*dz);
            min_separation = std::min(min_separation, r);
            if (r < pt.h[i]) close_pairs++;
        }
    }
    
    std::cout << "Min separation (sample): " << min_separation << ", Softening: " << pt.h[0] << "\n";
    if (close_pairs > 0) {
        std::cout << "WARNING: " << close_pairs << " particle pairs closer than softening length\n";
    }
    
    // Calculate test force to check magnitudes
    pt.calculate_a_dirnbody();
    preventOptimization(pt);
    
    float total_force = 0.0f, max_force = 0.0f;
    int valid_forces = 0;
    for (int i = 0; i < pt.N; i++) {
        float f = std::sqrt(pt._ax[i]*pt._ax[i] + pt._ay[i]*pt._ay[i] + pt._az[i]*pt._az[i]);
        if (std::isfinite(f)) {
            total_force += f;
            max_force = std::max(max_force, f);
            valid_forces++;
        }
    }
    float avg_force = (valid_forces > 0) ? total_force / valid_forces : 0.0f;
    
    std::cout << "Force stats: avg=" << std::scientific << std::setprecision(2) 
              << avg_force << ", max=" << max_force << ", valid=" << valid_forces << "/" << pt.N << "\n";
    
    if (avg_force < 1e-8f) {
        std::cout << "WARNING: Forces very small - particles may be too far apart!\n";
    }
    if (max_force > 1e4f) {
        std::cout << "WARNING: Forces very large - numerical instability possible!\n";
    }
    if (valid_forces < pt.N) {
        std::cout << "ERROR: " << (pt.N - valid_forces) << " particles have invalid forces!\n";
    }
    
    // Check if 2D
    float z_range = z_max - z_min;
    float xy_range = std::max(x_max - x_min, y_max - y_min);
    bool is_2d = pt.is2D();
    std::cout << "Z-range/XY-range ratio: " << (z_range / xy_range) 
              << " -> Detected as: " << (is_2d ? "2D" : "3D") << "\n";
}

ParticlesTable createTestSystem(int N, const std::string& type, int seed = 42) {
    UnitsTable units(1.0f, 1.0f);
    ParticlesTable pt(units, N);
    pt.SimulationTag = "bh_test_" + type;
    
    std::mt19937 gen(seed);
    std::uniform_real_distribution<float> uniform(0.0f, 1.0f);
    std::normal_distribution<float> normal(0.0f, 1.0f);
    
    float totalMass = 0.0f;
    
    if (type == "plummer") {
        // Compact Plummer sphere - realistic astrophysical distribution
        float plummer_a = 0.5f;  // Scale radius
        float virial_velocity = std::sqrt(1.0f / plummer_a);  // Rough virial velocity
        
        for (int i = 0; i < N; i++) {
            // Sample radius using inverse transform
            float u = uniform(gen);
            float r = plummer_a / std::sqrt(std::pow(u, -2.0f/3.0f) - 1.0f);
            r = std::min(r, 10.0f * plummer_a);  // Limit outliers
            
            // Isotropic angles
            float theta = std::acos(1.0f - 2.0f * uniform(gen));
            float phi = 2.0f * M_PI * uniform(gen);
            
            pt.x[i] = r * std::sin(theta) * std::cos(phi);
            pt.y[i] = r * std::sin(theta) * std::sin(phi);
            pt.z[i] = r * std::cos(theta);
            
            // Velocities from virial relation
            pt.vx[i] = normal(gen) * virial_velocity * 0.3f;
            pt.vy[i] = normal(gen) * virial_velocity * 0.3f;
            pt.vz[i] = normal(gen) * virial_velocity * 0.3f;
            
            pt.m[i] = 1.0f / N;
            totalMass += pt.m[i];
            pt.h[i] = 0.02f;  // Softening ~4% of scale radius
            pt.dt[i] = 0.001f;
        }
    } else if (type == "random_3d") {
        // Uniform sphere with reasonable density
        float sphere_radius = 1.0f;
        for (int i = 0; i < N; i++) {
            // Uniform distribution in sphere
            float u = uniform(gen);
            float r = sphere_radius * std::pow(u, 1.0f/3.0f);
            float theta = std::acos(1.0f - 2.0f * uniform(gen));
            float phi = 2.0f * M_PI * uniform(gen);
            
            pt.x[i] = r * std::sin(theta) * std::cos(phi);
            pt.y[i] = r * std::sin(theta) * std::sin(phi);
            pt.z[i] = r * std::cos(theta);
            
            // Random velocities
            pt.vx[i] = normal(gen) * 0.1f;
            pt.vy[i] = normal(gen) * 0.1f;
            pt.vz[i] = normal(gen) * 0.1f;
            
            pt.m[i] = 1.0f / N;
            totalMass += pt.m[i];
            pt.h[i] = 0.02f;
            pt.dt[i] = 0.001f;
        }
    } else if (type == "disk_2d") {
        // Exponential disk profile
        float disk_scale = 1.0f;
        for (int i = 0; i < N; i++) {
            // Exponential disk profile: P(r) ∝ r * exp(-r/R)
            float u1 = uniform(gen);
            float u2 = uniform(gen);
            float r = -disk_scale * std::log(u1 * u2);  // Approximate sampling
            r = std::min(r, 5.0f * disk_scale);  // Limit extent
            
            float phi = 2.0f * M_PI * uniform(gen);
            pt.x[i] = r * std::cos(phi);
            pt.y[i] = r * std::sin(phi);
            pt.z[i] = normal(gen) * 0.01f;  // Very thin disk
            
            // Circular velocity with dispersion
            float v_circ = std::sqrt(0.5f * r);  // Approximate for disk
            pt.vx[i] = -v_circ * std::sin(phi) + normal(gen) * 0.1f;
            pt.vy[i] = v_circ * std::cos(phi) + normal(gen) * 0.1f;
            pt.vz[i] = normal(gen) * 0.01f;
            
            pt.m[i] = 1.0f / N;
            totalMass += pt.m[i];
            pt.h[i] = 0.02f;
            pt.dt[i] = 0.001f;
        }
    }
    
    pt.Mtot = totalMass;
    return pt;
}

float calculateForceAccuracy(const ParticlesTable& pt, 
                           const std::vector<float>& ax_direct, 
                           const std::vector<float>& ay_direct, 
                           const std::vector<float>& az_direct,
                           std::vector<float>& errors) {
    errors.resize(pt.N);
    float rms_error = 0.0f;
    int valid_comparisons = 0;
    
    for (int i = 0; i < pt.N; i++) {
        // Calculate force magnitudes
        float f_direct = std::sqrt(ax_direct[i]*ax_direct[i] + ay_direct[i]*ay_direct[i] + az_direct[i]*az_direct[i]);
        float f_bh = std::sqrt(pt._ax[i]*pt._ax[i] + pt._ay[i]*pt._ay[i] + pt._az[i]*pt._az[i]);
        
        if (f_direct > 1e-10f && std::isfinite(f_direct) && std::isfinite(f_bh)) {
            // Relative error in force magnitude
            errors[i] = std::abs(f_bh - f_direct) / f_direct;
            rms_error += errors[i] * errors[i];
            valid_comparisons++;
        } else {
            errors[i] = 0.0f;  // Skip negligible forces
        }
    }
    
    return (valid_comparisons > 0) ? std::sqrt(rms_error / valid_comparisons) : 0.0f;
}

void savePerformanceData(const std::vector<PerformanceData>& data, const std::string& filename) {
    std::ofstream file(filename);
    file << "N,theta,direct_mean_us,direct_std_us,bh_mean_us,bh_std_us,speedup,rms_error,max_error,avg_force,method,system_type\n";
    
    for (const auto& d : data) {
        file << d.N << "," << d.theta << "," 
             << d.direct_timing.mean_time_us << "," << d.direct_timing.std_dev_us << ","
             << d.bh_timing.mean_time_us << "," << d.bh_timing.std_dev_us << ","
             << d.speedup << "," << d.rms_error << "," << d.max_error << "," 
             << d.avg_force_magnitude << "," << d.method << "," << d.system_type << "\n";
    }
    file.close();
}

void saveParticleForceData(const ParticleForceData& data, const std::string& filename) {
    std::ofstream file(filename);
    file << "x,y,z,m,fx_direct,fy_direct,fz_direct,fx_bh,fy_bh,fz_bh,force_error,theta,method,system_type\n";
    
    for (int i = 0; i < data.N; i++) {
        file << data.x[i] << "," << data.y[i] << "," << data.z[i] << "," << data.m[i] << ","
             << data.fx_direct[i] << "," << data.fy_direct[i] << "," << data.fz_direct[i] << ","
             << data.fx_bh[i] << "," << data.fy_bh[i] << "," << data.fz_bh[i] << ","
             << data.force_errors[i] << "," << data.theta << "," << data.method 
             << "," << data.system_type << "\n";
    }
    file.close();
}

int main() {
    std::cout << "\n=== Robust Barnes-Hut Performance & Accuracy Analysis ===\n\n";
    
    std::vector<PerformanceData> performance_results;
    std::vector<ParticleForceData> force_comparisons;
    
    // Test parameters
    std::vector<int> particle_counts = {100, 300, 500, 1000, 2000, 5000};
    std::vector<float> thetas = {0.3f, 0.5f, 0.7f, 1.0f};
    std::vector<std::string> system_types = {"plummer", "random_3d", "disk_2d"};
    
    int timing_runs = 3;  // Number of timing runs for statistics
    
    for (const std::string& sys_type : system_types) {
        std::cout << "\n" << std::string(60, '=') << "\n";
        std::cout << "Testing " << sys_type << " system\n";
        std::cout << std::string(60, '=') << "\n";
        
        for (int N : particle_counts) {
            std::cout << "\n--- N = " << N << " particles ---\n";
            
            // Create and validate test system
            ParticlesTable pt = createTestSystem(N, sys_type);
            validateParticleSystem(pt, sys_type);
            
            // Direct N-body timing
            auto direct_func = [&pt]() { 
                pt.calculate_a_dirnbody(); 
                preventOptimization(pt);
            };
            TimingResult direct_timing = benchmarkFunction(direct_func, timing_runs);
            direct_timing.print("Direct N-body");
            
            // Store reference forces
            std::vector<float> ax_direct(N), ay_direct(N), az_direct(N);
            for (int i = 0; i < N; i++) {
                ax_direct[i] = pt._ax[i];
                ay_direct[i] = pt._ay[i];
                az_direct[i] = pt._az[i];
            }
            
            // Calculate average force magnitude for reference
            float total_force = 0.0f;
            for (int i = 0; i < N; i++) {
                total_force += std::sqrt(ax_direct[i]*ax_direct[i] + ay_direct[i]*ay_direct[i] + az_direct[i]*az_direct[i]);
            }
            float avg_force_magnitude = total_force / N;
            
            std::cout << "\nBarnes-Hut Results:\n";
            std::cout << std::setw(8) << "Theta" 
                      << std::setw(18) << "Time (μs)"
                      << std::setw(12) << "Speedup"
                      << std::setw(15) << "RMS Error"
                      << std::setw(15) << "Max Error" << "\n";
            std::cout << std::string(68, '-') << "\n";
            
            for (float theta : thetas) {
                pt.setBHTreeTheta(theta);
                
                // Barnes-Hut timing
                auto bh_func = [&pt]() { 
                    pt.calculate_a_BHtree(); 
                    preventOptimization(pt);
                };
                TimingResult bh_timing = benchmarkFunction(bh_func, timing_runs);
                
                // Calculate accuracy
                std::vector<float> errors;
                float rms_error = calculateForceAccuracy(pt, ax_direct, ay_direct, az_direct, errors);
                float max_error = *std::max_element(errors.begin(), errors.end());
                
                float speedup = static_cast<float>(direct_timing.mean_time_us) / bh_timing.mean_time_us;
                
                std::cout << std::fixed << std::setprecision(2);
                std::cout << std::setw(8) << theta
                          << std::setw(12) << bh_timing.mean_time_us << " ± " << std::setw(3) << bh_timing.std_dev_us
                          << std::setw(12) << speedup << "x"
                          << std::setw(15) << std::scientific << rms_error
                          << std::setw(15) << max_error << "\n";
                
                // Store performance data
                PerformanceData perf;
                perf.N = N;
                perf.theta = theta;
                perf.direct_timing = direct_timing;
                perf.bh_timing = bh_timing;
                perf.speedup = speedup;
                perf.rms_error = rms_error;
                perf.max_error = max_error;
                perf.avg_force_magnitude = avg_force_magnitude;
                perf.method = pt.is2D() ? "2D" : "3D";
                perf.system_type = sys_type;
                performance_results.push_back(perf);
                
                // Store detailed force comparison for selected cases
                if ((N == 500 && theta == 0.5f) || (N == 2000 && theta == 0.5f)) {
                    ParticleForceData force_data;
                    force_data.N = N;
                    force_data.theta = theta;
                    force_data.method = pt.is2D() ? "2D" : "3D";
                    force_data.system_type = sys_type;
                    
                    force_data.x = pt.x;
                    force_data.y = pt.y;
                    force_data.z = pt.z;
                    force_data.m = pt.m;
                    
                    force_data.fx_direct = ax_direct;
                    force_data.fy_direct = ay_direct;
                    force_data.fz_direct = az_direct;
                    
                    force_data.fx_bh.resize(N);
                    force_data.fy_bh.resize(N);
                    force_data.fz_bh.resize(N);
                    for (int i = 0; i < N; i++) {
                        force_data.fx_bh[i] = pt._ax[i];
                        force_data.fy_bh[i] = pt._ay[i];
                        force_data.fz_bh[i] = pt._az[i];
                    }
                    
                    force_data.force_errors = errors;
                    force_comparisons.push_back(force_data);
                }
            }
        }
    }
    
    // Extended scaling test with larger N
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Extended Scaling Analysis\n";
    std::cout << std::string(60, '=') << "\n";
    
    std::vector<int> large_N = {1000, 3000, 5000, 10000, 20000};
    std::ofstream scaling_file("scaling_analysis.csv");
    scaling_file << "N,direct_mean_us,direct_std_us,bh_mean_us,bh_std_us,speedup,avg_force,system_type\n";
    
    for (int N : large_N) {
        std::cout << "\nTesting N = " << N << "...\n";
        
        ParticlesTable pt = createTestSystem(N, "plummer", 42);
        pt.setBHTreeTheta(0.5f);
        
        // Validate system
        pt.calculate_a_dirnbody();
        float total_force = 0.0f;
        for (int i = 0; i < N; i++) {
            total_force += std::sqrt(pt._ax[i]*pt._ax[i] + pt._ay[i]*pt._ay[i] + pt._az[i]*pt._az[i]);
        }
        float avg_force = total_force / N;
        
        // Time both methods
        auto direct_func = [&pt]() { pt.calculate_a_dirnbody(); preventOptimization(pt); };
        auto bh_func = [&pt]() { pt.calculate_a_BHtree(); preventOptimization(pt); };
        
        TimingResult direct_timing = benchmarkFunction(direct_func, 3);
        TimingResult bh_timing = benchmarkFunction(bh_func, 3);
        
        float speedup = static_cast<float>(direct_timing.mean_time_us) / bh_timing.mean_time_us;
        
        std::cout << "N=" << N << ": ";
        direct_timing.print("Direct");
        std::cout << "       ";
        bh_timing.print("Barnes-Hut");
        std::cout << "       Speedup: " << std::fixed << std::setprecision(2) << speedup << "x\n";
        
        scaling_file << N << "," << direct_timing.mean_time_us << "," << direct_timing.std_dev_us << ","
                    << bh_timing.mean_time_us << "," << bh_timing.std_dev_us << ","
                    << speedup << "," << avg_force << ",plummer\n";
    }
    scaling_file.close();
    
    // Save all analysis data
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Saving Analysis Data\n";
    std::cout << std::string(60, '=') << "\n";
    
    savePerformanceData(performance_results, "robust_bh_performance.csv");
    
    for (size_t i = 0; i < force_comparisons.size(); i++) {
        std::string filename = "robust_force_comparison_" + std::to_string(i) + ".csv";
        saveParticleForceData(force_comparisons[i], filename);
    }
    
    // Generate summary statistics
    std::ofstream summary("analysis_summary.txt");
    summary << "Barnes-Hut Analysis Summary\n";
    summary << "==========================\n\n";
    
    // Find best speedups
    auto best_speedup = *std::max_element(performance_results.begin(), performance_results.end(),
        [](const PerformanceData& a, const PerformanceData& b) { return a.speedup < b.speedup; });
    
    auto best_accuracy = *std::min_element(performance_results.begin(), performance_results.end(),
        [](const PerformanceData& a, const PerformanceData& b) { return a.rms_error < b.rms_error; });
    
    summary << "Best speedup: " << best_speedup.speedup << "x at N=" << best_speedup.N 
            << ", θ=" << best_speedup.theta << " (" << best_speedup.system_type << ", " << best_speedup.method << ")\n";
    summary << "Best accuracy: " << best_accuracy.rms_error << " at θ=" << best_accuracy.theta 
            << ", N=" << best_accuracy.N << " (" << best_accuracy.system_type << ")\n";
    
    // Count speedup achievements
    int speedup_count = std::count_if(performance_results.begin(), performance_results.end(),
        [](const PerformanceData& d) { return d.speedup > 1.0f; });
    
    summary << "\nSpeedup achieved in " << speedup_count << "/" << performance_results.size() 
            << " test cases (" << (100.0f * speedup_count / performance_results.size()) << "%)\n";
    
    summary.close();
    
    std::cout << "\nFiles created:\n";
    std::cout << "• robust_bh_performance.csv - Complete performance data\n";
    std::cout << "• robust_force_comparison_*.csv - Detailed force accuracy\n";
    std::cout << "• scaling_analysis.csv - Extended scaling analysis\n";
    std::cout << "• analysis_summary.txt - Summary statistics\n";
    std::cout << "\nAnalysis complete! Use Python script for visualization.\n\n";
    
    return 0;
}