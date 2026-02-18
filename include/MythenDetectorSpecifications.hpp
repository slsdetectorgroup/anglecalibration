#pragma once
#include <cmath>
#include <cstdint>

#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "CustomFiles.hpp"
#include "aare/NDArray.hpp"
#include "helpers/FileInterface.hpp"

using namespace aare;
namespace angcal {

struct MythenDetectorSpecifications {

    /// @brief number of strips/channels per module
    static constexpr size_t strips_per_module = 1280;

    /// @brief Strip/channel width of Mythen detector [mm]
    static constexpr double pitch = 0.05;

    /// @brief Strip/channel height of Mythen detector [mm]
    static constexpr double transverse_width = 8.0;

    /// @brief average euclidean distance between sample and pixel [mm]
    double average_distance_sample_pixel =
        2500.0 / M_PI; // TODO why are two values
                       // 4420.97064144153710469121564923651006_DP (R_std_H) in
                       // Antonios code

    /// @brief Offset between sample horizontal plane and detector [degrees]
    double sample_detector_offset = 1.4715;

    /// @brief additional offset to sample detector offset (can change in
    /// experimental setup) [degrees]
    double offset = 0.0;

    /// @brief elastic correction factor
    double elastic_correction_factor = 0.0;

    /// @brief offste of detector to vertical axis (used for elastic correction)
    /// [degrees]
    double detector_vertical_axis_offset = 0.0;

    /// @brief measured dead-time [s]
    double dead_time =
        76.08e-9; // measured dead-time 31.8E-9_DP, 95.E-9_DP other beam
                  // settings ! TODO: should it be configurable?

    /// @brief number of modules in detector
    size_t max_modules = 48; // TODO can be static

    /// @brief list of module indices which are unconnected
    std::vector<ssize_t> unconnected_modules{};

    ssize_t num_strips() const {
        return static_cast<ssize_t>(max_modules * strips_per_module);
    }
};

} // namespace angcal