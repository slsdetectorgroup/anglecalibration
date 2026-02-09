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

// TODO needs read_into_array
template <class CustomFile> struct badchannel_file_compatibility {
    static constexpr bool value =
        std::is_constructible<CustomFile, std::string>::value;
};

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

    /// @brief measured dead-time [s]
    double dead_time =
        76.08e-9; // measured dead-time 31.8E-9_DP, 95.E-9_DP other beam
                  // settings ! TODO: should it be configurable?

    /// @brief additional offset to sample detector offset (can change in
    /// experimental setup) [degrees]
    double offset = 0.0;

    /// @brief number of modules in detector
    size_t max_modules = 48; // TODO can be static

    std::vector<ssize_t> unconnected_modules{}; // list of unconnected modules

    ssize_t num_strips() const {
        return static_cast<ssize_t>(max_modules * strips_per_module);
    }
};

} // namespace angcal