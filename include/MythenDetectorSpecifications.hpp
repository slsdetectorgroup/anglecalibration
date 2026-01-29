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

class MythenDetectorSpecifications {

  public:
    // TODO: constructor that reads from a config file

    /**
     * @brief constructor for MythenDetectorSpecifications
     * @param custom_file_ptr (optional) pass Filereader to read bad channels
     * file default
     */
    MythenDetectorSpecifications() = default;

    /**
     * @brief constructor for MythenDetectorSpecifications
     * @param offset Additional offset to sample detector offset.
     * @param num_counters Number of counters active (default 1)
     * @param max_modules Number of modules in detector (default 48).
     */
    MythenDetectorSpecifications(const double offset,
                                 const size_t num_counters = 1,
                                 const size_t max_modules = 48)
        : offset_(offset), num_counters_(num_counters),
          max_modules_(max_modules) {}

    void
    set_unconnected_modules(const std::vector<ssize_t> &unconnected_modules) {
        m_unconnected_modules = unconnected_modules;
    }

    std::vector<ssize_t> &get_unconnected_modules() {
        return m_unconnected_modules;
    }

    static constexpr double pitch() { return pitch_; }

    static constexpr double transverse_width() { return transverse_width_; }

    static constexpr size_t strips_per_module() { return strips_per_module_; }

    /**
     * @brief number of modules in detector
     * (default '48')
     */
    size_t max_modules() const { return max_modules_; }

    size_t num_counters() const { return num_counters_; }

    static constexpr double sample_detector_offset() {
        return sample_detector_offset_;
    }

    double offset() const { return offset_; }

    static constexpr double min_angle() { return min_angle_; }

    static constexpr double max_angle() { return max_angle_; }

    /**
     * @brief total number of strips/channels in detector
     */
    ssize_t num_strips() const { return max_modules_ * strips_per_module_; }

  private:
    /// @brief number of strips/channels per module
    static constexpr size_t strips_per_module_ = 1280;

    /// @brief Strip/channel width of Mythen detector [mm]
    static constexpr double pitch_ = 0.05;

    /// @brief Strip/channel height of Mythen detector [mm]
    static constexpr double transverse_width_ = 8.0;

    /// @brief Minimum potential detector angle
    /// (measured as displacement of first strip) [degrees]
    static constexpr double min_angle_ =
        -180.0; // maybe shoudnt be static but configurable

    /// @brief Maximum potential detector angle
    /// (measured as displacement of first strip) [degrees]
    static constexpr double max_angle_ = 180.0;

    /// @brief Offset between sample horizontal plane and detector [degrees]
    static constexpr double sample_detector_offset_ = 1.4715;

    /// @brief additional offset to sample detector offset (can change in
    /// experimental setup) [degrees]
    double offset_ = 0.0;

    /// @brief num counters active in detector
    size_t num_counters_ = 1;

    /// @brief number of modules in detector
    size_t max_modules_ = 48; // TODO: could be static if bloffset is static

    std::vector<ssize_t> m_unconnected_modules{}; // list of unconnected modules
};

} // namespace angcal