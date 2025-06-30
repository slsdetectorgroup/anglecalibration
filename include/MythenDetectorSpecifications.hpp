#pragma once
#include <cmath>
#include <cstdint>

#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>

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

    MythenDetectorSpecifications(
        std::optional<std::shared_ptr<SimpleFileInterface>> custom_file_ptr =
            std::nullopt)
        : m_custom_file_ptr(std::move(custom_file_ptr)) {
        num_strips_ = max_modules_ * strips_per_module_;

        bad_channels =
            NDArray<bool, 1>(std::array<ssize_t, 1>{num_strips_}, false);
    }

    MythenDetectorSpecifications(
        const size_t max_modules, const double exposure_time,
        const double num_counters = 1, double bloffset = 1.532,
        std::optional<std::shared_ptr<SimpleFileInterface>> custom_file_ptr =
            std::nullopt)
        : max_modules_(max_modules), num_counters_(num_counters),
          exposure_time_(exposure_time), bloffset_(bloffset),
          m_custom_file_ptr(std::move(custom_file_ptr)) {
        num_strips_ = max_modules_ * strips_per_module_;

        bad_channels =
            NDArray<bool, 1>(std::array<ssize_t, 1>{num_strips_}, false);
    }

    void read_bad_channels_from_file(const std::string &filename) {
        if (m_custom_file_ptr.has_value()) {
            m_custom_file_ptr.value()->open(filename);
            m_custom_file_ptr.value()->read_into(
                reinterpret_cast<std::byte *>(bad_channels.data()));
        } else {
            throw std::runtime_error("provide ptr to CustomFile class");
        }
    }

    void
    set_unconnected_modules(const std::vector<ssize_t> &unconnected_modules) {
        m_unconnected_modules = NDArray<ssize_t, 1>(std::array<ssize_t, 1>{
            static_cast<ssize_t>(unconnected_modules.size())});

        std::copy(unconnected_modules.begin(), unconnected_modules.end(),
                  m_unconnected_modules.begin());

        // update bad_channels
        std::for_each(
            m_unconnected_modules.begin(), m_unconnected_modules.end(),
            [this](const auto &module_index) {
                for (size_t i = module_index * this->strips_per_module_;
                     i < (module_index + 1) * this->strips_per_module_; ++i)
                    this->bad_channels[i] = true;
            });
    }

    NDView<bool, 1> get_bad_channels() const { return bad_channels.view(); }

    void set_bad_channels(const NDArray<bool, 1> &&bad_channels_) {
        bad_channels = bad_channels_;
    }

    NDView<ssize_t, 1> get_unconnected_modules() const {
        return m_unconnected_modules.view();
    }

    static constexpr double pitch() { return pitch_; }

    static constexpr size_t strips_per_module() { return strips_per_module_; }

    size_t max_modules() const { return max_modules_; }

    size_t num_counters() const { return num_counters_; }

    double exposure_time() const { return exposure_time_; }

    double bloffset() const { return bloffset_; }

    double dtt0() const { return dtt0_; }

    static constexpr double min_angle() { return min_angle_; }

    static constexpr double max_angle() { return max_angle_; }

    ssize_t num_strips() const { return num_strips_; }

  private:
    static constexpr size_t strips_per_module_ = 1280;
    static constexpr double pitch_ = 0.05; // strip width [mm]
    static constexpr double min_angle_ =
        -180.0; // maybe shoudnt be static but configurable
    static constexpr double max_angle_ = 180.0;
    static constexpr double dtt0_ =
        0.0; // No idea what this is - probably configurable

    size_t max_modules_ = 48;

    size_t num_counters_ = 1;

    double exposure_time_ = 5.0; // TODO: could read from acquired file but
                                 // maybe should be configurable
    double bloffset_ = 1.532; // what is this? detector offset relative to what?

    ssize_t num_strips_{};

    NDArray<bool, 1> bad_channels{};
    NDArray<ssize_t, 1> m_unconnected_modules{}; // list of unconnected modules

    std::optional<std::shared_ptr<SimpleFileInterface>> m_custom_file_ptr{};
};

} // namespace angcal