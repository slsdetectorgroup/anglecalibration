#pragma once
#include <cmath>
#include <cstdint>

#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>

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
    MythenDetectorSpecifications() {

        bad_channels =
            NDArray<bool, 1>(std::array<ssize_t, 1>{num_strips()}, false);
    }

    // TODO: is this constructor still relevant? or does one read from the file?
    /**
     * @brief constructor for MythenDetectorSpecifications
     * @param max_modules Number of modules in detector (default 48).
     * @param exposure_time Exposure time [s].
     * @param num_counters Number of counters active.
     * @param bloffset
     */
    MythenDetectorSpecifications(const size_t max_modules,
                                 const double exposure_time,
                                 const double num_counters = 1,
                                 double bloffset = 1.4715)
        : max_modules_(max_modules), num_counters_(num_counters),
          exposure_time_(exposure_time), bloffset_(bloffset) {

        bad_channels =
            NDArray<bool, 1>(std::array<ssize_t, 1>{num_strips()}, false);
    }

    /**
     * @brief read bad channels from file
     * @param filename bad channels filename
     * @param file_reader file_reader to read bad channels file (default:
     * CustomBadChannelsFile bad channel file is expected to be a text file
     * where each line stores the channel index of a bad channel. Consecutive
     * bad channels can be stored in one line by seperating the first and last
     * channel index of the bad channel block e.g.
     * bad_channel_index0-bad_channel_index1.))
     */
    void read_bad_channels_from_file(
        const std::string &filename,
        std::shared_ptr<SimpleFileInterface> file_reader =
            std::make_shared<CustomBadChannelsFile>()) {

        file_reader->open(filename);
        file_reader->read_into(
            reinterpret_cast<std::byte *>(bad_channels.data()));
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

    /*
    NDArray<bool, 1> get_bad_channels() const { return bad_channels; }
    */

    NDView<bool, 1> get_bad_channels() const { return bad_channels.view(); }

    /**
     * @brief set bad channels
     * @param bad_channels_ boolean NDArray of size of total number of
     * strips/channels in detector, stores 'TRUE' if strip is a bad channel
     * otherwise 'FALSE'
     */
    void set_bad_channels(const NDArray<bool, 1> &bad_channels_) {
        bad_channels = bad_channels_;
    }

    NDArray<ssize_t, 1> get_unconnected_modules() const {
        return m_unconnected_modules;
    }

    static constexpr double pitch() { return pitch_; }

    static constexpr size_t strips_per_module() { return strips_per_module_; }

    /**
     * @brief number of modules in detector
     * (default '48')
     */
    size_t max_modules() const { return max_modules_; }

    size_t num_counters() const { return num_counters_; }

    double exposure_time() const { return exposure_time_; }

    double bloffset() const { return bloffset_; }

    double dtt0() const { return dtt0_; }

    static constexpr double min_angle() { return min_angle_; }

    static constexpr double max_angle() { return max_angle_; }

    /**
     * @brief total number of strips/channels in detector
     */
    ssize_t num_strips() const { return max_modules_ * strips_per_module_; }

  private:
    /**
     * number of strips/channels per module
     */
    static constexpr size_t strips_per_module_ = 1280;
    /**
     * Strip/channel width of Mythen detector [mm]
     */
    static constexpr double pitch_ = 0.05;
    /** Minimum potential detector angle
     *(measured as displacement of first strip) [degree]
     */
    static constexpr double min_angle_ =
        -180.0; // maybe shoudnt be static but configurable

    /** Maximum potential detector angle
     *(measured as displacement of first strip) [degree]
     */
    static constexpr double max_angle_ = 180.0;
    static constexpr double dtt0_ = 0.0;

    /**
     * number of modules in detector
     */
    size_t max_modules_ = 48;

    /**
     * num counters active in detector
     */
    size_t num_counters_ = 1;

    /**
     * exposure time [s]
     */
    double exposure_time_ = 5.0; // TODO: could read from acquired file but
                                 // maybe should be configurable
    double bloffset_ =
        1.4715; // what is this? detector offset relative to what?

    /** Array of size strips/channels in detector, stores 'TRUE' if strip is a
     * bad channel otherwise 'FALSE'
     */
    NDArray<bool, 1> bad_channels{};
    NDArray<ssize_t, 1> m_unconnected_modules{}; // list of unconnected modules
};

} // namespace angcal