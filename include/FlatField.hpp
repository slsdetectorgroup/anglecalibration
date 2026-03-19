
/**
 * stores flatfield for angle calibration
 */

#pragma once
#include <cmath>
#include <cstdint>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>

#include "CustomFiles.hpp"
#include "MythenDetectorSpecifications.hpp"
#include "MythenFileReader.hpp"
#include "Parameters.hpp"
#include "aare/File.hpp"
#include "aare/NDArray.hpp"
#include "helpers/FileInterface.hpp"
#include "helpers/LogFiles.hpp"

using namespace aare;
namespace angcal {

class FlatField {

  public:
    /**
     * @tparam custom_file_ptr: Fileclass that inherits from
     * DetectorFileInterface class needs to overload read_frame and constructor
     * needs to have the signature Fileclass(std::string filename, ssize_t rows,
     * ssize_t cols, std::string reading_mode) rows, and cols define the
     * intended image size to read
     */
    FlatField(std::shared_ptr<MythenDetectorSpecifications> mythen_detector_);

    void set_scale_factor(const double scale_factor_);

    double get_scale_factor() const;

    void set_soft_window(const std::pair<double, double> soft_window_);

    std::pair<double, double> get_soft_window() const;

    double elastic_correction(const double detector_angle) const;

    double diffraction_angle_from_DG_parameters(
        const size_t module_index, const double detector_angle,
        size_t strip_index, const double distance_to_strip) const;

    double transverse_width_correction_factor(const size_t module_index,
                                              size_t strip_index) const;

    double solid_angle_correction(const size_t module_index,
                                  const size_t strip_index) const;

    double
    Antonio_solid_angle_correction_factor(const size_t module_index,
                                          const size_t strip_index) const;

    void create_flatfield_from_filelist(
        const std::vector<std::filesystem::path> &filelist,
        std::shared_ptr<MythenFileReader> file_reader);

    // TODO: maybe dont support reading at all
    // caveat user has to pass by value or move otherwise copied twice from
    // filesystem and then into member variable !!!!
    /**
     * @brief read normalized flatfield from file
     * @param filename name of file
     * @param file_reader custom file_reader to read flatfield file default:
     * ascii file with //TODO finish documentation
     */
    void read_normalized_flatfield_from_file(
        const std::string &filename,
        const std::shared_ptr<SimpleFileInterface> file_reader =
            std::make_shared<CustomFlatFieldFile>());

    void read_module_parameters_from_file(
        const std::filesystem::path &filename,
        const std::shared_ptr<SimpleFileInterface> file_reader =
            std::make_shared<InitialAngCalParametersFile>());

    void read_bad_channels_from_file(
        const std::filesystem::path &filename,
        const std::shared_ptr<SimpleFileInterface> file_reader =
            std::make_shared<CustomBadChannelsFile>());

    void set_bad_channels(const NDArray<bool, 1> &bad_channels_);

    NDArray<bool, 1> get_bad_channels() const;

    /**
     * set flatfield from NDArray
     */
    void set_normalized_flatfield(const NDArray<double, 2> &flat_field_);

    void set_normalized_flatfield(NDArray<double, 2> &&flat_field_);

    // TODO: should it be view? or copy?
    NDArray<double, 2> get_normalized_flatfield() const;

    NDView<double, 2> get_normalized_flatfield_view() const;

  private:
    /// @brief scale_factor is used to scale incident intensity to ~ one
    double scale_factor{1.0};

    std::shared_ptr<MythenDetectorSpecifications> mythen_detector;

    DGParameters DGparameters{};

    NDArray<double, 2> flat_field;

    NDArray<bool, 1> bad_channels{};

    /// @brief  soft window bounds [degrees] to exclude strips with lower
    /// exposure.
    std::pair<double, double> soft_window{3.0, 34.0};
};

} // namespace angcal