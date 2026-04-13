
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
     * @brief Constructor for FlatField class
     * @param mythen_detector_ ptr to MythenDetectorSpecifications storing all
     mythen specific parameters
     */
    FlatField(std::shared_ptr<MythenDetectorSpecifications> mythen_detector_);

    /**
     * @brief set scale factor to scale incident intensity to reasonable values
     * e.g. can be set to incident intensity of first acquisition
     */
    void set_scale_factor(const double scale_factor_);

    /**
     * @brief get scale factor
     */
    double get_scale_factor() const;

    /**
     * @brief set soft window bounds [degrees] to exclude strips with lower
     * exposure.
     */
    void set_soft_window(const std::pair<double, double> soft_window_);

    /**
     * @brief get soft window bounds [degrees]
     */
    std::pair<double, double> get_soft_window() const;

    /**
     * @brief correction factor for solid angle of strips
     * @param module_index index of module
     * @param strip_index local strip index of module e.g. 0-1279
     * @return correction factor for solid angle of strips
     */
    double solid_angle_correction_factor(const size_t module_index,
                                         const size_t strip_index) const;

    /*
    double
    Antonio_solid_angle_correction_factor(const size_t module_index,
                                          const size_t strip_index) const;

    double transverse_width_correction_factor(const size_t module_index,
                                              size_t strip_index) const;
    */

    /**
     * @brief creates normalized flatfield from list of files containing
     * flatfield acquisitions
     * @param filelist list of files containing flatfield acquisitions
     * @param file_reader custom file reader to read mythen acquisition files
     */
    void create_normalized_flatfield_from_filelist(
        const std::vector<std::filesystem::path> &filelist,
        std::shared_ptr<MythenFileReader> file_reader);

    // TODO: maybe dont support reading at all
    // caveat user has to pass by value or move otherwise copied twice from
    // filesystem and then into member variable !!!!
    /**
     * @brief read normalized flatfield and std from file
     * @param filename name of file
     * @param file_reader custom file_reader to read flatfield file default:
     * ascii file with
     */
    void read_normalized_flatfield_from_file(
        const std::string &filename,
        const std::shared_ptr<SimpleFileInterface> file_reader =
            std::make_shared<CustomFlatFieldFile>());

    /**
     * @brief reads DG calibration parameters from file
     * @param filename name of file
     * @param file_reader custom file_reader to read DG parameters file
     * (default InitialAngCalParametersFile:
     * following format module [module_index] center [center] +- [error]
     * conversion [conversion] +- [error] offset [offset] +- [error])
     */
    void read_module_parameters_from_file(
        const std::filesystem::path &filename,
        const std::shared_ptr<SimpleFileInterface> file_reader =
            std::make_shared<InitialAngCalParametersFile>());

    /**
     * @brief reads bad channels from file
     * @param filename bad channels filename
     * @param file_reader file_reader to read bad channels file (default:
     * CustomBadChannelsFile: following format line by line either single strip
     * index or block of strip indices seperated by '-' e.g. 0-10)
     */
    void read_bad_channels_from_file(
        const std::filesystem::path &filename,
        const std::shared_ptr<SimpleFileInterface> file_reader =
            std::make_shared<CustomBadChannelsFile>());

    void set_bad_channels(const NDArray<bool, 1> &bad_channels_);

    NDArray<bool, 1> get_bad_channels() const;

    /**
     * @brief set flatfield from NDArray
     */
    void set_normalized_flatfield(const NDArray<double, 2> &flat_field_);

    void set_normalized_flatfield(NDArray<double, 2> &&flat_field_);

    // TODO: should it be view? or copy?
    NDArray<double, 2> get_normalized_flatfield() const;

    NDView<double, 2> get_normalized_flatfield_view() const;

  private:
    /**
     * @brief elastic correction to diffraction angle
     * @param detector_angle detector angle [degrees]
     * @return elastic correction to diffraction angle [degrees]
     */
    double elastic_correction(const double detector_angle) const;

    /**
     * @brief calculates diffraction angle from DG module parameters (used in
     * Beer's Law)
     * @param detector_angle detector position [degrees]
     * @param strip_index local strip index of module e.g. 0-1279
     * @param distance_to_strip distance to strip [given in strips]
     */
    double diffraction_angle_from_DG_parameters(
        const size_t module_index, const double detector_angle,
        size_t strip_index, const double distance_to_strip) const;

  private:
    /// @brief scale_factor is used to scale incident intensity to ~ one
    double scale_factor{1.0};

    /// @brief detector specifications of Mythen detector
    std::shared_ptr<MythenDetectorSpecifications> mythen_detector;

    /// @brief DG module parameters stored as 2D array of shape (num_modules, 3)
    /// where the three parameters are center, conversion and offset
    DGParameters DGparameters{};

    /// @brief normalized flatfield stored as 2D array of shape (num_strips, 2)
    /// where the two values are normalized flatfield value and standard
    /// deviation
    NDArray<double, 2> flat_field;

    /// @brief boolean array of size of total number of strips/channels in
    /// detector, stores 'TRUE' if strip is a bad channel otherwise 'FALSE'
    NDArray<bool, 1> bad_channels{};

    /// @brief  soft window bounds [degrees] to exclude strips with lower
    /// exposure.
    std::pair<double, double> soft_window{3.0, 34.0};
};

} // namespace angcal