#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "FlatField.hpp"
#include "MythenDetectorSpecifications.hpp"
#include "MythenFileReader.hpp"
#include "aare/NDArray.hpp"
#include "helpers/FileInterface.hpp"

using namespace aare;

namespace angcal {

using parameters =
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>;

// TODO: check if interpretation and units are correct
// historical DG parameters
/**
 * historical Detector Group parameters
 */
struct DGParameters {

    DGParameters() = default;

    DGParameters(const ssize_t n_modules) {
        parameters = NDArray<double, 2>(std::array<ssize_t, 2>{n_modules, 3});
    }

    double &operator()(const size_t module_index,
                       const size_t parameter_index) {
        return parameters(module_index, parameter_index);
    }

    double operator()(const size_t module_index,
                      const size_t parameter_index) const {
        return parameters(module_index, parameter_index);
    }

    /**
     * orthogonal projection of sample onto
     * detector (given in strip number) [mm]
     * D/pitch
     */
    double &centers(const size_t module_index) {
        return parameters(module_index, 0);
    }

    double centers(const size_t module_index) const {
        return parameters(module_index, 0);
    }

    /**
     * pitch/(normal distance from sample
     * to detector (R)) [mm]
     * used for easy conversion
     */
    double &conversions(const size_t module_index) {
        return parameters(module_index, 1);
    }

    double conversions(const size_t module_index) const {
        return parameters(module_index, 1);
    }

    /** position of strip zero relative to sample [degrees] phi
     * 180/pi*D/R TODO: expected an arcsin(D/R)?
     */
    double &offsets(const size_t module_index) {
        return parameters(module_index, 2);
    }

    double offsets(const size_t module_index) const {
        return parameters(module_index, 2);
    }

    NDArray<double, 2> parameters{};
};

/**
 * geometric parameters
 */
struct EEParameters {

    EEParameters(const ssize_t n_modules) {
        parameters = NDArray<double, 2>(std::array<ssize_t, 2>{n_modules, 3});
    }

    double &operator()(const size_t module_index,
                       const size_t parameter_index) {
        return parameters(module_index, parameter_index);
    }

    double operator()(const size_t module_index,
                      const size_t parameter_index) const {
        return parameters(module_index, parameter_index);
    }

    /**
     * normal distance between sample and detector (R)
     */
    double &normal_distances(const size_t module_index) {
        return parameters(module_index, 0);
    }

    double normal_distances(const size_t module_index) const {
        return parameters(module_index, 0);
    }

    /**
     * distances between intersection point of sample normal and module origin
     * (D)
     */
    double &module_center_distances(const size_t module_index) {
        return parameters(module_index, 1);
    }
    double module_center_distances(const size_t module_index) const {
        return parameters(module_index, 1);
    }

    /** angles between undiffracted beam and orthogonal sample projection on
     * detector (phi)
     */
    double &angles(const size_t module_index) {
        return parameters(module_index, 2);
    }

    double angles(const size_t module_index) const {
        return parameters(module_index, 2);
    }

    NDArray<double, 2> parameters{};
};

class AngleCalibration {

  public:
    AngleCalibration(
        std::shared_ptr<MythenDetectorSpecifications> mythen_detector_,
        std::shared_ptr<FlatField> flat_field_,
        std::shared_ptr<MythenFileReader> mythen_file_reader_,
        std::optional<std::shared_ptr<SimpleFileInterface>> custom_file_ptr_ =
            std::nullopt);

    /** set the histogram bin width [degrees] */
    void set_histogram_bin_width(double bin_width);

    double get_histogram_bin_width() const;

    ssize_t get_new_num_bins() const;

    /** reads the historical Detector Group (DG) parameters from file **/
    void read_initial_calibration_from_file(const std::string &filename);

    const DGParameters &get_DGparameters() const; // TODO: do i want a
                                                  // reference?

    NDArray<double, 1> get_new_photon_counts() const;

    NDArray<double, 1> get_new_statistical_errors() const;

    /** converts DG parameters to easy EE parameters e.g.geometric
     * parameters */
    EEParameters convert_to_EE_parameters() const;

    std::tuple<double, double, double>
    convert_to_EE_parameters(const size_t module_index) const;

    std::tuple<double, double, double>
    convert_to_EE_parameters(const double center, const double conversion,
                             const double offset) const;

    /*
    //converts DG parameters to BC parameters e.g. best computing
     parameters
    parameters convert_to_BC_parameters() const;
    */

    /**
     * calculates new histogram with fixed sized angle bins
     * for several acquisitions at different detector angles for given frame
     * indices
     * @param start_frame_index, end_frame_index gives range of frames
     */
    void
    calculate_fixed_bin_angle_width_histogram(const size_t start_frame_index,
                                              const size_t end_frame_index);

    void write_to_file(const std::string &filename,
                       const bool store_nonzero_bins = false,
                       const std::filesystem::path &filepath =
                           std::filesystem::current_path()) const;

  private:
    /** calculates diffraction angle from EE module parameters (used in
     * Beer's Law)
     * @param strip_index local strip index of module
     */
    double diffraction_angle_from_EE_parameters(
        const double module_center_distance, const double normal_distance,
        const double angle, const size_t strip_index,
        const double distance_to_strip = 0) const;

    /** calculates diffraction angle from EE module parameters (used in
     * Beer's Law)
     * @param center module center
     * @param conversion module conversion
     * @param offset module offset
     * @param strip_index local strip index of module
     * @param distance_to_strip distance to strip given by strip_index and
     * module -> note needs to be small enough to be in the respective module
     */
    double diffraction_angle_from_DG_parameters(
        const double center, const double conversion, const double offset,
        const size_t strip_index, const double distance_to_strip = 0) const;

    /** calculated the strip width expressed as angle [degrees]
     * @param strip_index local strip index of module
     */
    double angular_strip_width_from_DG_parameters(
        const double center, const double conversion, const double offset,
        const size_t local_strip_index) const;

    double angular_strip_width_from_EE_parameters(
        const double module_center_distance, const double normal_distance,
        const double angle, const size_t local_strip_index) const;

    /** converts global strip index to local strip index of that module */
    size_t global_to_local_strip_index_conversion(
        const size_t global_strip_index) const;

    /**
     * redistributes photon counts with of histogram using one bin per strip
     * to histogram with fixed size angle bins
     * @param frame MythenFrame storing data from image
     * @param bin_counts accumulate new photon counts
     * @param new_statistical_weights accumulate new statistical weights
     * @param new_errors accumulate new_errors
     */
    void redistribute_photon_counts_to_fixed_angle_bins(
        const MythenFrame &frame, NDView<double, 1> bin_counts,
        NDView<double, 1> new_statistical_weights, NDView<double, 1> new_errors,
        NDView<double, 1> inverse_nromalized_flatfield) const;

  private:
    DGParameters DGparameters;

    std::shared_ptr<MythenDetectorSpecifications> mythen_detector{};

    std::shared_ptr<FlatField> flat_field{};

    NDArray<double, 1> new_photon_counts{};
    NDArray<double, 1> new_photon_count_errors{};

    double histogram_bin_width = 0.0036; // [degrees]

    ssize_t num_bins{};

    std::shared_ptr<MythenFileReader>
        mythen_file_reader{}; // TODO replace by FileInterface ptr

    std::optional<std::shared_ptr<SimpleFileInterface>> custom_file_ptr{};
};

} // namespace angcal
