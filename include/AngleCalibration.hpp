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
    /**
     * ptr to MythenDetectorSpecifications storing all specific variables of
     mythen detector
     * ptr to FlatField class
     * @param file_path_ base_path to acquisition files
     * @param mythen_file_reader optional, pass if you use custom acquisition
     files - default: reads hdf5 files
     * @param custom_file_ptr optional, pass if you use custom files to store
     initial angle parameters - default: initial angle parameters supports
     * format module [module_index] center [center] +- [error] conversion
     [conversion] +- [error] offset [offset] +- [error]
     */
    AngleCalibration(
        std::shared_ptr<MythenDetectorSpecifications> mythen_detector_,
        std::shared_ptr<FlatField> flat_field_,
        const std::filesystem::path &file_path_,
        std::optional<std::shared_ptr<MythenFileReader>> mythen_file_reader_ =
            std::nullopt,
        std::optional<std::shared_ptr<SimpleFileInterface>> custom_file_ptr_ =
            std::nullopt);

    /** set the histogram bin width [degrees] */
    void set_histogram_bin_width(double bin_width);

    double get_histogram_bin_width() const;

    ssize_t get_base_peak_ROI_num_bins() const;

    ssize_t get_base_peak_ROI() const;

    /** reads the historical Detector Group (DG) parameters from file
     * @warning only works if member m_custom_detectot_file_ptr supports reading
     * the format
     */
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
     * @param file_list vector of file_names
     */
    void calculate_fixed_bin_angle_width_histogram(
        const std::vector<std::string> &file_list);

    /**
     * calculates new histogram with fixed sized angle bins
     * for one acquisition
     */
    NDArray<double, 1>
    calculate_fixed_bin_angle_width_histogram(const std::string &file_name);

    /**
     * calibrates the DG parameters
     * @param file_list vector of file_names of acquisition files - detector
     * position should change
     * @param base_peak_angle angle of the selcted base peak given in degrees
     * (generally base peaks are tabulated peaks based on sample) //TODO: in
     * radians or degrees?
     */
    void calibrate(const std::vector<std::string> &file_list,
                   const double base_peak_angle);

    void calibrate(const std::vector<std::string> &file_list_,
                   const double base_peak_angle_, const size_t module_index);

    void write_to_file(const std::string &filename,
                       const bool store_nonzero_bins = false,
                       const std::filesystem::path &filepath =
                           std::filesystem::current_path()) const;

    double calculate_similarity_of_peaks(const size_t module_index) const;

    /** check if base_peak is cituated in module */
    bool base_peak_is_in_module(const size_t module_index,
                                const double detector_angle) const;

    void set_base_peak_angle(const double base_peak_angle_);

    bool module_is_disconnected(const size_t module_index) const;

    /**
     * redistributes photon counts around region of interest of base peak to
     * fixed angle width bins
     * @param module_index index of module
     * @param frame MythenFrame storing data from acquisition
     * @param new_fixed_angle_width_bin_histogram accumulate new photon counts
     * (histogram covers region of interest around base peak)
     */
    template <bool base_peak_ROI_only = false>
    void redistribute_photon_counts_to_fixed_angle_width_bins(
        const size_t module_index, const MythenFrame &frame,
        NDView<double, 1> fixed_angle_width_bins_photon_counts,
        NDView<double, 1> fixed_angle_width_bins_photon_counts_variance,
        std::optional<NDView<double, 1>> S0 = std::nullopt,
        std::optional<NDView<double, 1>> S1 = std::nullopt,
        std::optional<NDView<double, 1>> S2 = std::nullopt) const;

    NDArray<double, 1> redistribute_photon_counts_to_fixed_angle_width_bins(
        const MythenFrame &frame) const;

    /** calculates diffraction angle from EE module parameters (used in
     * Beer's Law)
     * @param center module center
     * @param conversion module conversion
     * @param offset module offset
     * @param strip_index local strip index of module
     * @param distance_to_strip distance to strip given by strip_index and
     * module -> note needs to be small enough to be in the respective module
     * @return diffraction angle in degrees
     */
    double diffraction_angle_from_DG_parameters(
        const size_t module_index, const double detector_angle,
        const size_t strip_index, const double distance_to_strip = 0) const;

  private:
    /** calculates diffraction angle from EE module parameters (used in
     * Beer's Law)
     * @param strip_index local strip index of module
     */
    double diffraction_angle_from_EE_parameters(
        const double module_center_distance, const double normal_distance,
        const double angle, const double detector_angle,
        const size_t strip_index, const double distance_to_strip = 0) const;

    /** calculated the strip width expressed as angle [degrees]
     * @param strip_index local strip index of module
     */
    double angular_strip_width_from_DG_parameters(
        const size_t module_index, const size_t local_strip_index) const;

    double angular_strip_width_from_EE_parameters(
        const double module_center_distance, const double normal_distance,
        const double angle, const size_t local_strip_index) const;

    /** converts global strip index to local strip index of that module */
    size_t global_to_local_strip_index_conversion(
        const size_t global_strip_index) const;

    /**
     * compares multiple ROI around peaks and calculates similarity/variance
     * based on goodness_of_fit and weighted average
     * @param S0 photon_varaince over all runs
     * @param S1 photon_variance*photon_count over all runs
     * @param S2 photon_variance*photon_count² over all runs
     * @return similarity criterion
     */
    double similarity_criterion(const NDView<double, 1> S0,
                                const NDView<double, 1> S1,
                                const NDView<double, 1> S2) const;

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

    // TODO: also a bad design - make shift_parameter configurable - consider
    // having second class Optimization
    void optimization_algorithm(const size_t module_index,
                                const double shift_parameter1 = 0.01,
                                const double shift_parameter2 = 0.005);

  private:
    DGParameters DGparameters;

    std::shared_ptr<MythenDetectorSpecifications> mythen_detector{};

    std::shared_ptr<FlatField> flat_field{};

    NDArray<double, 1> new_photon_counts{};
    NDArray<double, 1> new_photon_count_errors{};

    double histogram_bin_width = 0.0036; // [degrees]

    // TODO add in constructor or setter
    ssize_t base_peak_roi = 50; // base_peak_as_bin_index +/- base_peak_roi -
                                // region of interest around base peak

    ssize_t num_bins{};

    double base_peak_angle{}; // center of base peak in angles to use for
                              // calibration

    std::vector<std::string>
        file_list{}; // list of acquisition files used for calibration

    std::shared_ptr<MythenFileReader>
        mythen_file_reader{}; // TODO replace by FileInterface ptr

    std::shared_ptr<SimpleFileInterface> custom_file_ptr =
        std::make_shared<InitialAngCalParametersFile>();
};

template <bool base_peak_ROI_only>
void AngleCalibration::redistribute_photon_counts_to_fixed_angle_width_bins(
    const size_t module_index, const MythenFrame &frame,
    NDView<double, 1> fixed_angle_width_bins_photon_counts,
    NDView<double, 1> fixed_angle_width_bins_photon_counts_variance,
    std::optional<NDView<double, 1>> S0, std::optional<NDView<double, 1>> S1,
    std::optional<NDView<double, 1>> S2) const {
    double left_boundary_roi_base_peak =
        mythen_detector->min_angle(); // dummy values
    double right_boundary_roi_base_peak =
        mythen_detector->max_angle(); // dummy values

    if constexpr (base_peak_ROI_only) {
        ssize_t base_peak_as_bin_index = static_cast<ssize_t>(
            (base_peak_angle) /
            histogram_bin_width); // TODO: in antonios code actually rounded
                                  // to nearest integer
        left_boundary_roi_base_peak =
            (base_peak_as_bin_index - base_peak_roi - 0.5) *
            histogram_bin_width; // in degrees
        right_boundary_roi_base_peak =
            (base_peak_as_bin_index + base_peak_roi + 0.5) *
            histogram_bin_width; // in degrees
    }
    for (size_t strip_index = 0;
         strip_index < mythen_detector->strips_per_module(); ++strip_index) {

        size_t global_strip_index =
            module_index * mythen_detector->strips_per_module() +
            strip_index; // TODO: is this really correct - check sign

        if (mythen_detector->get_bad_channels()(global_strip_index)) {
            continue; // skip bad channels
        }

        double left_strip_boundary_angle = diffraction_angle_from_DG_parameters(
            module_index, frame.detector_angle, strip_index, -0.5);

        double right_strip_boundary_angle =
            diffraction_angle_from_DG_parameters(
                module_index, frame.detector_angle, strip_index, 0.5);

        if (base_peak_ROI_only &&
            (left_strip_boundary_angle > right_boundary_roi_base_peak ||
             right_strip_boundary_angle < left_boundary_roi_base_peak)) {
            continue; // skip strip if not in ROI
        }

        double flatfield_normalized_photon_counts =
            (frame.photon_counts(global_strip_index) + 1) *
            flat_field->get_inverse_normalized_flatfield()(global_strip_index);

        double some_flatfield_error = 1.0; // TODO: some dummy value - implement

        // I guess it measures the
        // expcected noise - where is the formula -used as the variance
        double photon_counts_variance =
            1. / (std::pow(flatfield_normalized_photon_counts, 2) *
                  (1. / (frame.photon_counts(global_strip_index) + 1) +
                   std::pow(some_flatfield_error *
                                flat_field->get_inverse_normalized_flatfield()(
                                    global_strip_index),
                            2)));

        double strip_width_angle =
            angular_strip_width_from_DG_parameters(module_index, strip_index);

        double correction_factor = histogram_bin_width / strip_width_angle;

        ssize_t left_bin_index_covered_by_strip = static_cast<ssize_t>(
            left_strip_boundary_angle /
            histogram_bin_width); // -
                                  // mythen_detector->min_angle()/histogram_bin_width

        ssize_t right_bin_index_covered_by_strip = static_cast<ssize_t>(
            right_strip_boundary_angle /
            histogram_bin_width); // -
                                  // mythen_detector->min_angle()/histogram_bin_width

        if constexpr (base_peak_ROI_only) {

            left_bin_index_covered_by_strip = std::max(
                static_cast<ssize_t>(base_peak_angle / histogram_bin_width) -
                    base_peak_roi,
                left_bin_index_covered_by_strip); // -
                                                  // mythen_detector->min_angle()/histogram_bin_width
            right_bin_index_covered_by_strip = std::min(
                static_cast<ssize_t>(base_peak_angle / histogram_bin_width) +
                    base_peak_roi,
                right_bin_index_covered_by_strip);
        }

        size_t proper_bin_index =
            0; // the computed bin indices dont start at zero but arange around
               // zero depending on the sign if the diffraction angle
        for (ssize_t bin_index = left_bin_index_covered_by_strip;
             bin_index <= right_bin_index_covered_by_strip; ++bin_index) {

            double bin_coverage =
                std::abs(std::min(right_strip_boundary_angle,
                                  (bin_index + 0.5) * histogram_bin_width) -
                         std::max(left_strip_boundary_angle,
                                  (bin_index - 0.5) * histogram_bin_width));

            double bin_coverage_factor = bin_coverage / histogram_bin_width;

            // convert to bin index
            if constexpr (base_peak_ROI_only) {
                proper_bin_index =
                    bin_index - (static_cast<ssize_t>(base_peak_angle /
                                                      histogram_bin_width) -
                                 base_peak_roi); // bin index starts at zero
            } else {
                proper_bin_index =
                    bin_index -
                    static_cast<ssize_t>(
                        mythen_detector->min_angle() /
                        histogram_bin_width); // bin index starts at zero
            }

            fixed_angle_width_bins_photon_counts(proper_bin_index) +=
                flatfield_normalized_photon_counts * correction_factor *
                bin_coverage_factor;

            fixed_angle_width_bins_photon_counts_variance(proper_bin_index) +=
                photon_counts_variance /
                std::pow(correction_factor * bin_coverage_factor,
                         2); // TODO no idea if this is correct

            // S_index = sum_i^num_runs
            // photon_count^index*photon_variance
            if (S0.has_value()) {
                S0.value()(proper_bin_index) +=
                    fixed_angle_width_bins_photon_counts_variance(
                        proper_bin_index);
            }
            if (S1.has_value()) {
                S1.value()(proper_bin_index) +=
                    fixed_angle_width_bins_photon_counts(proper_bin_index) *
                    fixed_angle_width_bins_photon_counts_variance(
                        proper_bin_index);
            }
            if (S2.has_value()) {
                S2.value()(proper_bin_index) +=
                    fixed_angle_width_bins_photon_counts(proper_bin_index) *
                    fixed_angle_width_bins_photon_counts(proper_bin_index) *
                    fixed_angle_width_bins_photon_counts_variance(
                        proper_bin_index);
            }
        }
    }
}

} // namespace angcal
