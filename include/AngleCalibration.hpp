#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "FlatField.hpp"
#include "MythenDetectorSpecifications.hpp"
#include "MythenFileReader.hpp"
#include "Parameters.hpp"
#include "aare/NDArray.hpp"
#include "helpers/FileInterface.hpp"

#ifdef ANGCAL_PLOT
#include "plot_histogram.hpp"
using PlotHandle =
    std::shared_ptr<angcal::PlotCalibrationProcess>; // Or reference wrapper
#else
using PlotHandle = std::nullptr_t;
#endif

using namespace aare;

namespace angcal {

class AngleCalibration {

  public:
    /**
     * @param mythen_detector_ ptr to MythenDetectorSpecifications storing all
     mythen specific parameters
     * @param flat_field_ ptr to FlatField storing (Note calculate inverse
     normalized flatfield before usage)
     * @param mythen_file_reader optional, pass if you use custom acquisition
     files - default: reads hdf5 files
     * @param custom_file_ptr optional, pass if you use custom files to store
     initial angle parameters - default: initial angle parameters supports
     * following format module [module_index] center [center] +- [error]
     conversion [conversion] +- [error] offset [offset] +- [error]
     */
    AngleCalibration(
        std::shared_ptr<MythenDetectorSpecifications> mythen_detector_,
        std::shared_ptr<FlatField> flat_field_,
        std::optional<std::shared_ptr<MythenFileReader>> mythen_file_reader_ =
            std::nullopt,
        std::optional<std::shared_ptr<SimpleFileInterface>> custom_file_ptr_ =
            std::nullopt);

    /** set the histogram bin width [degrees]
     * default 0.0036
     */
    void set_histogram_bin_width(double bin_width);

    double get_histogram_bin_width() const;

    /** number of bins of fixed angle width bin histogram */
    ssize_t new_number_of_bins() const;

    /** number of bins covered of base peak region of interest
     * default 101 bins
     */
    ssize_t get_base_peak_ROI_num_bins() const;

    /** set base peak region of interest e.g. [base_peak - base_peak_ROI,
     * base_peak + base_peak_ROI] given in angles [degrees]
     * default (0.05)
     * */
    void set_base_peak_ROI(const double base_peak_roi_);

    /** base peak region of interest e.g. [base_peak - base_peak_ROI,
     * base_peak + base_peak_ROI] given in angles [degrees]
     * default (0.05)
     * */
    double get_base_peak_ROI() const;

    std::shared_ptr<MythenDetectorSpecifications>
    get_detector_specifications() const;

    /** reads the historical Detector Group (DG) parameters from file
     * @warning only works if member m_custom_detectot_file_ptr supports reading
     * the format
     */
    void read_initial_calibration_from_file(const std::string &filename);

    const DGParameters &get_DGparameters() const; // TODO do I need that?

    /**
     * @brief calibrates the BC (best computing) parameters
     * @param file_list vector of file_names of acquisition files - detector
     * position should change
     * @param base_peak_angle angle of the selcted base peak given in degrees
     * (generally base peaks are tabulated peaks based on sample)
     */
    void calibrate(const std::vector<std::string> &file_list,
                   const double base_peak_angle);

    /**
     * @brief calibrates the BC (best computing) parameters for one module
     * @param file_list vector of file_names of acquisition files - detector
     * position should change
     * @param base_peak_angle angle of the selcted base peak given in degrees
     * (generally base peaks are tabulated peaks based on sample)
     * @param module_index index of module to be calibrated
     */
    void calibrate(const std::vector<std::string> &file_list_,
                   const double base_peak_angle_, const size_t module_index);

    // deprecated
    template <typename T>
    void write_to_file(const std::string &filename,
                       const std::filesystem::path &filepath,
                       const NDView<T, 1> data) const;

    /** @brief check if base_peak ROI overlaps with module region
     * @param detector_angle: detector position at which acquisition is taken
     * @param bounds_in_angle: boundary that is acceptable [degrees]- per
     * default the base peak ROI width is used (adjustable for debugging) //TODO
     * deprecated
     * @return true if [base_peak - bound_in_angle, base_peak + bound_in_angle]
     * overlap with the module region in degrees
     */
    bool base_peak_is_in_module(
        const size_t module_index, const double detector_angle,
        std::optional<double> bounds_in_angles = std::nullopt) const;

    /**
     * @brief set angle of base_peak_module to calibrate
     */
    void set_base_peak_angle(const double base_peak_angle_);

    double get_base_peak_angle() const;

    /**
     * @brief check if a module only has bad channels
     */
    bool module_is_disconnected(const size_t module_index) const;

    /**
     * @brief redistribute photon counts to fixed angle width bins
     */
    NDArray<double, 1> redistribute_photon_counts_to_fixed_angle_width_bins(
        const MythenFrame &frame) const;

    /**
     * @brief redistribute photon counts to fixed angle width bins for the
     * specified module region
     * @param module_index: index of module region to redistribute
     */
    NDArray<double, 1> redistribute_photon_counts_to_fixed_angle_width_bins(
        const MythenFrame &frame, const size_t module_index) const;

    /**
     * @brief redistribute photon counts to fixed angle width bins for base peak
     * region
     * @param module_index: index of module region to redistribute
     */
    NDArray<double, 1> redistributed_photon_counts_in_base_peak_ROI(
        const MythenFrame &frame, const size_t module_index) const;

    /** @brief calculates diffraction angle from DG module parameters (used in
     * Beer's Law)
     * @param module_index
     * @param detector_angle detector position at which acquisition is taken
     * @param strip_index local strip index of module
     * @param distance_to_strip distance to strip [given in strip bins]
     * @return diffraction angle in degrees
     */
    double diffraction_angle_from_DG_parameters(
        const size_t module_index, const double detector_angle,
        size_t strip_index, const double distance_to_strip = 0) const;

    /** @brief calculates diffraction angle from BC (best computing) module
     * parameters (used in Beer's Law)
     * @param module_index
     * @param detector_angle detector position at which acquisition is taken
     * @param strip_index local strip index of module
     * @param distance_to_strip distance to strip [given in strip bins]
     * @return diffraction angle in degrees
     */
    double diffraction_angle_from_BC_parameters(
        const size_t module_index, const double detector_angle,
        const size_t strip_index, const double distance_to_strip = 0) const;

    /** @brief calculates diffraction angle from EE module parameters (used in
     * Beer's Law)
     * @param module_index
     * @param detector_angle detector position at which acquisition is taken
     * @param strip_index local strip index of module
     * @param distance_to_strip distance to strip [given in strip bins]
     * @return diffraction angle in degrees
     */
    double diffraction_angle_from_EE_parameters(
        const double module_center_distance, const double normal_distance,
        const double angle, const double detector_angle,
        const size_t strip_index, const double distance_to_strip = 0) const;

  private:
    /** @brief calculated the strip width expressed as angle [degrees] from DG
     * module parameters
     * @param moudle_index index of module
     * @param strip_index local strip index of module
     */
    double
    angular_strip_width_from_DG_parameters(const size_t module_index,
                                           size_t local_strip_index) const;

    /** @brief calculated the strip width expressed as angle [degrees] from BC
     * module parameters
     * @param moudle_index index of module
     * @param strip_index local strip index of module
     */
    double angular_strip_width_from_BC_parameters(
        const size_t module_index, const size_t local_strip_index) const;

    /** @brief calculated the strip width expressed as angle [degrees] from EE
     * module parameters
     * @param moudle_index index of module
     * @param strip_index local strip index of module
     */
    double angular_strip_width_from_EE_parameters(
        const double module_center_distance, const double normal_distance,
        const double angle, const size_t local_strip_index) const;

    // TODO deprecated
    /** @brief converts global strip index to local strip index of that module
     */
    size_t global_to_local_strip_index_conversion(
        const size_t global_strip_index) const;

    /**
     * @brief redistributes photon counts around designated region to
     * fixed angle width bins
     * @tparam base_peak_ROI_only: false (use module region), true: (only use
     * base peak ROI)
     * @param module_index index of module
     * @param frame data from acquisition
     * @param fixed_angle_width_bin_photon_counts stores redistributed photon
     * counts for fixed angle width bin
     * @param fixed_angle_width_bins_photon_counts_variance stores variance
     * @param S0, S1, S2 used to calculate similarity criterion between peaks
     */
    template <bool base_peak_ROI_only = false>
    void redistribute_photon_counts_to_fixed_angle_width_bins(
        const size_t module_index, const MythenFrame &frame,
        NDView<double, 1> fixed_angle_width_bins_photon_counts,
        NDView<double, 1> fixed_angle_width_bins_photon_counts_variance,
        std::optional<NDView<double, 1>> S0 = std::nullopt,
        std::optional<NDView<double, 1>> S1 = std::nullopt,
        std::optional<NDView<double, 1>> S2 = std::nullopt) const;

    /**
     * @brief redistributes photon counts to fixed angle width bins around base
     * peak for all acquisitions in file_list and calculates the similarity
     * between the found base peaks
     */
    double calculate_similarity_of_peaks(const size_t module_index,
                                         PlotHandle gp = nullptr) const;

    /**
     * compares multiple base peak ROIS from different acquisitions and
     * calculate similarity/variance based on goodness_of_fit and weighted
     * average
     * @param S0 photon_varaince over all runs
     * @param S1 photon_variance*photon_count over all runs
     * @param S2 photon_variance*photon_count² over all runs
     * @return similarity criterion
     */
    double similarity_criterion(const NDView<double, 1> S0,
                                const NDView<double, 1> S1,
                                const NDView<double, 1> S2,
                                const size_t num_runs) const;

    // TODO: also a bad design - make shift_parameter configurable - consider
    // having second class Optimization
    /**
     * optimizes BC parameters of given module based on similarity criterion
     * @param shift_parameter1 parameter step for parameter angle between module
     * center and module normal (angle_center_module_normal)
     * @param shift_parameter2 parameter step for parameter distance between
     * module center and sample (module_center_sample_distances)
     */
    void optimization_algorithm(const size_t module_index,
                                PlotHandle gp = nullptr,
                                const double shift_parameter1 = 0.01,
                                const double shift_parameter2 = 0.005);

  private:
    DGParameters DGparameters{};
    BCParameters BCparameters{};

    /**
     * detector specifications
     */
    std::shared_ptr<MythenDetectorSpecifications> mythen_detector{};

    std::shared_ptr<FlatField> flat_field{};

    /**
     * bin width of fixed angle bin histogram [degrees]
     */
    double histogram_bin_width = 0.0036;

    /**
     * region of interest of base peak given in angles
     * [base_peak_angle - base_peak_roi, base_peak_angle +
     * base_peak_roi]
     */
    double base_peak_roi = 0.05;

    // TODO maybe deprecated - only compute in member function
    ssize_t num_bins{};

    /**
     * center of base peak to use for calibration [degrees]
     */
    double base_peak_angle{};

    /**
     * list of acquisition files used for calibration
     */
    std::vector<std::string> file_list{};

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

    NDArray<double, 1> fraction_covered_by_strip(
        std::array<ssize_t, 1>{new_number_of_bins()},
        0.0); // fraction of strip

    if constexpr (base_peak_ROI_only) {

        left_boundary_roi_base_peak = (base_peak_angle - base_peak_roi +
                                       0.5 * histogram_bin_width); // in degrees
        right_boundary_roi_base_peak =
            (base_peak_angle + base_peak_roi -
             0.5 * histogram_bin_width); // in degrees
    }

    for (size_t strip_index = 0;
         strip_index < mythen_detector->strips_per_module(); ++strip_index) {

        size_t global_strip_index =
            module_index * mythen_detector->strips_per_module() +
            strip_index; // TODO: is this really correct - check sign

        if (mythen_detector->get_bad_channels()(global_strip_index)) {
            continue; // skip bad channels
        }

        double left_strip_boundary_angle = diffraction_angle_from_BC_parameters(
            module_index, frame.detector_angle, strip_index, -0.5); // -0.5

        double right_strip_boundary_angle =
            diffraction_angle_from_BC_parameters(
                module_index, frame.detector_angle, strip_index, +0.5);

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
        // expcected noise - where is the formula - used as the variance
        double photon_counts_variance =
            1. / (std::pow(flatfield_normalized_photon_counts, 2) *
                  (1. / (frame.photon_counts(global_strip_index) + 1) +
                   std::pow(some_flatfield_error *
                                flat_field->get_inverse_normalized_flatfield()(
                                    global_strip_index),
                            2)));

        double strip_width_angle =
            std::abs(right_strip_boundary_angle - left_strip_boundary_angle);
        // angular_strip_width_from_BC_parameters(module_index, strip_index);

        /*
        double correction_factor =
            histogram_bin_width >= strip_width_angle
                ? 1.0
                : histogram_bin_width /
                      strip_width_angle; // coverage factor of one bin of the
                                         // strip
        */

        double photon_counts_per_bin = flatfield_normalized_photon_counts *
                                       histogram_bin_width / strip_width_angle;

        double photon_counts_variance_per_bin =
            photon_counts_variance * (strip_width_angle / histogram_bin_width) *
            (strip_width_angle / histogram_bin_width);

        ssize_t left_bin_index_covered_by_strip{};

        ssize_t right_bin_index_covered_by_strip{};

        if constexpr (base_peak_ROI_only) {
            left_bin_index_covered_by_strip =
                static_cast<ssize_t>(round(std::max(left_boundary_roi_base_peak,
                                                    left_strip_boundary_angle) /
                                           histogram_bin_width));
            right_bin_index_covered_by_strip = static_cast<ssize_t>(
                round(std::min(right_boundary_roi_base_peak,
                               right_strip_boundary_angle) /
                      histogram_bin_width));
        } else {
            left_bin_index_covered_by_strip = static_cast<ssize_t>(
                round(left_strip_boundary_angle /
                      histogram_bin_width)); // in Antonios code its rounded
            right_bin_index_covered_by_strip = static_cast<ssize_t>(
                round(right_strip_boundary_angle / histogram_bin_width));
        }

        size_t proper_bin_index =
            0; // the computed bin indices dont start at zero but arange around
               // zero depending on the sign if the diffraction angle

        for (ssize_t bin_index = left_bin_index_covered_by_strip;
             bin_index < right_bin_index_covered_by_strip; ++bin_index) {

            // some strips dont cover an entire bin or extend over multiple bins
            double bin_coverage =
                std::abs(std::min(right_strip_boundary_angle,
                                  (bin_index + 0.5) * histogram_bin_width) -
                         std::max(left_strip_boundary_angle,
                                  (bin_index - 0.5) * histogram_bin_width));

            double bin_coverage_factor =
                bin_coverage / histogram_bin_width; // how much of the strip is
                                                    // covered by the bin
            // convert to bin index
            if constexpr (base_peak_ROI_only) {
                proper_bin_index =
                    bin_index -
                    (static_cast<ssize_t>(
                        (base_peak_angle - base_peak_roi) /
                        histogram_bin_width)); // bin index starts at zero
            } else {
                proper_bin_index =
                    bin_index -
                    static_cast<ssize_t>(
                        mythen_detector->min_angle() /
                        histogram_bin_width); // bin index starts at zero
            }

            double corrected_photon_counts =
                photon_counts_per_bin; //* bin_coverage_factor;

            fraction_covered_by_strip(proper_bin_index) += bin_coverage_factor;

            double corrected_photon_counts_variance =
                photon_counts_variance_per_bin; //*
            // std::pow(bin_coverage_factor,
            // 2); // this seems to be a correction factor but does not really
            // redistribute it as it is not multiplied but only multiplied once

            fixed_angle_width_bins_photon_counts(proper_bin_index) +=
                bin_coverage_factor * corrected_photon_counts *
                corrected_photon_counts_variance;

            fixed_angle_width_bins_photon_counts_variance(proper_bin_index) +=
                corrected_photon_counts_variance * bin_coverage_factor;
        }
    }

    // S_index = sum_i^num_runs
    // photon_count^index*photon_variance
    for (ssize_t i = 0; i < fixed_angle_width_bins_photon_counts.size(); ++i) {
        fixed_angle_width_bins_photon_counts(i) /=
            fixed_angle_width_bins_photon_counts_variance(
                i); // y_k what exactly is
                    // that?

        if (S0.has_value()) {
            S0.value()(i) += fixed_angle_width_bins_photon_counts_variance(i);
        }
        if (S1.has_value()) {
            S1.value()(i) += fixed_angle_width_bins_photon_counts(i) *
                             fixed_angle_width_bins_photon_counts_variance(i);
        }
        if (S2.has_value()) {
            S2.value()(i) += fixed_angle_width_bins_photon_counts(i) *
                             fixed_angle_width_bins_photon_counts(i) *
                             fixed_angle_width_bins_photon_counts_variance(i);
            ;
        }
    }

    /*
    auto data_file_path =
        std::filesystem::current_path().parent_path() / "data";
    if (!std::filesystem::exists(data_file_path))
        std::filesystem::create_directories(data_file_path);
    write_to_file("fraction_covered_by_strip_" + std::to_string(module_index) +
                      ".txt",
                  data_file_path,
                  fraction_covered_by_strip.view()); // TODO: remove this - only
                                                     // for debugging
    */
}

template <typename T>
void AngleCalibration::write_to_file(const std::string &filename,
                                     const std::filesystem::path &filepath,
                                     const NDView<T, 1> data) const {

    std::ofstream output_file(filepath / filename);

    if (!output_file) {
        LOG(angcal::TLogLevel::logERROR) << "Error opening file!" << std::endl;
    }

    output_file << std::fixed << std::setprecision(15);

    for (ssize_t i = 0; i < data.size(); ++i) {
        output_file << data[i] << std::endl;
    }
    output_file.close();
}

} // namespace angcal
