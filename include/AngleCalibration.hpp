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
     * @brief Constructor for AngleCalibration class
     * @param mythen_detector_ ptr to MythenDetectorSpecifications storing all
     mythen specific parameters
     * @param flat_field_ ptr to FlatField class storing inverse normalized flat
     field
     * @param mythen_file_reader file reader to read mythen acquisition files
     */
    AngleCalibration(
        std::shared_ptr<MythenDetectorSpecifications> mythen_detector_,
        std::shared_ptr<FlatField> flat_field_,
        std::shared_ptr<MythenFileReader> mythen_file_reader_);

    /** @brief set the histogram bin width [degrees]
     * default '0.0036°'
     */
    void set_histogram_bin_width(double bin_width = 0.0036);

    /**
     * @brief get histogram bin width [degrees]
     */
    double get_histogram_bin_width() const;

    /** @brief number of bins of fixed angle width bin histogram */
    ssize_t num_fixed_angle_width_bins() const;

    /** @brief number of bins covered by base peak region of interest
     * default '101' bins
     */
    ssize_t get_base_peak_ROI_num_bins() const;

    /** @brief set base peak region of interest
     * @param base_peak_roi_ width of base peak ROI [degrees] e.g. [base_peak -
     * base_peak_ROI, base_peak + base_peak_ROI] given in angles defaults to
     * '0.05°'
     * */
    void set_base_peak_ROI_width(const double base_peak_roi_ = 0.05);

    /** @brief half odf the width of base peak region of interest
     * @return width of base peak ROI [degrees] e.g. [base_peak - base_peak_ROI,
     * base_peak + base_peak_ROI] given in angles [degrees]
     * */
    double get_base_peak_ROI_width() const;

    std::shared_ptr<MythenDetectorSpecifications>
    get_detector_specifications() const;

    /** @brief reads the historical Detector Group (DG) parameters from file and
     * transforms them to Best Computing parameters
     * @param filename name of file
     * @param file_reader filereader to read file (default
     InitialAngCalParametersFile:
     * following format module [module_index] center [center] +- [error]
     conversion [conversion] +- [error] offset [offset] +- [error])
     */
    void read_initial_calibration_from_file(
        const std::string &filename,
        const std::shared_ptr<SimpleFileInterface> file_reader =
            std::make_shared<InitialAngCalParametersFile>());

    /**
     * get the histoic DG parameters
     */
    const DGParameters &get_DGparameters() const; // TODO do I need that?

    const BCParameters &get_BCparameters() const;

    /**
     * @brief calibrates the BC (best computing) parameters for all modules
     * @param file_list vector of file_names of acquisition files
     * @param base_peak_angle angle of the selcted base peak [degrees]
     * @param output_file optional, if provided the optimized BC parameters are
     * converted back to DG parameters and are appended to the file (given as
     * full file path)
     */
    void
    calibrate(const std::vector<std::string> &file_list,
              const double base_peak_angle,
              std::optional<std::filesystem::path> output_file = std::nullopt);

    /**
     * @brief calibrates the BC (best computing) parameters for one module
     * @param file_list vector of file_names of acquisition files
     * @param base_peak_angle angle of the selcted base peak given in degrees
     * @param module_index index of module to be calibrated
     */
    void calibrate(const std::vector<std::string> &file_list_,
                   const double base_peak_angle_, const size_t module_index);

    /** @brief check if base peak ROI is contained within module region
     * @param detector_angle: detector position (offset of first strip from
     * default detector position) [degrees]
     * @param bounds_in_angle: boundary that is acceptable [degrees]- per
     * default the base peak ROI width is used (adjustable for debugging) //TODO
     * deprecated
     * @return true if [base_peak - base_peak_ROI, base_peak + base_peak_ROI]
     * is fully contained within the module region
     */
    bool base_peak_is_in_module(
        const size_t module_index, const double detector_angle,
        std::optional<double> bounds_in_angles = std::nullopt) const;

    /**
     * @brief set angle of base peak
     * @param base_peak_angle_ base peak angle [degrees]
     */
    void set_base_peak_angle(const double base_peak_angle_);

    /**
     * @brief get angle of base peak
     * @return base peak angle [degrees]
     */
    double get_base_peak_angle() const;

    /**
     * @brief check if a module only has bad channels
     */
    bool module_is_disconnected(const size_t module_index) const;

    /**
     * @brief set the incident intensity of the first acquisition used for
     * incident intensity correction
     */
    void set_incident_intensity_of_first_acquisition(
        const double incident_intensity);

    /**
     * @brief Performs angular conversion e.g. calculates from raw photon counts
     * the resulting diffraction pattern
     * @param file_list vector of file_names of acquisition files
     * @return flatfield corrected and variance scaled photon counts
     * redistributed to fixed angle width bins given in the range [min_angle,
     * max_angle]
     */
    NDArray<double, 1> convert(const std::vector<std::string> &file_list) const;

    // TODO more useful for debugging
    /**
     * @brief redistribute photon counts to fixed angle width bins for the
     * specified module region
     * @param frame: raw photons from acquisition
     * @param module_index: index of module region to redistribute
     * @return flatfield corrected and variance scaled photon counts in the
     * module region redistributed to fixed angle width bins array covers the
     * range [min_angle, max_angle]
     */
    NDArray<double, 1> redistribute_photon_counts_to_fixed_angle_width_bins(
        const MythenFrame &frame, const size_t module_index) const;

    // TODO more useful for debugging
    /**
     * @brief redistribute photon counts to fixed angle width bins which are
     * within base peak region
     * @param module_index: index of module region to redistribute
     * @return flatfield corrected and variance scaled photon counts in the base
     * peak region redistributed to fixed angle width bins array only covers the
     * base peak region
     */
    NDArray<double, 1> redistributed_photon_counts_in_base_peak_ROI(
        const MythenFrame &frame, const size_t module_index) const;

    /** @brief calculates diffraction angle from DG module parameters (used in
     * Beer's Law)
     * @param detector_angle detector position [degrees]
     * @param strip_index local strip index of module e.g. 0-1279
     * @param distance_to_strip distance to strip [given in strips]
     * @return diffraction angle [degrees]
     */
    double diffraction_angle_from_DG_parameters(
        const size_t module_index, const double detector_angle,
        size_t strip_index, const double distance_to_strip = 0) const;

    /** @brief calculates diffraction angle from BC (best computing) module
     * parameters (used in Beer's Law)
     * @param detector_angle detector position [degrees]
     * @param strip_index local strip index of module e.g. 0-1279
     * @param distance_to_strip distance to strip [given in strips]
     * @return diffraction angle [degrees]
     */
    double diffraction_angle_from_BC_parameters(
        const size_t module_index, const double detector_angle,
        size_t strip_index, const double distance_to_strip = 0) const;

    /** @brief calculates diffraction angle from EE module parameters (used in
     * Beer's Law)
     * @param detector_angle detector position [degrees]
     * @param strip_index local strip index of module
     * @param distance_to_strip distance to strip [given in strips]
     * @return diffraction angle [degrees]
     */
    double diffraction_angle_from_EE_parameters(
        const double module_center_distance, const double normal_distance,
        const double angle, const double detector_angle, size_t strip_index,
        const double distance_to_strip = 0) const;

    /** @brief converts optimized best computing parameters back to DG
     * parameters and writes them to file
     * @param filename full path tooutput file
     */
    void write_to_file(const std::filesystem::path &filename) const;

  private:
    /** @brief calculated the strip width expressed as angle from DG
     * module parameters
     * @param moudle_index index of module
     * @param strip_index local strip index of module e.g. 0-1279
     * @return strip width in angles [degrees]
     */
    double
    angular_strip_width_from_DG_parameters(const size_t module_index,
                                           size_t local_strip_index) const;

    /** @brief calculated the strip width expressed as angle from BC
     * module parameters
     * @param moudle_index index of module
     * @param strip_index local strip index of module e.g. 0-1279
     * @return strip width in angles [degrees]
     */
    double angular_strip_width_from_BC_parameters(
        const size_t module_index, const size_t local_strip_index) const;

    /** @brief calculated the strip width expressed as angle from EE
     * module parameters
     * @param moudle_index index of module
     * @param strip_index local strip index of module e.g. 0-1279
     * @return strip width in angles [degrees]
     */
    double angular_strip_width_from_EE_parameters(
        const double module_center_distance, const double normal_distance,
        const double angle, const size_t local_strip_index) const;

    // TODO deprecated
    /** @brief converts global strip index to local strip index
     */
    size_t global_to_local_strip_index_conversion(
        const size_t global_strip_index) const;

    /**
     * @brief redistributes photon counts around designated region to
     * fixed angle width bins
     * @tparam base_peak_ROI_only: false (redistribute entire module region),
     * true: (only redistribute base peak ROI)
     * @param module_index index of module
     * @param frame data from acquisition (storing detector position and photon
     * counts)
     * @param fixed_angle_width_bin_photon_counts stores redistributed photon
     * counts for fixed angle width bin scaled by variance and flatfield
     * corrected
     * @param inverse_fixed_angle_width_bins_photon_counts_variance stores
     * inverse variance
     * @param sum_statistical_weights sum of statistical weights to normalize
     */
    template <bool base_peak_ROI_only = false>
    void redistribute_photon_counts_to_fixed_angle_width_bins(
        const size_t module_index, const MythenFrame &frame,
        NDView<double, 1> fixed_angle_width_bins_photon_counts,
        NDView<double, 1> inverse_fixed_angle_width_bins_photon_counts_variance,
        NDView<double, 1> sum_statistical_weights) const;

    /**
     * @brief redistributes photon counts to fixed angle width bins around base
     * peak for all acquisitions in file_list and calculates the similarity
     * between the found base peaks regions
     * @brief plot wrapper for gnuplot plot to visualize calibration process
     * (default nullptr)
     * @return similarity of peaks
     */
    double calculate_similarity_of_peaks(const size_t module_index,
                                         PlotHandle gp = nullptr) const;

    /**
     * @brief compares multiple base peak ROIS from different acquisitions and
     * calculate similarity/variance based on goodness_of_fit and weighted
     * average //TODO explain better
     * @param S0 inverse photon_varaince over all runs
     * @param S1 inverse_photon_variance*photon_count over all runs
     * @param S2 inverse_photon_variance*photon_count² over all runs
     * @param num_runs number of frames for which base peak ROI is covered by
     * the respective module
     * @return similarity criterion
     */
    double similarity_criterion(const NDView<double, 1> S0,
                                const NDView<double, 1> S1,
                                const NDView<double, 1> S2,
                                const size_t num_runs) const;

    // TODO: also a bad design - make shift_parameter configurable - consider
    // having second class Optimization
    /**
     * @brief optimizes BC parameters of given module based on similarity
     * criterion
     * @param gp plot wrapper for gnuplot plot to visualize calibration process
     * (default nullptr)
     * @param shift_parameter1 parameter step for parameter angle between module
     * center and module normal (angle_center_module_normal) (default '0.01')
     * @param shift_parameter2 parameter step for parameter distance between
     * module center and sample (module_center_sample_distances) (default
     * '0.005')
     */
    // TODO have PlotHandle as a member change title accordingly
    void optimization_algorithm(const size_t module_index,
                                PlotHandle gp = nullptr,
                                const double shift_parameter1 = 0.01,
                                const double shift_parameter2 = 0.005);

    /**
     * @brief calculated the flatfield corrected photon counts and the error
     * @param photon_counts photon counts
     * @param photon_counts_error error of photon counts (generally variance
     * error is used)
     * @param global_strip_index strip index of photon counts
     * @return pair {corrected photon counts, propagated_error}
     */
    std::pair<double, double>
    flatfield_correction(const double photon_counts,
                         const double photon_counts_error,
                         const size_t global_strip_index) const;

    /**
     * @brief calculate the rate corrected photon counts taking into account the
     * dead time
     * @param photon_counts photon counts
     * @param photon_count_error error of photon counts (generally variance
     * error is used)
     * @return pair {corrected photon counts, propagated_error}
     */
    std::pair<double, double>
    rate_correction(const double photon_counts,
                    const double photon_count_error) const;

    /** @brief calculate the incident intensity corrected photon counts
     * @param photon_counts photon counts
     * @param photon_counts_error error of photon counts (generally variance
     * error is used)
     * @param incident_intensity incident intensity
     * @return pair {corrected photon counts, propagated_error}
     */
    std::pair<double, double>
    incident_intensity_correction(const double photon_counts,
                                  const double photon_counts_error,
                                  const uint64_t incident_intensity) const;

    /** @brief calculate the corrected photon counts (flatfield, rate, incident
     * intensity)
     * @param photon_counts photon counts
     * @param global_strip_index strip index of photon counts
     * @param I0 incident intensity
     * @return pair {corrected photon counts, propagated_error}
     */
    std::pair<double, double>
    photon_count_correction(double photon_counts,
                            const size_t global_strip_index,
                            const uint64_t I0) const;

    /**
     * @brief appends given parameters to file
     */
    void append_to_file(std::ofstream &file, const size_t module_index,
                        const double center, const double conversion,
                        const double offset) const;

    void plot_last_calibration_step(const size_t module_index,
                                    PlotHandle gp) const;

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
     * half the width of region of interest of base peak e.g.
     * [base_peak_angle - base_peak_roi, base_peak_angle +
     * base_peak_roi] [degrees]
     */
    double base_peak_roi_width = 0.05;

    // TODO maybe deprecated - only compute in member function
    ssize_t num_bins{};

    /**
     * center of base peak to use for calibration [degrees]
     */
    double base_peak_angle{};

    /**
     * incident intensity of first acquisition used for incident intensity
     * correction
     */
    std::optional<uint64_t> I0_of_first_acquisition{};

    /**
     * list of acquisition files used for calibration
     */
    std::vector<std::string> file_list{};

    /**
     * file reader to read mythen acquisition files
     */
    std::shared_ptr<MythenFileReader>
        mythen_file_reader{}; // TODO replace by FileInterface ptr
};

template <bool base_peak_ROI_only>
void AngleCalibration::redistribute_photon_counts_to_fixed_angle_width_bins(
    const size_t module_index, const MythenFrame &frame,
    NDView<double, 1> fixed_angle_width_bins_photon_counts,
    NDView<double, 1> inverse_fixed_angle_width_bins_photon_counts_variance,
    NDView<double, 1> sum_statistical_weights) const {

    ssize_t number_of_bins{}; // number of bins using fixed angle width bins
    double left_boundary_roi_base_peak{};  // left boundary of base peak ROI
                                           // [degrees]
    double right_boundary_roi_base_peak{}; // right boundary of base peak ROI
                                           // [degrees]

    if constexpr (base_peak_ROI_only) {
        number_of_bins = get_base_peak_ROI_num_bins(); // fixed angle width bins
                                                       // in base peak ROI
        left_boundary_roi_base_peak = (base_peak_angle - base_peak_roi_width +
                                       0.5 * histogram_bin_width); // in degrees
        right_boundary_roi_base_peak =
            (base_peak_angle + base_peak_roi_width -
             0.5 * histogram_bin_width); // in degrees
    } else {
        number_of_bins = num_fixed_angle_width_bins();
        left_boundary_roi_base_peak =
            mythen_detector->min_angle(); // dummy values
        right_boundary_roi_base_peak =
            mythen_detector->max_angle(); // dummy values
    }

    // iterate over all strips of module
    for (size_t strip_index = 0;
         strip_index < mythen_detector->strips_per_module(); ++strip_index) {

        size_t global_strip_index =
            module_index * mythen_detector->strips_per_module() +
            strip_index; // TODO: is this really correct - check sign

        if (mythen_detector->get_bad_channels()(global_strip_index)) {
            continue; // skip bad channels
        }

        double left_strip_boundary_angle = diffraction_angle_from_BC_parameters(
            module_index, frame.detector_angle, strip_index,
            -0.5); // left strip boundary in angles [degrees]

        double right_strip_boundary_angle =
            diffraction_angle_from_BC_parameters(
                module_index, frame.detector_angle, strip_index,
                +0.5); // right strip boundary in angles [degrees]

        // skip strip if not in ROI
        if constexpr (base_peak_ROI_only) {
            if (left_strip_boundary_angle > right_boundary_roi_base_peak ||
                right_strip_boundary_angle < left_boundary_roi_base_peak) {
                continue;
            }
        }

        // get corrected photon counts and propagate error
        auto [corrected_photon_counts, corrected_photon_counts_variance] =
            photon_count_correction(frame.photon_counts(global_strip_index) + 1,
                                    global_strip_index,
                                    frame.incident_intensity);

        double inverse_corrected_photon_counts_variance =
            1.0 / corrected_photon_counts_variance;

        double strip_width_angle =
            std::abs(right_strip_boundary_angle -
                     left_strip_boundary_angle); // strip width in angles
        // angular_strip_width_from_BC_parameters(module_index, strip_index);

        // redistribute to fixed angle width bins
        double photon_counts_per_bin =
            corrected_photon_counts * histogram_bin_width / strip_width_angle;

        double inverse_photon_counts_variance_per_bin =
            inverse_corrected_photon_counts_variance *
            std::pow(strip_width_angle / histogram_bin_width,
                     2); // Var(aX) = a^2 Var(X)

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
            0; // the computed bin indices dont start at zero but are relative
               // to detector angle and strip index

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
                        (base_peak_angle - base_peak_roi_width) /
                        histogram_bin_width)); // bin index starts at zero
            } else {
                proper_bin_index =
                    bin_index -
                    static_cast<ssize_t>(
                        mythen_detector->min_angle() /
                        histogram_bin_width); // bin index starts at zero
            }

            double statistical_weight =
                bin_coverage_factor * inverse_photon_counts_variance_per_bin;

            // std::cout << "statistical weight: " << statistical_weight <<
            // "\n";

            fixed_angle_width_bins_photon_counts(proper_bin_index) +=
                statistical_weight * photon_counts_per_bin;

            inverse_fixed_angle_width_bins_photon_counts_variance(
                proper_bin_index) +=
                inverse_photon_counts_variance_per_bin *
                bin_coverage_factor; // statistical_weigth² *
                                     // corrected_photon_counts_variance

            sum_statistical_weights(proper_bin_index) += statistical_weight;
        }
    }
}

} // namespace angcal
