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
#include "MythenDetector.hpp"
#include "MythenFileReader.hpp"
#include "Parameters.hpp"
#include "PlotCalibrationProcess.hpp"
#include "aare/NDArray.hpp"
#include "helpers/FileInterface.hpp"
#include "helpers/LogFiles.hpp"
#include "helpers/custom_errors.hpp"

using namespace aare;

namespace angcal {

class AngleCalibration {

  public:
    /**
     * @brief Constructor for AngleCalibration class
     * @param mythen_detector_ ptr to MythenDetector storing all mythen specific
     parameters and methods to calculate diffraction angles
     * @param flat_field_ ptr to FlatField class storing inverse normalized flat
     field
     */
    AngleCalibration(std::shared_ptr<MythenDetector> mythen_detector_,
                     std::shared_ptr<FlatField> flat_field_);

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

    /// @brief get base peak angle
    /// @return angle of base peak [degrees]
    double get_base_peak_angle() const;

    /// @brief set base peak angle
    /// @param base_peak_angle_ angle of base peak [degrees]
    void set_base_peak_angle(const double base_peak_angle_);

    /// @brief get detector geometry used for angular conversion
    /// @return shared pointer to MythenDetector
    const std::shared_ptr<MythenDetector> get_mythen_detector() const {
        return m_mythen_detector;
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
            std::make_shared<CustomBadChannelsFile>());

    NDView<bool, 1> get_bad_channels();

    /**
     * @brief set bad channels
     * @param bad_channels_ boolean NDView of size of total number of
     * strips/channels in detector, stores 'TRUE' if strip is a bad channel
     * otherwise 'FALSE'
     */
    void set_bad_channels(const NDView<bool, 1> bad_channels_);

    /**
     * @brief get the current DG parameters
     */
    const DGParameters &get_DGparameters() const;

    /**
     * @brief get the current BC parameters
     */
    const BCParameters &get_BCparameters() const;

    /**
     * @brief calibrates the BC (best computing) parameters for all modules
     * @tparam visualize_calibration if true visualize calibration process
     * (default: false)
     * @param frames 2D NDView of frames (rows: acquisitions, cols: channels) -
     * need to be sorted by detector angle
     * @param detector_angles 1D NDView of monitor positions for acquisitions
     * @param base_peak_angle angle of the selcted base peak [degrees]
     */
    template <bool visualize_calibration = false>
    void calibrate(const NDView<uint32_t, 2> frames,
                   const NDView<double, 1> detector_angles,
                   const double base_peak_angle);

    /**
     * @brief calibrates the BC (best computing) parameters for one module
     * @tparam visualize_calibration if true visualize calibration process
     * (default: false)
     * @param frames 2D NDView of frames (rows: acquisitions, cols: channels) -
     * need to be sorted by detector angle
     * @param detector_angles 1D NDView of monitor positions for the
     * acquisitions
     * @param base_peak_angle angle of the selcted base peak given in degrees
     * @param module_index index of module to be calibrated
     */
    template <bool visualize_calibration = false>
    void calibrate(const NDView<uint32_t, 2> frames,
                   const NDView<double, 1> detector_angles,
                   const double base_peak_angle_, const size_t module_index);

    /** @brief check if base peak ROI is contained within module region
     * @param detector_angle: detector position (offset of first strip from
     * default detector position) [degrees]
     * @return true if [base_peak - base_peak_ROI, base_peak + base_peak_ROI]
     * is fully contained within the module region
     */
    bool base_peak_is_in_module(const size_t module_index,
                                const double detector_angle) const;

    /**
     * @brief check if a module only has bad channels
     */
    bool module_is_disconnected(const size_t module_index) const;

    /**
     * @brief sets the angular range for the diffraction pattern e.g.
     * diffraction pattern calculated for [min_angle, max_angle]
     * @param min_angle minimum angle [degrees]
     * @param max_angle maximum angle [degrees]
     */
    void set_angular_range(const double min_angle, const double max_angle);

    /**
     * @brief gets the angular range for the diffraction pattern e.g.
     * diffraction pattern calculated for [min_angle, max_angle]
     * @return pair of (min_angle, max_angle) [degrees]
     */
    std::pair<double, double> get_angular_range() const;

  private:
    /**
     * @brief calculates the frames for which the base peak ROI is covered by
     * the respective module region
     * @param frames 2D NDView of frames (rows: acquisitions, cols: channels) -
     * need to be sorted by detector angle
     * @param detector_angles 1D NDView of monitor positions for acquisition
     * @param module_index index of module
     */
    std::tuple<NDView<uint32_t, 2>, NDView<double, 1>>
    get_frames_with_base_peak_overlap(const NDView<uint32_t, 2> frames,
                                      const NDView<double, 1> detector_angles,
                                      const size_t module_index) const;

    /**
     * @brief applies flatfield correction for a given module and redistributes
     * them to fixed angle width bins
     * @param module_index index of module
     * @param frame data from acquisition (storing detector position and photon
     * counts)
     * @param weighted_moments weighted photon counts of order 0-2 eg. order k:
     * sum of photon_count_i^k*weight_i counts for fixed angle width bin and the
     * corresponding variance
     * @param num_contributing_photons number of photons contributing to the
     * fixed angle width bin
     * @param detector_angle detector axis motor position [degrees]
     */
    void process_photon_counts(const size_t module_index,
                               const NDView<uint32_t, 1> frame,
                               const double detector_angle,
                               NDView<double, 2> weighted_moments,
                               NDView<size_t, 1> num_contributing_photons);

    /**
     * @brief redistributes photon counts to fixed angle width bins around base
     * peak for all acquisitions in passed frames and calculates the similarity
     * between the found base peaks regions
     * @param frames_with_base_peak_overlap vector of frames from acquisitions
     * for which the base peak ROI is covered by the respective module region
     * @param detector_angles_with_base_peak_overlap 1D NDView of monitor
     * positions for acquisitions
     * @brief plot wrapper to visualize calibration process
     * (default nullptr)
     * @return similarity of peaks
     */
    template <bool visualize_calibration = false>
    double calculate_similarity_of_peaks_between_acquisitions(
        const size_t module_index,
        const NDView<uint32_t, 2> frames_with_base_peak_overlap,
        const NDView<double, 1> detector_angles_with_base_peak_overlap,
        std::shared_ptr<PlotCalibrationProcess> plot = nullptr);

    /**
     * @brief redistributes photon counts to fixed angle width bins around base
     * peak for all acquisitions in passed frames and calculates the similarity
     * between the found base peaks regions of two modules
     * @param frames_with_base_peak_overlap_module vector of frames from
     * acquisitions for which the base peak ROI is covered by the respective
     * module region
     * @param detector_angles_with_base_peak_overlap_module 1D NDView of monitor
     * positions for acquisitions for the module
     * @param frames_with_base_peak_overlap_prev_module vector of frames from
     * acquisitions for which the base peak ROI is covered by the respective
     * previous module region
     * @param detector_angles_with_base_peak_overlap_prev_module 1D NDView of
     * monitor positions for acquisitions for the previous module
     * @brief plot wrapper to visualize calibration process
     * (default nullptr)
     * @return similarity of peaks
     */
    template <bool visualize_calibration = false>
    double calculate_similarity_of_peaks_between_modules(
        const size_t module_index,
        const NDView<uint32_t, 2> frames_with_base_peak_overlap_module,
        const NDView<double, 1> detector_angles_with_base_peak_overlap_module,
        const NDView<uint32_t, 2> frames_with_base_peak_overlap_prev_module,
        const NDView<double, 1>
            detector_angles_with_base_peak_overlap_prev_module,
        std::shared_ptr<PlotCalibrationProcess> plot = nullptr);

    /**
     * @brief calibrates distance to module center and angle of module center
     * and normal (L and delta) for a specific module
     * @tparam visualize_calibration if true visualize calibration process
     * (default: false)
     * @param module_index index of module to be calibrated
     * @param frames_width_base_peak_overlap frames with base peak overlap for
     * the module
     * @param detector_angles 1D NDView of monitor positions for acquisitions
     */
    template <bool visualize_calibration = false>
    void calibrate_coupled_parameters(
        const size_t module_index,
        const NDView<uint32_t, 2> frames_width_base_peak_overlap,
        const NDView<double, 1> detector_angles);

    /**
     * @brief calibrates the angle between center of module and beam (phi) for a
     * specific module
     * @tparam visualize_calibration if true visualize calibration process
     * (default: false)
     * @param module_index index of module to be calibrated
     * @param frames_with_base_peak_overlap frames with base peak overlap module
     * to be calibrated
     * @param detector_angles_with_base_peak_overlap 1D NDView of monitor
     * positions for acquisitions
     * @param frames_with_base_peak_overlap_prev_module frames with base peak
     * overlap prev module
     * @param detector_angles_with_base_peak_overlap_prev_module 1D NDView of
     * monitor positions for acquisitions for the previous module
     */
    template <bool visualize_calibration = false>
    void calibrate_offset(
        const size_t module_index,
        const NDView<uint32_t, 2> frames_with_base_peak_overlap,
        const NDView<double, 1> detector_angles_with_base_peak_overlap,
        const NDView<uint32_t, 2> frames_with_base_peak_overlap_prev_module,
        const NDView<double, 1>
            detector_angles_with_base_peak_overlap_prev_module);

    /**
     * @brief calculates the average angle of the base peak center over several
     * acquisitions for a given module index
     * @param module_index index of module to calculate base peak center for
     * @param frames_with_base_peak_overlap frames with base peak overlap for
     * the given module index
     * @param detector_angles_with_base_peak_overlap 1D NDView of monitor
     * positions for acquisitions
     * @return average angle of base peak center [degrees]
     */
    double center_base_peak(
        const size_t module_index,
        const NDView<uint32_t, 2> frames_with_base_peak_overlap,
        const NDView<double, 1> detector_angles_with_base_peak_overlap);

    /**
     * @brief compares multiple base peak ROIS from different acquisitions and
     * calculate error based on goodness_of_fit and weighted
     * average
     * @param weighted_sums weighted photon counts of order 0-2 eg. order k: sum
     * of photon_count_i^k*weight_i
     * @param num_contributing_photons number of photons contributing to the
     * fixed angle width bin
     * @return similarity criterion
     */
    double chi_similarity_criterion(
        const NDView<double, 2> weighted_sums,
        const NDView<size_t, 1> num_contributing_photons) const;

    /**
     * @brief averages redistributed photon counts and calculates error
     * based on weighted average and goodness of fit
     * @param weighted_sums weighted photon counts of order 0-2 eg. order k: sum
     * of photon_count_i^k*weight_i
     * @param num_contributing_photons number of photons contributing to the
     * fixed angle width bin
     * @return ndarray of size (num_bins_in_ROI, 2) where first dimension
     * corresponds to average redistributed photon counts and second dimension
     * corresponds to error of average redistributed photon counts
     */
    NDArray<double, 2>
    computeErrorandMean(const NDView<double, 2> weighted_sums,
                        const NDView<size_t, 1> num_contributing_photons) const;

    /**
     * @brief optimizes module center distance L and angle center module
     * normal delta of BC parameters of given module based on chi similarity
     * criterion
     * @param module_index index of module to be calibrated
     * @param frames_with_base_peak_overlap vector of frames from
     * acquisitions for which the base peak ROI is covered by the respective
     * module region
     * @param detector_angles_with_base_peak_overlap 1D NDView of monitor
     * positions for acquisitions
     * @param plot plot wrapper to visualize calibration process
     * (default nullptr)
     * @param delta_parameter1 parameter step for parameter module center
     * distance L (default '0.5')
     * @param delta_parameter2 parameter step for parameter angle between
     * center of module and module normal (delta) (default '1.0e-6')
     */
    void optimize_coupled_parameters(
        const size_t module_index,
        const NDView<uint32_t, 2> frames_with_base_peak_overlap,
        const NDView<double, 1> detector_angles_with_base_peak_overlap,
        std::shared_ptr<PlotCalibrationProcess> plot = nullptr,
        const double delta_parameter1 = 0.5,
        const double delta_parameter2 = 1.0e-6);

    /**
     * @brief optimizes angle between center of module and beam (\psi) of BC
     * parameters of given module based on chi similarity criterion and center
     * of base peak
     * @param module_index index of module to be calibrated
     * @param frames_with_base_peak_overlap vector of frames with base peak
     * overlap for module
     * @param detector_angles_with_base_peak_overlap 1D NDView of monitor
     * positions for acquisitions for the module
     * @param frames_with_base_peak_overlap_prev_module vector of frames with
     * base peak overlap for previous module
     * @param detector_angles_with_base_peak_overlap_prev_module 1D NDView of
     * monitor positions for acquisitions for the previous module
     * @param plot plot wrapper to visualize calibration process
     * (default nullptr)
     * @param delta_parameter parameter step for angle between center of module
     * and beam (default '0.005')
     */
    void optimize_offset_parameter(
        const size_t module_index,
        const NDView<uint32_t, 2> frames_with_base_peak_overlap,
        const NDView<double, 1> detector_angles_with_base_peak_overlap,
        const NDView<uint32_t, 2> frames_with_base_peak_overlap_prev_module,
        const NDView<double, 1>
            detector_angles_with_base_peak_overlap_prev_module,
        std::shared_ptr<PlotCalibrationProcess> plot = nullptr,
        double delta_parameter = 0.005);

    /*
    void plot_calibration_step_coupled_parameters(
        const size_t module_index,
        const std::vector<MythenFrame> &frames_with_base_peak_overlap,
        PlotHandle gp);
    */

    /*
    void plot_calibration_step_offset_parameter(
        const size_t module_index,
        const std::vector<MythenFrame> &frames_with_base_peak_overlap_module,
        const std::vector<MythenFrame>
            &frames_with_base_peak_overlap_prev_module,
        PlotHandle gp);
    */

  private:
    /// @brief Array of size strips/channels in detector, stores 'TRUE' if strip
    /// is a bad channel otherwise 'FALSE'
    NDArray<bool, 1> bad_channels{};

    /**
     * detector geometry
     */
    std::shared_ptr<MythenDetector> m_mythen_detector{};

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

    /**
     * center of base peak to use for calibration [degrees]
     */
    double base_peak_angle{};

    /// @brief  @brief minimum angle for angular conversion [degrees]
    double m_min_angle = -180.0;

    /// @brief  maximum angle for angular conversion [degrees]
    double m_max_angle = 180.0;
};

template <bool visualize_calibration>
double AngleCalibration::calculate_similarity_of_peaks_between_acquisitions(
    const size_t module_index,
    const NDView<uint32_t, 2> frames_with_base_peak_overlap,
    const NDView<double, 1> detector_angles_bp_overlap,
    std::shared_ptr<PlotCalibrationProcess> plot) {

    ssize_t num_bins_in_ROI = get_base_peak_ROI_num_bins();

    NDArray<double, 2> total_weighted_sums(
        std::array<ssize_t, 2>{num_bins_in_ROI, 3},
        -1.0); // // weighted moments of order 0-2 for each bin

    NDArray<size_t, 1> num_contributing_photons(
        std::array<ssize_t, 1>{num_bins_in_ROI},
        0); // number of contributing photons for each bin

    size_t num_runs{};
    for (ssize_t frame_index = 0;
         frame_index < frames_with_base_peak_overlap.shape(0); ++frame_index) {

        // base peak angle is in module
        if (base_peak_is_in_module(module_index,
                                   detector_angles_bp_overlap[frame_index])) {

            LOG(TLogLevel::logDEBUG2) << fmt::format(
                "detector angle: {}", detector_angles_bp_overlap[frame_index]);

            if constexpr (visualize_calibration) {
                NDArray<double, 2> weighted_sums(
                    std::array<ssize_t, 2>{num_bins_in_ROI, 3},
                    -1.0); // // weighted moments of order 0-2 for each bin

                NDArray<double, 2> fixed_angle_width_bins_photon_counts(
                    std::array<ssize_t, 2>{num_bins_in_ROI, 2}, -1.0);

                process_photon_counts(
                    module_index,
                    frames_with_base_peak_overlap.sub_view(frame_index),
                    detector_angles_bp_overlap[frame_index],
                    weighted_sums.view(), num_contributing_photons.view());

                for (ssize_t i = 0; i < weighted_sums.shape(0); ++i) {

                    if (weighted_sums(i, 0) == -1.0) {
                        continue;
                    } else {

                        if (total_weighted_sums(i, 0) == -1.0) {
                            total_weighted_sums(i, 0) = 0.0;
                            total_weighted_sums(i, 1) = 0.0;
                            total_weighted_sums(i, 2) = 0.0;
                        }

                        total_weighted_sums(i, 0) += weighted_sums(i, 0);
                        total_weighted_sums(i, 1) += weighted_sums(i, 1);
                        total_weighted_sums(i, 2) += weighted_sums(i, 2);

                        fixed_angle_width_bins_photon_counts(i, 0) =
                            weighted_sums(i, 0) <
                                    std::numeric_limits<double>::epsilon()
                                ? 0.0
                                : weighted_sums(i, 1) /
                                      weighted_sums(
                                          i,
                                          0); // TODO: two data for now easier
                    }
                }

                plot->add_curve(fixed_angle_width_bins_photon_counts.view());
            } else {

                process_photon_counts(
                    module_index,
                    frames_with_base_peak_overlap.sub_view(frame_index),
                    detector_angles_bp_overlap[frame_index],
                    total_weighted_sums.view(),
                    num_contributing_photons.view());
            }

            ++num_runs;
        }
    }

    LOG(TLogLevel::logDEBUG1)
        << fmt::format("num runs with base peak overlap: {}", num_runs);

    if (num_runs == 0) {
        throw NoBasePeakOverLapError();
    }

    double similarity_of_peaks =
        num_runs == 0
            ? std::numeric_limits<double>::infinity()
            : chi_similarity_criterion(total_weighted_sums.view(),
                                       num_contributing_photons.view());

    LOG(TLogLevel::logDEBUG1) << "similarity_of_peaks: " << similarity_of_peaks;

    if (std::isnan(similarity_of_peaks) || std::isinf(similarity_of_peaks)) {
        throw std::runtime_error(LOCATION +
                                 "Similarity of peaks is NaN or "
                                 "infinity - something went wrong in "
                                 "calculation of similarity of peaks");
    }

    return similarity_of_peaks;
}

template <bool visualize_calibration>
double AngleCalibration::calculate_similarity_of_peaks_between_modules(
    const size_t module_index,
    const NDView<uint32_t, 2> frames_bp_overlap_module,
    const NDView<double, 1> detector_angles_bp_overlap_module,
    const NDView<uint32_t, 2> frames_bp_overlap_prev_module,
    const NDView<double, 1> detector_angles_bp_overlap_prev_module,
    std::shared_ptr<PlotCalibrationProcess> plot) {

    if (module_index == 0) {
        throw std::runtime_error(LOCATION +
                                 "Module index is zero - offset of module 0 is "
                                 "fixed to zero and cannot be optimized");
    }

    ssize_t num_bins_in_ROI = get_base_peak_ROI_num_bins();

    NDArray<double, 2> total_weighted_sums(
        std::array<ssize_t, 2>{num_bins_in_ROI, 3},
        -1.0); // weighted moments of order 0-2 for each bin and module

    NDArray<size_t, 1> num_contributing_photons(
        std::array<ssize_t, 1>{num_bins_in_ROI},
        0); // number of contributing photons for each bin

    size_t num_modules_with_base_peak_overlap = 2;

    for (const size_t mod_index : {module_index - 1, module_index}) {

        size_t num_runs{};

        const auto &frames_with_base_peak_overlap =
            mod_index == module_index ? frames_bp_overlap_module
                                      : frames_bp_overlap_prev_module;

        const auto &detector_angles_with_base_peak_overlap =
            mod_index == module_index ? detector_angles_bp_overlap_module
                                      : detector_angles_bp_overlap_prev_module;

        if (!module_is_disconnected(mod_index)) {

            for (ssize_t frame_index = 0;
                 frame_index < frames_with_base_peak_overlap.shape(0);
                 ++frame_index) {

                if (base_peak_is_in_module(
                        mod_index,
                        detector_angles_with_base_peak_overlap[frame_index])) {

                    if constexpr (visualize_calibration) {
                        NDArray<double, 2> weighted_sums(
                            std::array<ssize_t, 2>{num_bins_in_ROI, 3},
                            -1.0); // // weighted moments of order 0-2 for each
                                   // bin

                        NDArray<double, 2> fixed_angle_width_bins_photon_counts(
                            std::array<ssize_t, 2>{num_bins_in_ROI, 2}, -1.0);

                        process_photon_counts(
                            mod_index,
                            frames_with_base_peak_overlap.sub_view(frame_index),
                            detector_angles_with_base_peak_overlap[frame_index],
                            weighted_sums.view(),
                            num_contributing_photons.view());

                        for (ssize_t i = 0; i < weighted_sums.shape(0); ++i) {

                            if (weighted_sums(i, 0) == -1.0) {
                                continue;
                            } else {

                                if (total_weighted_sums(i, 0) == -1.0) {
                                    total_weighted_sums(i, 0) = 0.0;
                                    total_weighted_sums(i, 1) = 0.0;
                                    total_weighted_sums(i, 2) = 0.0;
                                }

                                total_weighted_sums(i, 0) +=
                                    weighted_sums(i, 0);
                                total_weighted_sums(i, 1) +=
                                    weighted_sums(i, 1);
                                total_weighted_sums(i, 2) +=
                                    weighted_sums(i, 2);

                                fixed_angle_width_bins_photon_counts(i, 0) =
                                    weighted_sums(i, 0) < std::numeric_limits<
                                                              double>::epsilon()
                                        ? 0.0
                                        : weighted_sums(i, 1) /
                                              weighted_sums(
                                                  i,
                                                  0); // TODO: two data for now
                                                      // easier
                            }
                        }

                        plot->add_curve(
                            fixed_angle_width_bins_photon_counts.view());
                    } else {
                        process_photon_counts(
                            mod_index,
                            frames_with_base_peak_overlap.sub_view(frame_index),
                            detector_angles_with_base_peak_overlap[frame_index],
                            total_weighted_sums.view(),
                            num_contributing_photons.view());
                    }
                    ++num_runs;
                }
            }
            if (num_runs == 0) {
                num_modules_with_base_peak_overlap--;
                throw NoBasePeakOverLapError();
            }
        } else {
            num_modules_with_base_peak_overlap--;
            LOG(TLogLevel::logINFO) << fmt::format(
                "Module {} is disconnected - skipping module", module_index);
        }
    }

    if (num_modules_with_base_peak_overlap != 2) {
        throw std::runtime_error(
            LOCATION + "Base peak does not overlap with both modules - cannot "
                       "calculate similarity of peaks between modules");
    }

    double similarity_of_peaks = chi_similarity_criterion(
        total_weighted_sums.view(), num_contributing_photons.view());

    if (std::isnan(similarity_of_peaks) || std::isinf(similarity_of_peaks)) {
        throw std::runtime_error(LOCATION +
                                 "Similarity of peaks is NaN or "
                                 "infinity - something went wrong in "
                                 "calculation of similarity of peaks");
    }

    return similarity_of_peaks;
}

template <bool visualize_calibration>
void AngleCalibration::calibrate_offset(
    const size_t module_index, const NDView<uint32_t, 2> frames_bp_overlap,
    const NDView<double, 1> detector_angles_bp_overlap,
    const NDView<uint32_t, 2> frames_bp_overlap_prev_module,
    const NDView<double, 1> detector_angles_bp_overlap_prev_module) {

    if (module_index == 0) {
        throw std::runtime_error(LOCATION +
                                 "Module index is zero - offset of module 0 is "
                                 "fixed to zero and cannot be optimized");
    }

    // skip if module is not connected
    if (module_is_disconnected(module_index) or
        module_is_disconnected(module_index - 1)) {
        LOG(TLogLevel::logINFO)
            << fmt::format("module {} or module {} is disconnected",
                           module_index, module_index - 1);
        return;
    }

    std::shared_ptr<PlotCalibrationProcess> plot = nullptr;
    if constexpr (visualize_calibration) {
        std::string plot_title =
            fmt::format("Base Peaks between modules {} and {} ",
                        module_index - 1, module_index);
        plot = std::make_shared<PlotCalibrationProcess>(this, plot_title);
        plot->initializematplotlib();
    }

    try {
        optimize_offset_parameter(
            module_index, frames_bp_overlap, detector_angles_bp_overlap,
            frames_bp_overlap_prev_module,
            detector_angles_bp_overlap_prev_module, plot, 1.0e-3);
    } catch (const NoBasePeakOverLapError &e) {
        LOG(angcal::TLogLevel::logINFO)
            << e.what() << fmt::format(" Skipping module {}.", module_index);
        return;
    }

    LOG(TLogLevel::logDEBUG)
        << fmt::format("calibration for module {} finished", module_index);
}

template <bool visualize_calibration>
void AngleCalibration::calibrate_coupled_parameters(
    const size_t module_index, const NDView<uint32_t, 2> frames,
    const NDView<double, 1> detector_angles) {

    // skip if module is not connected

    LOG(angcal::TLogLevel::logINFO)
        << "starting calibration for module " << module_index;

    std::shared_ptr<PlotCalibrationProcess> plot = nullptr;
    if constexpr (visualize_calibration) {
        std::string plot_title =
            fmt::format("Base Peaks for module {} ", module_index);
        plot = std::make_shared<PlotCalibrationProcess>(this, plot_title);

        plot->initializematplotlib();
    }

    try {
        optimize_coupled_parameters(module_index, frames, detector_angles, plot,
                                    0.5, 1.0e-6);
    } catch (const NoBasePeakOverLapError &e) {
        LOG(angcal::TLogLevel::logINFO)
            << e.what() << fmt::format(" Skipping module {}.", module_index);
        return;
    }

    LOG(TLogLevel::logDEBUG)
        << fmt::format("calibration for module {} finished", module_index);
}

template <bool visualize_calibration>
void AngleCalibration::calibrate(const NDView<uint32_t, 2> frames,
                                 const NDView<double, 1> detector_angles,
                                 const double base_peak_angle_) {

    base_peak_angle = base_peak_angle_;

    NDView<uint32_t, 2> prev_frames_bp_overlap;
    NDView<double, 1> prev_detector_angles_bp_overlap;
    NDView<uint32_t, 2> frames_bp_overlap;
    NDView<double, 1> detector_angles_bp_overlap;

    try {
        std::tie(prev_frames_bp_overlap, prev_detector_angles_bp_overlap) =
            get_frames_with_base_peak_overlap(frames, detector_angles, 0);
    } catch (std::runtime_error &e) {

        throw std::runtime_error(
            std::string(e.what()) +
            fmt::format(" No base peak overlap for module {}.", 0));
    }

    calibrate_coupled_parameters<visualize_calibration>(
        0, prev_frames_bp_overlap, prev_detector_angles_bp_overlap);

    for (size_t module_idx = 1; module_idx < m_mythen_detector->num_modules();
         ++module_idx) {

        try {
            std::tie(frames_bp_overlap, detector_angles_bp_overlap) =
                get_frames_with_base_peak_overlap(frames, detector_angles,
                                                  module_idx);
        } catch (std::runtime_error &e) {

            LOG(angcal::TLogLevel::logINFO)
                << std::string(e.what()) +
                       fmt::format(" No base peak overlap for module {}.",
                                   module_idx);
            continue; // skip module if no base peak overlap
        }

        calibrate_coupled_parameters<visualize_calibration>(
            module_idx, frames_bp_overlap, detector_angles_bp_overlap);

        calibrate_offset<visualize_calibration>(
            module_idx, frames_bp_overlap, detector_angles_bp_overlap,
            prev_frames_bp_overlap, prev_detector_angles_bp_overlap);

        prev_frames_bp_overlap = frames_bp_overlap;
        prev_detector_angles_bp_overlap = detector_angles_bp_overlap;
    }

    PlotCalibrationProcess::kill_python_interpreter();
}

template <bool visualize_calibration>
void AngleCalibration::calibrate(const NDView<uint32_t, 2> frames,
                                 const NDView<double, 1> detector_angles,
                                 const double base_peak_angle_,
                                 const size_t module_index) {

    if (module_is_disconnected(module_index) ||
        (module_index > 0 && module_is_disconnected(module_index - 1))) {
        throw std::runtime_error(
            LOCATION +
            fmt::format("Module {} or module {} is disconnected - cannot "
                        "calibrate module {}",
                        module_index, module_index > 0 ? module_index - 1 : 0,
                        module_index));
    }
    base_peak_angle = base_peak_angle_;

    NDView<uint32_t, 2> frames_with_base_peak_overlap;
    NDView<double, 1> detector_angles_with_base_peak_overlap;
    NDView<uint32_t, 2> frames_with_base_peak_overlap_prev_module;
    NDView<double, 1> detector_angles_with_base_peak_overlap_prev_module;

    try {
        std::tie(frames_with_base_peak_overlap,
                 detector_angles_with_base_peak_overlap) =
            get_frames_with_base_peak_overlap(frames, detector_angles,
                                              module_index);
    } catch (std::runtime_error &e) {
        throw std::runtime_error(
            std::string(e.what()) +
            fmt::format(" No base peak overlap for module {}.", module_index));
    }

    const auto t0 = std::chrono::steady_clock::now();

    calibrate_coupled_parameters<visualize_calibration>(
        module_index, frames_with_base_peak_overlap,
        detector_angles_with_base_peak_overlap);

    if (module_index > 0) {

        try {
            std::tie(frames_with_base_peak_overlap_prev_module,
                     detector_angles_with_base_peak_overlap_prev_module) =
                get_frames_with_base_peak_overlap(frames, detector_angles,
                                                  module_index - 1);
        } catch (std::runtime_error &e) {
            throw std::runtime_error(
                std::string(e.what()) +
                fmt::format(" No base peak overlap for module {}.",
                            module_index - 1));
        }

        calibrate_coupled_parameters<visualize_calibration>(
            module_index - 1, frames_with_base_peak_overlap_prev_module,
            detector_angles_with_base_peak_overlap_prev_module);

        calibrate_offset<visualize_calibration>(
            module_index, frames_with_base_peak_overlap,
            detector_angles_with_base_peak_overlap,
            frames_with_base_peak_overlap_prev_module,
            detector_angles_with_base_peak_overlap_prev_module);
    }

    const auto t1 = std::chrono::steady_clock::now();
    const auto elapsed_time =
        std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    LOG(TLogLevel::logINFO)
        << fmt::format("calibration for module {} took {} milliseconds",
                       module_index, elapsed_time);

    PlotCalibrationProcess::kill_python_interpreter();
}

} // namespace angcal
