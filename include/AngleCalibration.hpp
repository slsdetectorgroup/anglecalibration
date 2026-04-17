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
#include "PlotCalibrationProcess.hpp"
#include "aare/NDArray.hpp"
#include "helpers/FileInterface.hpp"
// #include "helpers/LogFiles.hpp"
#include "helpers/custom_errors.hpp"

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

    // void plot_all_base_peaks(PlotHandle gp);

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

    /// @brief set vector of acquisitions files for calibration or conversion
    /// @param file_list
    void set_calibration_files(const std::vector<std::string> &file_list);

    /// @brief get vector of acquisition files for calibration or conversion
    /// @return vector of file names
    const std::vector<std::string> &get_calibration_files() const;

    /// @brief  get detector specifications
    std::shared_ptr<MythenDetectorSpecifications>
    get_detector_specifications() const;

    /// @brief get file reader to read mythen acquisition files
    /// @return shared pointer to MythenFileReader
    const std::shared_ptr<MythenFileReader> get_file_reader() const {
        return mythen_file_reader;
    }

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

    NDView<bool, 1> get_bad_channels() const;

    /**
     * @brief set bad channels
     * @param bad_channels_ boolean NDArray of size of total number of
     * strips/channels in detector, stores 'TRUE' if strip is a bad channel
     * otherwise 'FALSE'
     */
    void set_bad_channels(const NDArray<bool, 1> &bad_channels_);

    /**
     * @brief get the histoic DG parameters
     */
    const DGParameters &get_DGparameters() const; // TODO do I need that?

    /**
     * @brief get the current BC parameters
     */
    const BCParameters &get_BCparameters() const;

    /**
     * @brief calibrates the BC (best computing) parameters for all modules
     * @tparam visualize_calibration if true visualize calibration process
     * (default: false)
     * @param file_list vector of file_names of acquisition files
     * @param base_peak_angle angle of the selcted base peak [degrees]
     * @param output_file optional, if provided the optimized BC parameters are
     * converted back to DG parameters and are appended to the file (given as
     * full file path)
     */
    template <bool visualize_calibration = false>
    void
    calibrate(const std::vector<std::string> &file_list,
              const double base_peak_angle,
              std::optional<std::filesystem::path> output_file = std::nullopt);

    /**
     * @brief calibrates distance to module center and angle of module center
     * and normal (L and delta) for all modules
     * @tparam visualize_calibration if true visualize calibration process
     * (default: false)
     * @param detector_angles vector of monitor positions for acquisition files
     */
    template <bool visualize_calibration = false>
    void
    calibrate_coupled_parameters(const std::vector<double> &detector_angles);

    /**
     * @brief calibrates distance to module center and angle of module center
     * and normal (L and delta) for a specific module
     * @tparam visualize_calibration if true visualize calibration process
     * (default: false)
     * @param module_index index of module to be calibrated
     * @param frames_width_base_peak_overlap frames with base peak overlap for
     * the module
     */
    template <bool visualize_calibration = false>
    void calibrate_coupled_parameters(
        const size_t module_index,
        const std::vector<MythenFrame> &frames_width_base_peak_overlap);

    /**
     * @brief calibrates the angle between center of module and beam (phi) for
     * all modules
     * @tparam visualize_calibration if true visualize calibration process
     * (default: false)
     * @param detector_angles vector of monitor positions for acquisition files
     */
    template <bool visualize_calibration = false>
    void calibrate_offset(const std::vector<double> &detector_angles);

    /**
     * @brief calibrates the angle between center of module and beam (phi) for a
     * specific module
     * @tparam visualize_calibration if true visualize calibration process
     * (default: false)
     * @param module_index index of module to be calibrated
     * @param frames_with_base_peak_overlap frames with base peak overlap module
     * to be calibrated
     * @param frames_with_base_peak_overlap_prev_module frames with base peak
     * overlap prev module
     */
    template <bool visualize_calibration = false>
    void calibrate_offset(
        const size_t module_index,
        const std::vector<MythenFrame> &frames_with_base_peak_overlap,
        const std::vector<MythenFrame>
            &frames_with_base_peak_overlap_prev_module);

    /**
     * @brief calibrates the BC (best computing) parameters for one module
     * @tparam visualize_calibration if true visualize calibration process
     * (default: false)
     * @param file_list vector of file_names of acquisition files
     * @param base_peak_angle angle of the selcted base peak given in degrees
     * @param module_index index of module to be calibrated
     */
    template <bool visualize_calibration = false>
    void calibrate(const std::vector<std::string> &file_list_,
                   const double base_peak_angle_, const size_t module_index);

    /**
     * @brief calculates the average angle of the base peak center over several
     * acquisitions for a given module index
     * @param module_index index of module to calculate base peak center for
     * @param frames_with_base_peak_overlap frames with base peak overlap for
     * the given module index
     * @return average angle of base peak center [degrees]
     */
    double center_base_peak(
        const size_t module_index,
        const std::vector<MythenFrame> &frames_with_base_peak_overlap);

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
     * @brief set the scale factor to scale everything to a reasonable scale
     * e.g. can be incident_intensity of first acquisition
     */
    void set_scale_factor(const double scale_factor_);

    /**
     * @brief get the currently configured scale factor (if any)
     */
    double get_scale_factor() const;

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

    /**
     * @brief Performs angular conversion e.g. calculates from raw photon counts
     * the resulting diffraction pattern
     * @param file_list vector of file_names of acquisition files
     * @return flatfield corrected and variance scaled photon counts
     * redistributed to fixed angle width bins given in the range [min_angle,
     * max_angle]
     */
    NDArray<double, 1> convert(const std::vector<std::string> &file_list);

    /**
     * @brief Performs angular conversion e.g. calculates from raw photon counts
     * the resulting diffraction pattern
     * @param file_list vector of file_names of acquisition files
     * @param module_index index of module
     * @return flatfield corrected and variance scaled photon counts
     * redistributed to fixed angle width bins given in the range [min_angle,
     * max_angle]
     */
    NDArray<double, 1> convert(const std::vector<std::string> &file_list,
                               const size_t module_index);

    /** @brief calculates diffraction angle from DG module parameters (used in
     * Beer's Law)
     * @param detector_angle detector position [degrees]
     * @param strip_index local strip index of module e.g. 0-1279
     * @param distance_to_strip distance to strip (if 0.0 calculates diffraction
     * angle at center of strip) [given in strips]
     * @return diffraction angle [degrees]
     */
    double diffraction_angle_from_DG_parameters(
        const size_t module_index, const double detector_angle,
        size_t strip_index, const double distance_to_strip = 0) const;

    /** @brief calculates diffraction angle from BC (best computing) module
     * parameters (used in Beer's Law)
     * @param detector_angle detector position [degrees]
     * @param strip_index local strip index of module e.g. 0-1279
     * @param distance_to_strip distance to strip (if 0.0 calculates diffraction
     * angle at center of strip) [given in strips]
     * @return diffraction angle [degrees]
     */
    double diffraction_angle_from_BC_parameters(
        const size_t module_index, const double detector_angle,
        size_t strip_index, const double distance_to_strip = 0) const;

    /** @brief calculates diffraction angle from EE module parameters (used in
     * Beer's Law)
     * @param detector_angle detector position [degrees]
     * @param strip_index local strip index of module
     * @param distance_to_strip distance to strip (if 0.0 calculates diffraction
     * angle at center of strip) [given in strips]
     * @return diffraction angle [degrees]
     */
    double diffraction_angle_from_EE_parameters(
        const double module_center_distance, const double normal_distance,
        const double angle, const double detector_angle, size_t strip_index,
        const double distance_to_strip = 0) const;

    /** @brief writes DG parameters to file
     * @param filename full path tooutput file
     */
    static void
    write_DG_parameters_to_file(const std::filesystem::path &filename,
                                const DGParameters &parameters);

    /**
     * @brief calculate the rate corrected photon counts taking into account the
     * dead time
     * @param photon_counts photon counts
     * @param photon_count_error error of photon counts (generally variance
     * error is used)
     * @param exposure_time exposure time of acquisition [s]
     * @return pair {corrected photon counts, propagated_error}
     */
    std::pair<double, double> rate_correction(const double photon_counts,
                                              const double photon_count_error,
                                              const double exposure_time) const;

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
     * @brief calculates the elastic correction factor for given detector angle
     * based on torsional compliance correction model
     * @param detector_angle detector axis motor position [degrees]
     */
    double elastic_correction(const double detector_angle) const;

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
        NDView<double, 1> sum_statistical_weights);

    /**
     * @brief redistributes photon counts to fixed angle width bins around base
     * peak for all acquisitions in file_list and calculates the similarity
     * between the found base peaks regions
     * @param frames_with_base_peak_overlap vector of frames from acquisitions
     * for which the base peak ROI is covered by the respective module region
     * @brief plot wrapper to visualize calibration process
     * (default nullptr)
     * @return similarity of peaks
     */
    double calculate_similarity_of_peaks_between_acquisitions(
        const size_t module_index,
        const std::vector<MythenFrame> &frames_with_base_peak_overlap,
        std::shared_ptr<PlotCalibrationProcess> plot = nullptr);

    double calculate_similarity_of_peaks_between_modules(
        const size_t module_index,
        const std::vector<MythenFrame> &frames_with_base_peak_overlap_module,
        const std::vector<MythenFrame>
            &frames_with_base_peak_overlap_prev_module,
        std::shared_ptr<PlotCalibrationProcess> plot = nullptr);

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
    double chi_similarity_criterion(const NDView<double, 1> S0,
                                    const NDView<double, 1> S1,
                                    const NDView<double, 1> S2,
                                    const size_t num_runs) const;

    /**
     * @brief optimizes module center distance L and angle center module normal
     * delta of BC parameters of given module based on chi similarity criterion
     * @param frames_with_base_peak_overlap vector of frames from acquisitions
     * for which the base peak ROI is covered by the respective module region
     * @param plot plot wrapper to visualize calibration process
     * (default nullptr)
     * @param delta_parameter1 parameter step for parameter module center
     * distance L (default '0.5')
     * @param delta_parameter2 parameter step for parameter angle between center
     * of module and module normal (delta) (default '1.0e-6')
     */
    void optimize_coupled_parameters(
        const size_t module_index,
        const std::vector<MythenFrame> &frames_with_base_peak_overlap,
        std::shared_ptr<PlotCalibrationProcess> plot = nullptr,
        const double delta_parameter1 = 0.5,
        const double delta_parameter2 = 1.0e-6);

    /**
     * @brief optimizes angle between center of module and beam (\psi) of BC
     * parameters of given module based on chi similarity criterion and center
     * of base peak
     * @param frames_with_base_peak_overlap vector of frames with base peak
     * overlap for module
     * @param plot plot wrapper to visualize calibration process
     * (default nullptr)
     * @param delta_parameter parameter step for angle between center of module
     * and beam (default '0.005')
     */
    void optimize_offset_parameter(
        const size_t module_index,
        const std::vector<MythenFrame> &frames_with_base_peak_overlap,
        const std::vector<MythenFrame>
            &frames_with_base_peak_overlap_prev_module,
        std::shared_ptr<PlotCalibrationProcess> plot = nullptr,
        double delta_parameter = 0.005);

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

    std::pair<double, double> transverse_width_correction(
        const double photon_counts, const double photon_counts_error,
        const size_t module_index, const size_t strip_index) const;

    /** @brief calculate the corrected photon counts (flatfield, rate,
     * incident intensity)
     * @param photon_counts photon counts
     * @param global_strip_index strip index of photon counts
     * @param I0 incident intensity
     * @param exposure_time exposure time of acquisition [s]
     * @return pair {corrected photon counts, propagated_error}
     */
    std::pair<double, double>
    photon_count_correction(double photon_counts,
                            const size_t global_strip_index, const uint64_t I0,
                            const double exposure_time) const;

    /**
     * @brief appends given parameters to file
     */
    static void append_to_file(std::ofstream &file, const size_t module_index,
                               const double center, const double conversion,
                               const double offset);

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
    DGParameters DGparameters{};

    BCParameters BCparameters{};

    /// @brief Array of size strips/channels in detector, stores 'TRUE' if strip
    /// is a bad channel otherwise 'FALSE'
    NDArray<bool, 1> bad_channels{};

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

    /**
     * center of base peak to use for calibration [degrees]
     */
    double base_peak_angle{};

    /**
     * @brief scale factor to scale everything to a reasonable scale
     * e.g. can be incident_intensity of first acquisition
     * @default 1.0
     */
    double m_scale_factor{1.0};

    /// @brief  @brief minimum angle for angular conversion [degrees]
    double m_min_angle = -180.0;

    /// @brief  maximum angle for angular conversion [degrees]
    double m_max_angle = 180.0;

    /**
     * list of acquisition files used for calibration
     */
    std::vector<std::string> file_list{};

    /**
     * file reader to read mythen acquisition files
     */
    std::shared_ptr<MythenFileReader> mythen_file_reader{};
};

template <bool base_peak_ROI_only>
void AngleCalibration::redistribute_photon_counts_to_fixed_angle_width_bins(
    const size_t module_index, const MythenFrame &frame,
    NDView<double, 1> fixed_angle_width_bins_photon_counts,
    NDView<double, 1> fixed_angle_width_bins_photon_counts_variance,
    NDView<double, 1> sum_statistical_weights) {

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
        left_boundary_roi_base_peak = m_min_angle;  // dummy values
        right_boundary_roi_base_peak = m_max_angle; // dummy values
    }

    // iterate over all strips of module
    for (size_t strip_index = 0;
         strip_index < MythenDetectorSpecifications::strips_per_module;
         ++strip_index) {

        size_t global_strip_index =
            module_index * MythenDetectorSpecifications::strips_per_module +
            strip_index;

        if (bad_channels(global_strip_index)) {
            continue; // skip bad channels
        }

        double left_strip_boundary_angle = diffraction_angle_from_BC_parameters(
            module_index, frame.detector_angle, strip_index,
            -0.5); // left strip boundary in angles [degrees] // -0.5

        double right_strip_boundary_angle =
            diffraction_angle_from_BC_parameters(
                module_index, frame.detector_angle, strip_index,
                +0.5); // right strip boundary in angles [degrees] // +0.5

        if (left_strip_boundary_angle < m_min_angle ||
            right_strip_boundary_angle > m_max_angle) {
            LOG(TLogLevel::logDEBUG1) << fmt::format(
                "strip {} of module {} covers an angle of [{}, {}] and is "
                "outside "
                "of the given angular range of [{}, "
                "{}] and will be skipped",
                strip_index, module_index, left_strip_boundary_angle,
                right_strip_boundary_angle, m_min_angle, m_max_angle);
            continue;
        }

        // skip strip if not in ROI
        if constexpr (base_peak_ROI_only) {
            if (left_strip_boundary_angle > right_boundary_roi_base_peak ||
                right_strip_boundary_angle < left_boundary_roi_base_peak) {
                continue;
            }
        }

        // get corrected photon counts and propagate error
        auto [corrected_photon_counts, corrected_photon_counts_variance] =
            photon_count_correction(
                frame.photon_counts(global_strip_index), global_strip_index,
                frame.incident_intensity, frame.exposure_time);

        double inverse_corrected_photon_counts_variance =
            corrected_photon_counts_variance <
                    std::numeric_limits<double>::epsilon()
                ? 0.0
                : 1.0 / corrected_photon_counts_variance;

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

        if (photon_counts_per_bin == 0.0 ||
            inverse_photon_counts_variance_per_bin == 0.0) {

            bad_channels(global_strip_index) =
                true; // mark channel as bad if no
            continue;
        }

        LOG(TLogLevel::logDEBUG1)
            << fmt::format("corrected photon count for strip {} of module {}",
                           strip_index, module_index);

        ssize_t left_bin_index_covered_by_strip{};

        ssize_t right_bin_index_covered_by_strip{};

        if constexpr (base_peak_ROI_only) {
            left_bin_index_covered_by_strip =
                static_cast<ssize_t>((std::max(left_boundary_roi_base_peak,
                                               left_strip_boundary_angle) /
                                      histogram_bin_width));
            right_bin_index_covered_by_strip =
                static_cast<ssize_t>((std::min(right_boundary_roi_base_peak,
                                               right_strip_boundary_angle) /
                                      histogram_bin_width));
        } else {

            left_bin_index_covered_by_strip = static_cast<ssize_t>(
                (left_strip_boundary_angle / histogram_bin_width));
            right_bin_index_covered_by_strip = static_cast<ssize_t>(
                (right_strip_boundary_angle / histogram_bin_width));

            /*
            // Antonio stuff
            left_bin_index_covered_by_strip = std::floor(
                (left_strip_boundary_angle / histogram_bin_width) + 0.5);
            right_bin_index_covered_by_strip = std::ceil(
                (right_strip_boundary_angle / histogram_bin_width) - 0.5);
            */
        }

        size_t proper_bin_index =
            0; // the computed bin indices dont start at zero but are relative
               // to detector angle and strip index

        for (ssize_t bin_index = left_bin_index_covered_by_strip;
             bin_index <= right_bin_index_covered_by_strip; ++bin_index) {

            // some strips dont cover an entire bin or extend over multiple bins
            /*
            double bin_coverage =
                std::abs(std::min(right_strip_boundary_angle,
                                  (bin_index + 0.5) * histogram_bin_width) -
                         std::max(left_strip_boundary_angle,
                                  (bin_index - 0.5) * histogram_bin_width));
            */

            double bin_coverage =
                std::abs(std::min(right_strip_boundary_angle,
                                  (bin_index + 1) * histogram_bin_width) -
                         std::max(left_strip_boundary_angle,
                                  (bin_index)*histogram_bin_width));

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
                        m_min_angle /
                        histogram_bin_width); // bin index starts at zero
            }

            double statistical_weight =
                bin_coverage_factor * inverse_photon_counts_variance_per_bin;

            fixed_angle_width_bins_photon_counts(proper_bin_index) +=
                statistical_weight * photon_counts_per_bin;

            // TODO: need to fix !!!
            fixed_angle_width_bins_photon_counts_variance(proper_bin_index) +=
                inverse_photon_counts_variance_per_bin *
                std::pow(bin_coverage_factor,
                         2); // variance * statistical_weight^2 = variance⁻2 *
                             // bin_coverage_factor^2

            sum_statistical_weights(proper_bin_index) += statistical_weight;
        }
    }
}

template <bool visualize_calibration>
void AngleCalibration::calibrate_offset(
    const std::vector<double> &detector_angles) {

    for (size_t module_index = 1; module_index < mythen_detector->max_modules;
         ++module_index) {

        LOG(TLogLevel::logINFO)
            << fmt::format("starting calibration for module {} ", module_index);

        const auto t0 = std::chrono::steady_clock::now();

        // skip if module is not connected
        if (module_is_disconnected(module_index) or
            module_is_disconnected(module_index - 1)) {
            LOG(TLogLevel::logINFO)
                << fmt::format("module {} or module {} is disconnected",
                               module_index, module_index - 1);
            continue;
        }

        std::vector<MythenFrame> frames_with_base_peak_overlap;
        frames_with_base_peak_overlap.reserve(file_list.size());

        double potential_change_in_diffraction_angle = base_peak_roi_width;
        for (size_t i = 0; i < detector_angles.size(); ++i) {
            if (base_peak_is_in_module(
                    module_index, detector_angles[i] -
                                      potential_change_in_diffraction_angle) ||
                base_peak_is_in_module(
                    module_index, detector_angles[i] +
                                      potential_change_in_diffraction_angle) ||
                base_peak_is_in_module(module_index, detector_angles[i])) {

                frames_with_base_peak_overlap.push_back(
                    mythen_file_reader->read_frame(file_list[i]));
            }
        }

        std::vector<MythenFrame> frames_with_base_peak_overlap_prev_module;
        frames_with_base_peak_overlap_prev_module.reserve(file_list.size());

        for (size_t i = 0; i < detector_angles.size(); ++i) {
            if (base_peak_is_in_module(
                    module_index - 1,
                    detector_angles[i] -
                        potential_change_in_diffraction_angle) ||
                base_peak_is_in_module(
                    module_index - 1,
                    detector_angles[i] +
                        potential_change_in_diffraction_angle) ||
                base_peak_is_in_module(module_index - 1, detector_angles[i])) {
                frames_with_base_peak_overlap_prev_module.push_back(
                    mythen_file_reader->read_frame(file_list[i]));
            }
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
                module_index, frames_with_base_peak_overlap,
                frames_with_base_peak_overlap_prev_module, plot, 1.0e-3);
        } catch (const NoBasePeakOverLapError &e) {
            LOG(angcal::TLogLevel::logINFO)
                << e.what()
                << fmt::format(" Skipping module {}.", module_index);
            continue;
        }

        const auto t1 = std::chrono::steady_clock::now();
        const auto elapsed_time =
            std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
                .count();
        LOG(TLogLevel::logINFO)
            << fmt::format("calibration for module {} took {} milliseconds",
                           module_index, elapsed_time);
    }
}

template <bool visualize_calibration>
void AngleCalibration::calibrate_offset(
    const size_t module_index,
    const std::vector<MythenFrame> &frames_with_base_peak_overlap,
    const std::vector<MythenFrame> &frames_with_base_peak_overlap_prev_module) {

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
        optimize_offset_parameter(module_index, frames_with_base_peak_overlap,
                                  frames_with_base_peak_overlap_prev_module,
                                  plot, 1.0e-3);
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
    const std::vector<double> &detector_angles) {

    for (size_t module_index = 0; module_index < mythen_detector->max_modules;
         ++module_index) {

        LOG(TLogLevel::logINFO)
            << fmt::format("starting calibration for module {} ", module_index);

        const auto t0 = std::chrono::steady_clock::now();

        std::vector<MythenFrame> frames_with_base_peak_overlap;
        frames_with_base_peak_overlap.reserve(file_list.size());

        double potential_change_in_diffraction_angle = base_peak_roi_width;
        for (size_t i = 0; i < detector_angles.size(); ++i) {
            if (base_peak_is_in_module(
                    module_index, detector_angles[i] -
                                      potential_change_in_diffraction_angle) ||
                base_peak_is_in_module(
                    module_index, detector_angles[i] +
                                      potential_change_in_diffraction_angle) ||
                base_peak_is_in_module(module_index, detector_angles[i])) {

                frames_with_base_peak_overlap.push_back(
                    mythen_file_reader->read_frame(file_list[i]));
            }
        }

        // skip if module is not connected
        if (!module_is_disconnected(module_index)) {

            std::shared_ptr<PlotCalibrationProcess> plot = nullptr;
            if constexpr (visualize_calibration) {
                std::string plot_title =
                    fmt::format("Base Peaks for module {} ", module_index);
                plot =
                    std::make_shared<PlotCalibrationProcess>(this, plot_title);

                plot->initializematplotlib();
            }
            try {
                optimize_coupled_parameters(module_index,
                                            frames_with_base_peak_overlap, plot,
                                            0.5, 1.0e-6);
            } catch (const NoBasePeakOverLapError &e) {
                LOG(angcal::TLogLevel::logINFO)
                    << e.what()
                    << fmt::format(" Skipping module {}.", module_index);
                continue;
            }

        } else {
            LOG(TLogLevel::logINFO)
                << fmt::format("module {} is disconnected", module_index);
        }

        const auto t1 = std::chrono::steady_clock::now();
        const auto elapsed_time =
            std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
                .count();
        LOG(TLogLevel::logINFO)
            << fmt::format("calibration for module {} took {} milliseconds",
                           module_index, elapsed_time);
    }
}

template <bool visualize_calibration>
void AngleCalibration::calibrate_coupled_parameters(
    const size_t module_index,
    const std::vector<MythenFrame> &frames_with_base_peak_overlap) {

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
        optimize_coupled_parameters(module_index, frames_with_base_peak_overlap,
                                    plot, 0.5, 1.0e-6);
    } catch (const NoBasePeakOverLapError &e) {
        LOG(angcal::TLogLevel::logINFO)
            << e.what() << fmt::format(" Skipping module {}.", module_index);
        return;
    }

    LOG(TLogLevel::logDEBUG)
        << fmt::format("calibration for module {} finished", module_index);
}

template <bool visualize_calibration>
void AngleCalibration::calibrate(
    const std::vector<std::string> &file_list_, const double base_peak_angle_,
    std::optional<std::filesystem::path> output_file) {

    file_list = file_list_;

    base_peak_angle = base_peak_angle_;

    std::vector<double> detector_angles;
    detector_angles.reserve(file_list.size());
    // iterate over all files to get detector angles
    for (const auto &file : file_list) {
        double detector_angle =
            mythen_file_reader->read_detector_angle(file); // TODO: read angle
        detector_angles.push_back(detector_angle);
    }

    std::ofstream file;
    if (output_file.has_value()) {
        file.open(output_file.value(), std::ios::out | std::ios::app);
        if (file.tellp() == 0) {
            file << "module,center,conversion,offset\n"; // only write
                                                         // header if
                                                         // file is empty
        }
    }

    calibrate_coupled_parameters<visualize_calibration>(detector_angles);

    if (output_file.has_value()) {
        BCparameters.convert_to_DGParameters(DGparameters);
        write_DG_parameters_to_file(output_file.value(), DGparameters);
    }

    calibrate_offset<visualize_calibration>(detector_angles);

    if (output_file.has_value()) {
        BCparameters.convert_to_DGParameters(DGparameters);
        write_DG_parameters_to_file(output_file.value(), DGparameters);
    }

    if (output_file.has_value()) {
        file.close();
    }

    PlotCalibrationProcess::kill_python_interpreter();
}

template <bool visualize_calibration>
void AngleCalibration::calibrate(const std::vector<std::string> &file_list_,
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
    file_list = file_list_;

    std::vector<MythenFrame> frames_with_base_peak_overlap;
    frames_with_base_peak_overlap.reserve(file_list.size());

    // iterate over all files to get detector angles
    double potential_change_in_diffraction_angle = base_peak_roi_width;
    for (const auto &file : file_list) {
        double detector_angle =
            mythen_file_reader->read_detector_angle(file); // TODO: read angle
        if (base_peak_is_in_module(module_index,
                                   detector_angle -
                                       potential_change_in_diffraction_angle) ||
            base_peak_is_in_module(module_index,
                                   detector_angle +
                                       potential_change_in_diffraction_angle) ||
            base_peak_is_in_module(module_index, detector_angle)) {

            LOG(TLogLevel::logDEBUG1)
                << fmt::format("file {} with detector angle {} has base peak "
                               "overlap with module {}",
                               file, detector_angle, module_index);
            frames_with_base_peak_overlap.push_back(
                mythen_file_reader->read_frame(file));
        }
    }

    const auto t0 = std::chrono::steady_clock::now();

    calibrate_coupled_parameters<visualize_calibration>(
        module_index, frames_with_base_peak_overlap);

    if (module_index > 0) {
        // iterate over all files to get detector angles
        std::vector<MythenFrame> frames_with_base_peak_overlap_prev_module;
        frames_with_base_peak_overlap_prev_module.reserve(file_list.size());

        for (const auto &file : file_list) {
            double detector_angle = mythen_file_reader->read_detector_angle(
                file); // TODO: read angle
            if (base_peak_is_in_module(
                    module_index - 1,
                    detector_angle - potential_change_in_diffraction_angle) ||
                base_peak_is_in_module(
                    module_index - 1,
                    detector_angle + potential_change_in_diffraction_angle) ||
                base_peak_is_in_module(module_index - 1, detector_angle)) {
                frames_with_base_peak_overlap_prev_module.push_back(
                    mythen_file_reader->read_frame(file));
            }
        }

        calibrate_coupled_parameters<visualize_calibration>(
            module_index - 1, frames_with_base_peak_overlap_prev_module);

        calibrate_offset<visualize_calibration>(
            module_index, frames_with_base_peak_overlap,
            frames_with_base_peak_overlap_prev_module);
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
