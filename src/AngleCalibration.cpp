#include "AngleCalibration.hpp"
#include "Optimization.hpp"
#include "helpers/ErrorPropagation.hpp"
// #include "helpers/LogFiles.hpp"
#include "logger.hpp"

#include <csignal>
#include <cstdlib>
#include <filesystem>

namespace angcal {

/*
inline void handle_sigint(int) {
    LOG(TLogLevel::logINFO) << "\nSIGINT caught — shutting down...\n";
    std::exit(EXIT_SUCCESS);
}
struct SignalHandler {
    SignalHandler() { std::signal(SIGINT, handle_sigint); }
};
*/

AngleCalibration::AngleCalibration(
    std::shared_ptr<MythenDetectorSpecifications> mythen_detector_,
    std::shared_ptr<FlatField> flat_field_,
    std::shared_ptr<MythenFileReader> mythen_file_reader_)
    : mythen_detector(mythen_detector_), flat_field(flat_field_),
      mythen_file_reader(mythen_file_reader_) {

    DGparameters = DGParameters(mythen_detector->max_modules);
    BCparameters = BCParameters(mythen_detector->max_modules);

    bad_channels = NDArray<bool, 1>(std::array<ssize_t, 1>{static_cast<ssize_t>(
                                        mythen_detector->num_strips())},
                                    false);
}

void AngleCalibration::set_histogram_bin_width(double bin_width) {
    histogram_bin_width = bin_width;
}

double AngleCalibration::get_histogram_bin_width() const {
    return histogram_bin_width;
}

ssize_t AngleCalibration::get_base_peak_ROI_num_bins() const {

    return static_cast<ssize_t>((base_peak_angle + base_peak_roi_width) /
                                histogram_bin_width) -
           static_cast<ssize_t>((base_peak_angle - base_peak_roi_width) /
                                histogram_bin_width) +
           1;
}

void AngleCalibration::set_base_peak_ROI_width(
    const double base_peak_roi_width_) {
    base_peak_roi_width = base_peak_roi_width_;
}

double AngleCalibration::get_base_peak_ROI_width() const {
    return base_peak_roi_width;
}

double AngleCalibration::get_base_peak_angle() const { return base_peak_angle; }

void AngleCalibration::set_base_peak_angle(const double base_peak_angle_) {
    base_peak_angle = base_peak_angle_;
}

void AngleCalibration::set_calibration_files(
    const std::vector<std::string> &file_list_) {
    file_list = file_list_;
}

const std::vector<std::string> &
AngleCalibration::get_calibration_files() const {
    return file_list;
}

std::shared_ptr<MythenDetectorSpecifications>
AngleCalibration::get_detector_specifications() const {
    return mythen_detector;
}

ssize_t AngleCalibration::num_fixed_angle_width_bins() const {
    ssize_t num_fixed_angle_width_bins =
        std::ceil(m_max_angle / histogram_bin_width) -
        std::floor(m_min_angle / histogram_bin_width);
    return num_fixed_angle_width_bins;
}

const DGParameters &AngleCalibration::get_DGparameters() const {
    return DGparameters;
}

const BCParameters &AngleCalibration::get_BCparameters() const {
    return BCparameters;
}

// TODO: command is connected - maybe use that
bool AngleCalibration::module_is_disconnected(const size_t module_index) const {

    // all channels in the module are bad
    return std::all_of(bad_channels.begin() +
                           module_index * mythen_detector->strips_per_module,
                       bad_channels.begin() +
                           (module_index + 1) *
                               mythen_detector->strips_per_module,
                       [](const auto &elem) { return elem; });
}

void AngleCalibration::set_scale_factor(const double scale_factor_) {
    m_scale_factor = scale_factor_;
}

double AngleCalibration::get_scale_factor() const { return m_scale_factor; }

void AngleCalibration::set_angular_range(const double min_angle,
                                         const double max_angle) {
    m_min_angle = min_angle;
    m_max_angle = max_angle;
}

std::pair<double, double> AngleCalibration::get_angular_range() const {
    return {m_min_angle, m_max_angle};
}

void AngleCalibration::read_initial_calibration_from_file(
    const std::string &filename,
    const std::shared_ptr<SimpleFileInterface> file_reader) {

    file_reader->open(filename);
    file_reader->read_into(DGparameters.parameters.buffer(), 8);

    file_reader->close();

    DGparameters.convert_to_BCParameters(
        BCparameters); // initialize BC parameters
}

void AngleCalibration::read_bad_channels_from_file(
    const std::string &filename,
    std::shared_ptr<SimpleFileInterface> file_reader) {

    file_reader->open(filename);
    file_reader->read_into(reinterpret_cast<std::byte *>(bad_channels.data()));

    // update bad_channels
    std::for_each(
        mythen_detector->unconnected_modules.begin(),
        mythen_detector->unconnected_modules.end(),
        [this](const auto &module_index) {
            for (size_t i = module_index *
                            MythenDetectorSpecifications::strips_per_module;
                 i < (module_index + 1) *
                         MythenDetectorSpecifications::strips_per_module;
                 ++i)
                bad_channels[i] = true;
        });
}

NDView<bool, 1> AngleCalibration::get_bad_channels() const {
    return bad_channels.view();
}

void AngleCalibration::set_bad_channels(const NDArray<bool, 1> &bad_channels_) {
    bad_channels = bad_channels_;

    // update bad_channels
    std::for_each(
        mythen_detector->unconnected_modules.begin(),
        mythen_detector->unconnected_modules.end(),
        [this](const auto &module_index) {
            for (size_t i = module_index *
                            MythenDetectorSpecifications::strips_per_module;
                 i < (module_index + 1) *
                         MythenDetectorSpecifications::strips_per_module;
                 ++i)
                this->bad_channels[i] = true;
        });
}

size_t AngleCalibration::global_to_local_strip_index_conversion(
    const size_t global_strip_index) const {
    const size_t module_index =
        global_strip_index / MythenDetectorSpecifications::strips_per_module;
    // local strip index in module
    size_t local_strip_index =
        global_strip_index -
        module_index * MythenDetectorSpecifications::strips_per_module;
    // switch if indexing is in clock-wise direction
    local_strip_index = std::signbit(DGparameters.conversions(module_index))
                            ? MythenDetectorSpecifications::strips_per_module -
                                  1 - local_strip_index
                            : local_strip_index;

    return local_strip_index;
}

double AngleCalibration::sample_displacement_correction(
    const double diffraction_angle) const {

    /*
    double sample_displacement_correction =
        std::atan(std::tan(M_PI / 180.0 * diffraction_angle) *
                  std::sqrt(std::pow(mythen_detector->sample_x_displacement, 2)
    + std::pow(mythen_detector->sample_y_displacement, 2)) /
                  mythen_detector->average_distance_sample_pixel);
    */

    double whatever1 = std::sin(M_PI / 180.0 * diffraction_angle) +
                       mythen_detector->sample_y_displacement /
                           mythen_detector->average_distance_sample_pixel;

    double whatever2 = std::cos(M_PI / 180.0 * diffraction_angle) -
                       mythen_detector->sample_x_displacement /
                           mythen_detector->average_distance_sample_pixel;

    if (std::abs(whatever2) < std::numeric_limits<double>::epsilon()) {
        LOG(TLogLevel::logWARNING) << fmt::format(
            "Sample displacement correction is unstable for diffraction "
            "angle {} degrees due to small denominator in arctan.",
            diffraction_angle);
        return diffraction_angle;
    }

    double sample_displacement_correction =
        180.0 / M_PI * std::atan(whatever1 / whatever2);

    // preserve original angle information in range [-180, 180]
    const bool diffraction_angle_in_second_quadrant =
        diffraction_angle / 90.0 > 1.0;
    const bool diffraction_angle_in_third_quadrant =
        diffraction_angle / 90.0 < -1.0;

    if (diffraction_angle_in_second_quadrant) {
        return 180.0 + sample_displacement_correction;
    } else if (diffraction_angle_in_third_quadrant) {
        return -180.0 + sample_displacement_correction;
    } else {
        return sample_displacement_correction;
    }
}

double AngleCalibration::elastic_correction(const double detector_angle) const {
    return detector_angle +
           mythen_detector->elastic_correction_factor *
               (-std::sin(M_PI / 180.0 *
                          (detector_angle -
                           mythen_detector->detector_vertical_axis_offset)));
}

double AngleCalibration::diffraction_angle_from_DG_parameters(
    const size_t module_index, const double detector_angle, size_t strip_index,
    const double distance_to_strip) const {

    double offset = DGparameters.offsets(module_index);
    double center = DGparameters.centers(module_index);
    double conversion = DGparameters.conversions(module_index);

    strip_index =
        std::signbit(conversion)
            ? MythenDetectorSpecifications::strips_per_module - strip_index - 1
            : strip_index; // TODO: are the values sored in reserve?

    double diffraction_angle =
        offset +
        180.0 / M_PI *
            (center * std::abs(conversion) -
             std::atan((center - (strip_index + distance_to_strip)) *
                       std::abs(conversion))) +
        elastic_correction(detector_angle) + mythen_detector->offset +
        mythen_detector->sample_detector_offset;

    return sample_displacement_correction(diffraction_angle);
}

double AngleCalibration::diffraction_angle_from_BC_parameters(
    const size_t module_index, const double detector_angle, size_t strip_index,
    const double distance_to_strip) const {

    const double angle_module_center_normal =
        BCparameters.angle_center_module_normal(module_index);
    const double distance_center_sample =
        BCparameters.module_center_sample_distances(module_index);
    const double angle_module_center_beam =
        BCparameters.angle_center_beam(module_index);

    strip_index =
        std::signbit(distance_center_sample)
            ? MythenDetectorSpecifications::strips_per_module - strip_index - 1
            : strip_index; // TODO: are the values sored in reserve?

    double diffraction_angle =
        angle_module_center_beam + angle_module_center_normal -
        180.0 / M_PI *
            std::atan((std::abs(distance_center_sample) *
                           std::sin(M_PI / 180.0 * angle_module_center_normal) +
                       (MythenDetectorSpecifications::strips_per_module * 0.5 -
                        strip_index - distance_to_strip) *
                           MythenDetectorSpecifications::pitch) /
                      (std::abs(distance_center_sample) *
                       std::cos(M_PI / 180.0 * angle_module_center_normal))) +
        elastic_correction(detector_angle) + mythen_detector->offset +
        mythen_detector->sample_detector_offset;

    return sample_displacement_correction(diffraction_angle);
}

double AngleCalibration::diffraction_angle_from_EE_parameters(
    const double module_center_distance, const double normal_distance,
    const double angle, const double detector_angle, size_t strip_index,
    const double distance_to_strip) const {

    strip_index =
        std::signbit(normal_distance)
            ? MythenDetectorSpecifications::strips_per_module - strip_index - 1
            : strip_index;

    double diffraction_angle =
        angle -
        180.0 / M_PI *
            std::atan((module_center_distance -
                       MythenDetectorSpecifications::pitch *
                           (strip_index + distance_to_strip)) /
                      std::abs(normal_distance)) +
        elastic_correction(detector_angle) + mythen_detector->offset +
        mythen_detector->sample_detector_offset;

    return sample_displacement_correction(diffraction_angle);
}

double AngleCalibration::angular_strip_width_from_DG_parameters(
    const size_t module_index, const size_t local_strip_index) const {

    return std::abs(diffraction_angle_from_DG_parameters(
                        module_index, 0.0, local_strip_index, 0.0) -
                    diffraction_angle_from_DG_parameters(
                        module_index, 0.0, local_strip_index, 1.0));
}

double AngleCalibration::angular_strip_width_from_BC_parameters(
    const size_t module_index, const size_t local_strip_index) const {

    return std::abs(diffraction_angle_from_BC_parameters(
                        module_index, 0.0, local_strip_index, 0.0) -
                    diffraction_angle_from_BC_parameters(
                        module_index, 0.0, local_strip_index, 1.0));
}

double AngleCalibration::angular_strip_width_from_EE_parameters(
    const double module_center_distance, const double normal_distance,
    const double angle, const size_t local_strip_index) const {

    return std::abs(diffraction_angle_from_EE_parameters(
                        module_center_distance, normal_distance, angle, 0.0,
                        local_strip_index, 0.0) -
                    diffraction_angle_from_EE_parameters(
                        module_center_distance, normal_distance, angle, 0.0,
                        local_strip_index, 1.0));
}

double AngleCalibration::chi_similarity_criterion(
    const NDView<double, 2> weighted_sums,
    const NDView<size_t, 1> num_contributing_photons) const {

    double similarity_criterion = 0;

    /*
    if (num_runs == 0) {
        throw std::runtime_error(
            LOCATION +
            "Number of runs is zero - cannot calculate similarity criterion!");
    }
    */
    for (ssize_t bin_index = 0; bin_index < weighted_sums.shape(0);
         ++bin_index) {

        double weighted_average =
            weighted_sums(bin_index, 0) < std::numeric_limits<double>::epsilon()
                ? 0.0
                : 1. / weighted_sums(bin_index, 0);

        double goodness_of_fit =
            num_contributing_photons(bin_index) == 1
                ? 0.0
                : (weighted_sums(bin_index, 2) -
                   weighted_sums(bin_index, 1) * weighted_sums(bin_index, 1) *
                       weighted_average) /
                      (num_contributing_photons(bin_index) -
                       1); // calculates chi value for optimal
                           // parameter a over all runs chi_bin

        similarity_criterion +=
            goodness_of_fit * weighted_average; // scale by variance
    }

    return similarity_criterion; // TODO: do we have to didvide by num bins? -
                                 // probably?
}

// TODO maybe inline this
bool AngleCalibration::base_peak_is_in_module(
    const size_t module_index, const double detector_angle) const {

    double left_module_boundary_angle =
        diffraction_angle_from_BC_parameters(module_index, detector_angle, 0);
    double right_module_boundary_angle = diffraction_angle_from_BC_parameters(
        module_index, detector_angle, mythen_detector->strips_per_module - 1);

    LOG(TLogLevel::logDEBUG2) << fmt::format(
        "module_boundaries_in_angle for module {} [{}, {}]\n", module_index,
        left_module_boundary_angle, right_module_boundary_angle);

    // TODO check bounds with base_peak_hwid - TODO: I changed it such that it
    // needs to be fully in module - otherwise similarity criterion gets
    // distorted
    return (base_peak_angle + base_peak_roi_width <
                right_module_boundary_angle &&
            left_module_boundary_angle < base_peak_angle - base_peak_roi_width);
}

std::pair<double, double>
AngleCalibration::rate_correction(const double photon_count,
                                  const double photon_count_error,
                                  const double exposure_time) const {

    if (exposure_time == 0) {
        throw std::runtime_error(LOCATION +
                                 "Exposure time is zero - cannot perform rate "
                                 "correction!");
    }

    const double maximum_count_rate =
        std::exp(-1); // theoretical maximum count rate for
                      // dead_time*measured_photon_counts_per_second

    double photon_counts_per_second = photon_count / exposure_time;

    photon_counts_per_second =
        std::min(maximum_count_rate,
                 photon_counts_per_second *
                     mythen_detector->dead_time); // multiply with dead time -
                                                  // for numerical algorithm

    double error_photon_counts_per_second =
        photon_count_error *
        std::pow(mythen_detector->dead_time / exposure_time, 2);

    // -actual_count_rate*dead_time = W -> -dead_time*measured_count_rate =
    // We^{W} -> W: Lambert W function

    // numerical calculation of inverse of Lambart W - function
    // actually calculate dead_time*count_rate =
    // dead_time*actual_count_rate*exp(-dead_time*actual_count_rate)
    double W_prev_iter = photon_counts_per_second; // initial guess
    double W_next_iter{};
    double propagated_error = error_photon_counts_per_second;
    bool method_converged = false;
    while (!method_converged) {
        W_next_iter = photon_counts_per_second * std::exp(W_prev_iter);
        method_converged = std::abs(1 - W_next_iter / W_prev_iter) <
                           10 * std::numeric_limits<double>::epsilon();
        W_prev_iter = W_next_iter;
        propagated_error =
            error_propagation(std::exp(W_prev_iter),
                              photon_counts_per_second * std::exp(W_prev_iter),
                              error_photon_counts_per_second, propagated_error);
    }

    double actual_count_rate = W_prev_iter; // actual count rate * dead_time

    double rate_correction_factor =
        actual_count_rate / photon_counts_per_second;

    /*
    double error_rate_correction_factor = error_propagation(
        1. / photon_counts_per_second,
        -actual_count_rate / std::pow(photon_counts_per_second, 2),
        propagated_error, error_photon_counts_per_second);
    */

    double rate_corrected_photon_counts = photon_count * rate_correction_factor;

    // double rate_corrected_photon_counts_error =
    // error_propagation(rate_correction_factor, photon_count,
    // photon_count_error, error_rate_correction_factor);

    double rate_corrected_photon_counts_error =
        photon_count_error *
        std::pow(rate_correction_factor, 2); // TODO: In Antonios code like this

    return std::pair<double, double>{rate_corrected_photon_counts,
                                     rate_corrected_photon_counts_error};
}

std::pair<double, double> AngleCalibration::incident_intensity_correction(
    const double photon_counts, const double photon_count_error,
    const uint64_t incident_intensity) const {

    if (incident_intensity == 0) {
        throw std::runtime_error(LOCATION +
                                 "Incident intensity is zero - cannot perform "
                                 "incident intensity correction");
    }

    double I0_correction_factor = get_scale_factor() / (incident_intensity);

    double I0_corrected_photon_counts = photon_counts * I0_correction_factor;
    double I0_corrected_photon_counts_error =
        photon_count_error * std::pow(I0_correction_factor, 2);

    return std::pair<double, double>{I0_corrected_photon_counts,
                                     I0_corrected_photon_counts_error};
}

std::pair<double, double>
AngleCalibration::flatfield_correction(const double photon_counts,
                                       const double photon_counts_error,
                                       const size_t global_strip_index) const {

    if ((flat_field->get_normalized_flatfield_view()(global_strip_index, 0)) <=
        0.1) { // std::numeric_limits<double>::epsilon() in Antonios code 0.1

        LOG(TLogLevel::logDEBUG)
            << fmt::format("Flatfield value for strip {}", global_strip_index);
        return std::pair(0.0, 0.0);
    }

    // flatfield normalization
    const double flatfield_normalized_photon_counts =
        photon_counts /
        flat_field->get_normalized_flatfield_view()(global_strip_index, 0);

    const double normalized_flatfield_variance = std::pow(
        flat_field->get_normalized_flatfield_view()(global_strip_index, 1), 2);

    const double flatfield_normalized_photon_counts_error = error_propagation(
        1. / flat_field->get_normalized_flatfield_view()(global_strip_index, 0),
        -1.0 * photon_counts /
            std::pow(flat_field->get_normalized_flatfield_view()(
                         global_strip_index, 0),
                     2),
        photon_counts_error,
        normalized_flatfield_variance); // Variance // TODO: change should be
                                        // squared per default always talk about
                                        // variance error

    return std::pair(flatfield_normalized_photon_counts,
                     flatfield_normalized_photon_counts_error);
}

std::pair<double, double> AngleCalibration::transverse_width_correction(
    const double photon_counts, const double photon_counts_error,
    const size_t module_index, size_t strip_index) const {

    // convert to EE parameters
    const auto [normal_distance, module_center_distance, angle] =
        BCparameters.convert_to_EEParameters(module_index);

    strip_index =
        std::signbit(normal_distance)
            ? MythenDetectorSpecifications::strips_per_module - 1 - strip_index
            : strip_index;

    const double distance_sample_pixel = std::sqrt(
        std::pow(normal_distance, 2) +
        std::pow(module_center_distance - mythen_detector->pitch * strip_index,
                 2));

    const double transverse_width_correction_factor =
        (2 * std::atan(mythen_detector->transverse_width /
                       (2 * mythen_detector->average_distance_sample_pixel))) /
        (2 * std::atan(mythen_detector->transverse_width /
                       (2 * distance_sample_pixel)));

    double transverse_width_normalized_photon_counts =
        photon_counts * transverse_width_correction_factor;
    double transverse_width_normalized_photon_counts_error =
        photon_counts_error * std::pow(transverse_width_correction_factor, 2);

    return std::pair(transverse_width_normalized_photon_counts,
                     transverse_width_normalized_photon_counts_error);
}

std::pair<double, double> AngleCalibration::photon_count_correction(
    double photon_counts, const size_t global_strip_index, const uint64_t I0,
    const double exposure_time) const {

    if (photon_counts < std::numeric_limits<double>::epsilon()) {
        return std::pair(0.0, 0.0); // sanity check
    }

    // mighells statistics
    photon_counts += 1;

    auto [rate_corrected_photon_counts, rate_corrected_photon_counts_error] =
        rate_correction(photon_counts, photon_counts,
                        exposure_time); // same variance as poisson distributed
                                        // - maybe sqaured?

    auto [incident_intensity_corrected_photon_counts,
          incident_intensity_corrected_photon_counts_variance] =
        incident_intensity_correction(rate_corrected_photon_counts,
                                      rate_corrected_photon_counts_error, I0);

    const size_t module_index =
        global_strip_index / MythenDetectorSpecifications::strips_per_module;
    const size_t strip_index =
        global_strip_index % MythenDetectorSpecifications::strips_per_module;

    auto [transverse_width_corrected_photon_counts,
          transverse_width_corrected_photon_counts_variance] =
        transverse_width_correction(
            incident_intensity_corrected_photon_counts,
            incident_intensity_corrected_photon_counts_variance, module_index,
            strip_index);

    auto [flatfield_normalized_photon_counts,
          flatfield_normalized_photon_counts_variance] =
        flatfield_correction(
            transverse_width_corrected_photon_counts,
            transverse_width_corrected_photon_counts_variance,
            global_strip_index); // make sure to include bad chanenls of
                                 // flatfield into bad channel list

    return std::pair<double, double>{
        flatfield_normalized_photon_counts,
        flatfield_normalized_photon_counts_variance};
}

double AngleCalibration::center_base_peak(
    const size_t module_index,
    const std::vector<MythenFrame> &frames_with_base_peak_overlap) {

    if (module_is_disconnected(module_index)) {
        throw std::runtime_error(
            LOCATION + "Module is disconnected - cannot calculate center base "
                       "peak in module");
    }

    ssize_t num_bins_in_ROI = get_base_peak_ROI_num_bins();

    size_t num_runs = 0;

    double average_peak_center{};

    for (const auto &frame : frames_with_base_peak_overlap) {

        double sum_intensities_in_ROI = 0.0;
        double weighted_sum = 0.0;

        // base peak angle is in module
        if (base_peak_is_in_module(module_index, frame.detector_angle)) {

            NDArray<double, 2> weighted_sums(
                std::array<ssize_t, 2>{num_bins_in_ROI, 3},
                -1.0); // // weighted moments of order 0-2 for each bin

            NDArray<size_t, 1> num_contributing_photons(
                std::array<ssize_t, 1>{num_bins_in_ROI},
                0); // number of contributing photons for each bin

            // calculates flatfield normalized photon counts and
            // photon_count_variance for ROI around base_peak and
            // redistributes to fixed angle width bins
            redistribute_photon_counts_to_fixed_angle_width_bins<true>(
                module_index, frame, weighted_sums.view(),
                num_contributing_photons.view());

            for (ssize_t i = 0; i < weighted_sums.shape(0); ++i) {
                if (weighted_sums(i, 0) == -1.0) {
                    continue;
                } else {
                    double
                        photon_counts_redistributed_to_fixed_angle_width_bin =
                            weighted_sums(i, 0) <
                                    std::numeric_limits<double>::epsilon()
                                ? 0.0
                                : weighted_sums(i, 1) /
                                      weighted_sums(i, 0); // y_k
                    sum_intensities_in_ROI +=
                        photon_counts_redistributed_to_fixed_angle_width_bin;
                    weighted_sum +=
                        i *
                        photon_counts_redistributed_to_fixed_angle_width_bin;
                }
            }

            // centerof mass - peak center converted to angles
            average_peak_center +=
                weighted_sum / sum_intensities_in_ROI * histogram_bin_width -
                base_peak_roi_width + base_peak_angle;

            ++num_runs;
        }
    }
    if (num_runs == 0) {
        throw NoBasePeakOverLapError();
    }

    return average_peak_center / num_runs;
}

NDArray<double, 2> AngleCalibration::computeErrorandMean(
    const NDView<double, 2> weighted_sums,
    const NDView<size_t, 1> num_contributing_photons) const {

    const ssize_t new_num_bins = num_fixed_angle_width_bins();

    NDArray<double, 2> fixed_angle_width_bins_photon_counts(
        std::array<ssize_t, 2>{new_num_bins, 2},
        -1.0); // first dimension photon counts, second dimension error

    // divide by statistial weight
    for (ssize_t i = 0; i < weighted_sums.shape(0); ++i) {

        if (weighted_sums(i, 1) > -1.0 && num_contributing_photons(i) > 0) {
            fixed_angle_width_bins_photon_counts(i, 0) =
                weighted_sums(i, 0) < std::numeric_limits<double>::epsilon()
                    ? 0.0
                    : weighted_sums(i, 1) /
                          weighted_sums(i, 0); // weighted mean

            double inverse_weights = 1. / weighted_sums(i, 0);

            double goodness_of_fit =
                num_contributing_photons(i) == 1
                    ? 0.0
                    : (weighted_sums(i, 2) - weighted_sums(i, 1) *
                                                 weighted_sums(i, 1) *
                                                 inverse_weights) /
                          (num_contributing_photons(i) -
                           1); // calculates chi value for optimal

            fixed_angle_width_bins_photon_counts(i, 1) = std::sqrt(
                goodness_of_fit * inverse_weights); // scale by variance
        }
    }

    return fixed_angle_width_bins_photon_counts;
}

NDArray<double, 2>
AngleCalibration::convert(const std::vector<std::string> &file_list_) {

    ssize_t new_num_bins = num_fixed_angle_width_bins();

    LOG(TLogLevel::logDEBUG)
        << fmt::format("num_fixed_angle_width_bins: {}", new_num_bins);

    NDArray<double, 2> weighted_moments(
        std::array<ssize_t, 2>{new_num_bins, 3},
        -1.0); // weighted moments of order 0-2 for each bin
    NDArray<size_t, 1> num_contributing_photons(
        std::array<ssize_t, 1>{new_num_bins},
        0); // number of contributing photons for each bin

    for (const auto &file : file_list_) {
        MythenFrame frame = mythen_file_reader->read_frame(file);

        LOG(TLogLevel::logDEBUG)
            << fmt::format("redistributing photon counts for file: {}", file);

        for (size_t module_index = 0;
             module_index < mythen_detector->max_modules; ++module_index) {

            if (module_is_disconnected(module_index)) {
                continue;
            }

            redistribute_photon_counts_to_fixed_angle_width_bins<false>(
                module_index, frame, weighted_moments.view(),
                num_contributing_photons.view());
        }

        LOG(TLogLevel::logDEBUG) << fmt::format(
            "finished redistributing photon counts for file: {}", file);
    }

    LOG(TLogLevel::logDEBUG)
        << "redistributed photon counts to fixed angle width bins, now "
           "normalizing by statistical weight";

    return computeErrorandMean(weighted_moments.view(),
                               num_contributing_photons.view());
}

NDArray<double, 2>
AngleCalibration::convert(const std::vector<std::string> &file_list_,
                          const size_t module_index) {
    ssize_t new_num_bins = num_fixed_angle_width_bins();

    LOG(TLogLevel::logDEBUG)
        << fmt::format("num_fixed_angle_width_bins: {}", new_num_bins);

    NDArray<double, 2> weighted_moments(
        std::array<ssize_t, 2>{new_num_bins, 3},
        -1.0); // weighted moments of order 0-2 for each bin
    NDArray<size_t, 1> num_contributing_photons(
        std::array<ssize_t, 1>{new_num_bins},
        0); // number of contributing photons for each bin

    if (module_is_disconnected(module_index)) {
        LOG(TLogLevel::logINFO)
            << fmt::format("module {} is disconnected", module_index);
        return NDArray<double, 2>(std::array<ssize_t, 2>{new_num_bins, 2},
                                  -1.0);
    }

    for (const auto &file : file_list_) {

        LOG(TLogLevel::logDEBUG1) << fmt::format(
            "finished redistributing photon counts for file: {}", file);

        MythenFrame frame = mythen_file_reader->read_frame(file);

        redistribute_photon_counts_to_fixed_angle_width_bins<false>(
            module_index, frame, weighted_moments.view(),
            num_contributing_photons.view());
    }

    return computeErrorandMean(weighted_moments.view(),
                               num_contributing_photons.view());
}

void AngleCalibration::optimize_offset_parameter(
    const size_t module_index,
    const std::vector<MythenFrame> &frames_with_base_peak_overlap,
    const std::vector<MythenFrame> &frames_with_base_peak_overlap_prev_module,
    std::shared_ptr<PlotCalibrationProcess> plot, double delta_parameter) {

    // old reference parameters: // parameters can only change +,-
    // base_peak_roi_width
    double old_reference_left_strip_boundary_angle =
        diffraction_angle_from_BC_parameters(module_index, 0.0, 0.0, 0.0);

    LOG(TLogLevel::logDEBUG) << fmt::format(
        "initial parameter: {}", BCparameters.angle_center_beam(module_index));

    auto objective_function = [this, &module_index,
                               &frames_with_base_peak_overlap,
                               &frames_with_base_peak_overlap_prev_module,
                               &old_reference_left_strip_boundary_angle](
                                  const double x) {
        const double prev_value = BCparameters.angle_center_beam(module_index);
        BCparameters.angle_center_beam(module_index) = x;
        double similarity_of_peaks{};

        double new_reference_left_strip_boundary_angle =
            diffraction_angle_from_BC_parameters(module_index, 0.0, 0.0, 0.0);

        if (std::abs(new_reference_left_strip_boundary_angle -
                     old_reference_left_strip_boundary_angle) >
            base_peak_roi_width) {
            similarity_of_peaks = std::numeric_limits<double>::infinity();
        } else {

            try {
                similarity_of_peaks =
                    calculate_similarity_of_peaks_between_modules(
                        module_index, frames_with_base_peak_overlap,
                        frames_with_base_peak_overlap_prev_module, nullptr);
            } catch (const NoBasePeakOverLapError &e) {
                similarity_of_peaks = std::numeric_limits<
                    double>::infinity(); // if there is no overlap of base
                                         // peaks return infinity for
                                         // similarity of peaks
                                         // - this should lead to a shift
                                         // back into direction of overlap
                                         // in line search
            }
        }
        BCparameters.angle_center_beam(module_index) =
            prev_value; // revert shift
        return similarity_of_peaks;
    };

    const double peak_center_reference_module = center_base_peak(
        module_index - 1, frames_with_base_peak_overlap_prev_module);

    auto difference_peak_centers =
        [this, &module_index, &frames_with_base_peak_overlap,
         &peak_center_reference_module](const double x) {
            const double prev_value =
                BCparameters.angle_center_beam(module_index);
            BCparameters.angle_center_beam(module_index) = x;

            double peak_center_current_module =
                center_base_peak(module_index, frames_with_base_peak_overlap);

            BCparameters.angle_center_beam(module_index) =
                prev_value; // revert shift

            return peak_center_current_module - peak_center_reference_module;
        };

    double previous_similarity_of_peaks =
        calculate_similarity_of_peaks_between_modules<false>(
            module_index, frames_with_base_peak_overlap,
            frames_with_base_peak_overlap_prev_module, nullptr);

    LOG(TLogLevel::logDEBUG) << fmt::format("initial similarity of peaks: {}",
                                            previous_similarity_of_peaks);
    if (plot) {
        calculate_similarity_of_peaks_between_modules<true>(
            module_index, frames_with_base_peak_overlap,
            frames_with_base_peak_overlap_prev_module, plot);

        plot->show();
    }

    double next_similarity_of_peaks{};

    // define search direction
    double diff_peak_centers =
        difference_peak_centers(BCparameters.angle_center_beam(module_index));

    double previous_difference_in_peaks = diff_peak_centers;

    double next_difference_in_peaks{};

    double relative_change_difference_in_peaks{};

    size_t iteration_index = 0;

    if (diff_peak_centers < 0.0) {
        // line search in positive direction
        LOG(TLogLevel::logDEBUG) << "moving peak to the right: ";

        next_difference_in_peaks = difference_peak_centers(
            BCparameters.angle_center_beam(module_index) + delta_parameter);

        // reduce step size if overshoot reference peak
        while (next_difference_in_peaks > 0.0) {
            LOG(TLogLevel::logDEBUG) << fmt::format("reducing step size");

            delta_parameter *= 0.5; // reduce step size
            next_difference_in_peaks = difference_peak_centers(
                BCparameters.angle_center_beam(module_index) + delta_parameter);

            relative_change_difference_in_peaks =
                std::abs(std::abs(next_difference_in_peaks) -
                         std::abs(previous_difference_in_peaks)) /
                (std::abs(previous_difference_in_peaks) +
                 std::numeric_limits<double>::epsilon());
        }

        next_similarity_of_peaks = objective_function(
            BCparameters.angle_center_beam(module_index) + delta_parameter);

        relative_change_difference_in_peaks =
            std::abs(next_difference_in_peaks - previous_difference_in_peaks) /
            (std::abs(previous_difference_in_peaks) +
             std::numeric_limits<double>::epsilon());

        while (next_similarity_of_peaks < previous_similarity_of_peaks ||
               ((std::abs(next_difference_in_peaks) <
                 std::abs(previous_difference_in_peaks)) &&
                (relative_change_difference_in_peaks > 1e-4) &&
                (std::abs(next_difference_in_peaks) > 1e-10))) {
            BCparameters.angle_center_beam(module_index) += delta_parameter;

            LOG(TLogLevel::logDEBUG)
                << fmt::format("Iteration {} with peak similarity of {}, "
                               "difference in peak centers: {}",
                               iteration_index, next_similarity_of_peaks,
                               next_difference_in_peaks);

            ++iteration_index;

            if (plot) {
                calculate_similarity_of_peaks_between_modules<true>(
                    module_index, frames_with_base_peak_overlap,
                    frames_with_base_peak_overlap_prev_module, plot);
                plot->show();
            }

            previous_similarity_of_peaks = next_similarity_of_peaks;
            previous_difference_in_peaks = next_difference_in_peaks;

            next_difference_in_peaks = difference_peak_centers(
                BCparameters.angle_center_beam(module_index) + delta_parameter);

            relative_change_difference_in_peaks =
                std::abs(next_difference_in_peaks -
                         previous_difference_in_peaks) /
                (std::abs(previous_difference_in_peaks) +
                 std::numeric_limits<double>::epsilon());

            // reduce step size if overshoot reference peak
            while (next_difference_in_peaks > 0.0) {
                LOG(TLogLevel::logDEBUG) << fmt::format("reducing step size");
                delta_parameter *= 0.5; // reduce step size
                next_difference_in_peaks = difference_peak_centers(
                    BCparameters.angle_center_beam(module_index) +
                    delta_parameter);
            }

            next_similarity_of_peaks = objective_function(
                BCparameters.angle_center_beam(module_index) + delta_parameter);

            LOG(TLogLevel::logDEBUG1) << fmt::format(
                "previous difference peak centers: {}, next differnce peak "
                "centers: {}",
                previous_difference_in_peaks, next_difference_in_peaks);

            LOG(TLogLevel::logDEBUG1) << fmt::format(
                "previous similarity: {}, next similarity: {}",
                previous_similarity_of_peaks, next_similarity_of_peaks);

            relative_change_difference_in_peaks =
                std::abs(std::abs(next_difference_in_peaks) -
                         std::abs(previous_difference_in_peaks)) /
                (std::abs(previous_difference_in_peaks) +
                 std::numeric_limits<double>::epsilon());
        }
    } else {
        // line search in negative direction
        LOG(TLogLevel::logDEBUG) << "moving peak to the left: ";

        next_difference_in_peaks = difference_peak_centers(
            BCparameters.angle_center_beam(module_index) - delta_parameter);

        // reduce step size if overshoot reference peak
        while (next_difference_in_peaks < 0.0) {
            LOG(TLogLevel::logDEBUG) << fmt::format("reducing step size");
            delta_parameter *= 0.5; // reduce step size
            next_difference_in_peaks = difference_peak_centers(
                BCparameters.angle_center_beam(module_index) - delta_parameter);

            relative_change_difference_in_peaks =
                std::abs(std::abs(next_difference_in_peaks) -
                         std::abs(previous_difference_in_peaks)) /
                (std::abs(previous_difference_in_peaks) +
                 std::numeric_limits<double>::epsilon());
        }

        next_similarity_of_peaks = objective_function(
            BCparameters.angle_center_beam(module_index) - delta_parameter);

        relative_change_difference_in_peaks =
            std::abs(next_difference_in_peaks - previous_difference_in_peaks) /
            (std::abs(previous_difference_in_peaks) +
             std::numeric_limits<double>::epsilon());
        // TODO: probably need convergence criterion
        while (next_similarity_of_peaks < previous_similarity_of_peaks ||
               (std::abs(next_difference_in_peaks) <
                    std::abs(previous_difference_in_peaks) &&
                relative_change_difference_in_peaks > 1e-4 &&
                (std::abs(next_difference_in_peaks) > 1e-10))) {
            BCparameters.angle_center_beam(module_index) -= delta_parameter;

            LOG(TLogLevel::logDEBUG)
                << fmt::format("Iteration {} with peak similarity of {}, "
                               "difference in peak centers: {}",
                               iteration_index, next_similarity_of_peaks,
                               next_difference_in_peaks);

            ++iteration_index;

            if (plot) {
                calculate_similarity_of_peaks_between_modules<true>(
                    module_index, frames_with_base_peak_overlap,
                    frames_with_base_peak_overlap_prev_module, plot);
                plot->show();
            }

            previous_similarity_of_peaks = next_similarity_of_peaks;
            previous_difference_in_peaks = next_difference_in_peaks;

            next_difference_in_peaks = difference_peak_centers(
                BCparameters.angle_center_beam(module_index) - delta_parameter);

            // reduce step size if overshoot reference peak
            while (next_difference_in_peaks < 0.0) {
                LOG(TLogLevel::logDEBUG) << fmt::format("reducing step size");
                delta_parameter *= 0.5; // reduce step size

                next_difference_in_peaks = difference_peak_centers(
                    BCparameters.angle_center_beam(module_index) -
                    delta_parameter);
            }

            next_similarity_of_peaks = objective_function(
                BCparameters.angle_center_beam(module_index) - delta_parameter);

            LOG(TLogLevel::logDEBUG1) << fmt::format(
                "previous difference peak centers: {}, next differnce peak "
                "centers: {}",
                previous_difference_in_peaks, next_difference_in_peaks);

            LOG(TLogLevel::logDEBUG1) << fmt::format(
                "previous similarity: {}, next similarity: {}",
                previous_similarity_of_peaks, next_similarity_of_peaks);

            relative_change_difference_in_peaks =
                std::abs(std::abs(next_difference_in_peaks) -
                         std::abs(previous_difference_in_peaks)) /
                (std::abs(previous_difference_in_peaks) +
                 std::numeric_limits<double>::epsilon());
        }
    }
}

// actually used to optimize BC parameters, first parameters used to
// optimize Lm second parameter used to optimize phi
void AngleCalibration::optimize_coupled_parameters(
    const size_t module_index,
    const std::vector<MythenFrame> &frames_with_base_peak_overlap,
    std::shared_ptr<PlotCalibrationProcess> plot, const double delta_parameter1,
    const double delta_parameter2) {

    constexpr size_t max_iterations = 20;
    bool convergence_criterion = false;
    size_t iteration_index = 0;

    LOG(TLogLevel::logDEBUG) << fmt::format(
        "starting optimization for module {} with initial "
        "parameters: center_sample_distance: {}, "
        "angle_center_module_normal: "
        "{}, angle_center_beam: {}",
        module_index, BCparameters.module_center_sample_distances(module_index),
        BCparameters.angle_center_module_normal(module_index),
        BCparameters.angle_center_beam(module_index));

    // old reference parameters: // parameters can only change +,-
    // base_peak_roi_width
    double old_reference_left_strip_boundary_angle =
        diffraction_angle_from_BC_parameters(module_index, 0.0, 0.0, 0.0);

    double next_similarity_of_peaks{};

    double previous_similarity_of_peaks =
        calculate_similarity_of_peaks_between_acquisitions(
            module_index, frames_with_base_peak_overlap, plot);

    LOG(logDEBUG) << fmt::format("initial similarity of peaks: {}",
                                 previous_similarity_of_peaks);

    auto objective_function1 = [this, &module_index,
                                &frames_with_base_peak_overlap,
                                &old_reference_left_strip_boundary_angle](
                                   const std::array<double, 2> &parameters) {
        const double prev_parameter1 =
            BCparameters.module_center_sample_distances(module_index);
        const double prev_parameter2 =
            BCparameters.angle_center_module_normal(module_index);

        BCparameters.module_center_sample_distances(module_index) =
            parameters[0];
        BCparameters.angle_center_module_normal(module_index) = parameters[1];

        double new_reference_left_strip_boundary_angle =
            diffraction_angle_from_BC_parameters(module_index, 0.0, 0.0, 0.0);

        double similarity_of_peaks{};

        if (std::abs(new_reference_left_strip_boundary_angle -
                     old_reference_left_strip_boundary_angle) >
            base_peak_roi_width) {
            similarity_of_peaks = std::numeric_limits<double>::infinity();
        } else {
            try {
                similarity_of_peaks =
                    calculate_similarity_of_peaks_between_acquisitions(
                        module_index, frames_with_base_peak_overlap, nullptr);
            } catch (const NoBasePeakOverLapError &e) {

                similarity_of_peaks = std::numeric_limits<
                    double>::infinity(); // if there is no overlap of
                                         // base peaks return infinity
                                         // for similarity of peaks
                                         // - this should lead to a
                                         // shift back into direction of
                                         // overlap in line search
            }
        }

        BCparameters.module_center_sample_distances(module_index) =
            prev_parameter1; // revert shift
        BCparameters.angle_center_module_normal(module_index) =
            prev_parameter2; // revert shift

        return similarity_of_peaks;
    };

    while (!convergence_criterion && (iteration_index < max_iterations)) {

        // optimize first two parameters together using newtons method
        auto [optimized_parameters, optimized_similarity,
              optimality_criterion] =
            newton_method(
                delta_parameter1, delta_parameter2,
                {BCparameters.module_center_sample_distances(module_index),
                 BCparameters.angle_center_module_normal(module_index)},
                previous_similarity_of_peaks, objective_function1);

        LOG(TLogLevel::logDEBUG)
            << fmt::format("Iteration {} with peak similarity of {}",
                           iteration_index, next_similarity_of_peaks);

        if (plot) {
            calculate_similarity_of_peaks_between_acquisitions<true>(
                module_index, frames_with_base_peak_overlap, plot);
            plot->show();
        }

        double relative_change =
            std::abs(previous_similarity_of_peaks - next_similarity_of_peaks) /
            previous_similarity_of_peaks;

        double relative_change_in_parameters = std::max(
            std::abs(
                optimized_parameters[0] -
                BCparameters.module_center_sample_distances(module_index)) /
                std::abs(
                    BCparameters.module_center_sample_distances(module_index)),
            std::abs(optimized_parameters[1] -
                     BCparameters.angle_center_module_normal(module_index)) /
                std::abs(
                    BCparameters.angle_center_module_normal(module_index)));

        BCparameters.module_center_sample_distances(module_index) =
            optimized_parameters[0];
        BCparameters.angle_center_module_normal(module_index) =
            optimized_parameters[1];

        next_similarity_of_peaks = optimized_similarity;

        previous_similarity_of_peaks = next_similarity_of_peaks;

        convergence_criterion = (relative_change < 1.0e-5 &&
                                 relative_change_in_parameters < 1.0e-5) ||
                                optimality_criterion < 1.0e-4;

        ++iteration_index;
    }
}

// maybe better L-BFGS !!!!

// TODO: it just writes to a csv/txt file maybe consider using same file
// format as in m_custom_file_ptr but then user needs to implement more
// e.g. custom append
void AngleCalibration::write_DG_parameters_to_file(
    const std::filesystem::path &filename, const DGParameters &parameters) {

    std::ofstream output_file(filename);

    if (!output_file) {
        LOG(angcal::TLogLevel::logERROR) << "Error opening file!" << std::endl;
    }

    output_file.precision(15);

    // output_file << "module, center, conversion, offset" << std::endl;
    // // header

    for (ssize_t module_index = 0; module_index < parameters.num_modules();
         ++module_index) {
        /*
        if constexpr (std::is_same_v<decltype(parameters),
        DGParameters>) { auto [center, conversion, offset] =
                parameters.convert_to_DGParameters(module_index);
        }
        */
        append_to_file(output_file, module_index,
                       parameters.centers(module_index),
                       parameters.conversions(module_index),
                       parameters.offsets(module_index));
    }

    output_file.close();
}

void AngleCalibration::append_to_file(std::ofstream &os,
                                      const size_t module_index,
                                      const double center,
                                      const double conversion,
                                      const double offset) {
    os << "module " << module_index << " center " << center << " +- 0.0000"
       << " conversion " << conversion << " +- 0.0000"
       << " offset " << offset << " +- 0.0000"
       << std::endl; // dummy errors for now
}

} // namespace angcal
