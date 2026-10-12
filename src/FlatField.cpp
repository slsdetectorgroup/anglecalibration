#include "FlatField.hpp"
#include "logger.hpp"

using namespace aare;
namespace angcal {

FlatField::FlatField(std::shared_ptr<MythenDetector> mythen_detector_)
    : mythen_detector(mythen_detector_) {

    bad_channels =
        NDArray<bool, 1>(std::array<ssize_t, 1>{mythen_detector->num_strips()},
                         false); // default initialized to false

    flat_field = NDArray<double, 2>(
        std::array<ssize_t, 2>{
            static_cast<ssize_t>(mythen_detector->num_strips()), 2},
        -1.0);
}

void FlatField::set_scale_factor(const double scale_factor_) {
    scale_factor = scale_factor_;
}

double FlatField::get_scale_factor() const { return scale_factor; }

void FlatField::set_soft_window(const std::pair<double, double> soft_window_) {
    soft_window = soft_window_;
}

std::pair<double, double> FlatField::get_soft_window() const {
    return soft_window;
}

/*
double FlatField::transverse_width_correction_factor(const size_t module_index,
                                                     size_t strip_index) const {

    // convert to EE parameters
    const auto [normal_distance, module_center_distance, angle] =
        DGparameters.convert_to_EEParameters(module_index);

    strip_index =
        std::signbit(normal_distance)
            ? MythenDetector::strips_per_module - strip_index - 1
            : strip_index;

    const double distance_sample_pixel = std::sqrt(
        std::pow(normal_distance, 2) +
        std::pow(module_center_distance - mythen_detector->pitch * strip_index,
                 2));

    const double transverse_width_correction_factor =
        (2 * std::atan(mythen_detector->strip_transverse_width /
                       (2 * distance_sample_pixel)));

    return transverse_width_correction_factor;
}

double FlatField::Antonio_solid_angle_correction_factor(
    const size_t module_index, const size_t strip_index) const {

    double transwidth_correction_factor =
        transverse_width_correction_factor(module_index, strip_index);

    double angular_strip_width =
        std::abs(diffraction_angle_from_DG_parameters(module_index, 0.0,
                                                      strip_index, 0.5) -
                 diffraction_angle_from_DG_parameters(module_index, 0.0,
                                                      strip_index, -0.5));

    double solid_angle =
        M_PI / 180.0 * angular_strip_width * transwidth_correction_factor;

    // scale solid angle to unit
    double solid_angle_correction =
        mythen_detector->average_solid_angle / solid_angle;

    return solid_angle_correction;
}
*/

// TODO: should be 2d if multiple counters enabled
void FlatField::add_flatfield(const NDView<uint32_t, 1> photon_counts,
                              const uint32_t incident_intensity,
                              const double detector_angle) {

    if (incident_intensity == 0) {
        throw std::runtime_error("incident intensity is zero for file {}. "
                                 "Cannot apply I0 correction");
    }

    const double I0_correction_factor =
        scale_factor / static_cast<double>(incident_intensity);

    constexpr uint32_t mighells_correction_constant = 1;

    for (ssize_t strip_index = 0; strip_index < mythen_detector->num_strips();
         ++strip_index) {

        if (bad_channels(strip_index)) {
            continue; // skip bad channels
        }

        const size_t module_index =
            strip_index / MythenDetector::strips_per_module();
        const size_t local_strip_index =
            strip_index % MythenDetector::strips_per_module();

        const double left_strip_boundary_angle =
            mythen_detector->diffraction_angle_from_DG_parameters(
                module_index, detector_angle, local_strip_index,
                0.0); // left strip boundary in angles [degrees]

        const double right_strip_boundary_angle =
            mythen_detector->diffraction_angle_from_DG_parameters(
                module_index, detector_angle, local_strip_index,
                1.0); // right strip boundary in angles [degrees]

        // skip if channel is completely outside of soft window
        if (right_strip_boundary_angle < soft_window.first ||
            left_strip_boundary_angle > soft_window.second) {
            continue;
        }

        // mighells correction
        double photon_count =
            photon_counts(strip_index) != 0
                ? photon_counts(strip_index) + mighells_correction_constant
                : 0.0;

        double variance_photon_count = photon_count; // poisson statistics

        if (photon_count <= std::numeric_limits<double>::epsilon()) {
            // bad_channels(strip_index) = true; // mark as bad channel
            flat_field(strip_index, 0) = 0.0;
            flat_field(strip_index, 1) = 0.0;
            continue; // skip bad channels
        }

        // incident intensity correction
        photon_count *= I0_correction_factor;
        variance_photon_count *= I0_correction_factor * I0_correction_factor;

        // scale by soft_window coverage
        double soft_window_coverage =
            (std::min(right_strip_boundary_angle, soft_window.second) -
             std::max(left_strip_boundary_angle, soft_window.first)) /
            (right_strip_boundary_angle - left_strip_boundary_angle);

        LOG(TLogLevel::logDEBUG1)
            << fmt::format("strip {}, soft_window_coverage {}", strip_index,
                           soft_window_coverage);

        if (flat_field(strip_index, 0) == -1.0) {
            flat_field(strip_index, 0) =
                0.0; // initialize to zero if not set yet
            flat_field(strip_index, 1) =
                0.0; // initialize to zero if not set yet
        }
        flat_field(strip_index, 0) += photon_count * soft_window_coverage;
        flat_field(strip_index, 1) +=
            variance_photon_count * soft_window_coverage * soft_window_coverage;
    }
}

void FlatField::normalize_flatfield() {

    double sum_good_strips = 0.0; // for normalization
    size_t num_good_strips = 0;   // for normalization

    for (ssize_t strip_index = 0; strip_index < mythen_detector->num_strips();
         ++strip_index) {

        if (bad_channels(strip_index) ||
            flat_field(strip_index, 0) <
                std::numeric_limits<double>::epsilon()) {
            continue;
        } else {

            // solid angle correction
            double solid_angle_ = mythen_detector->solid_angle(
                strip_index / MythenDetector::strips_per_module(),
                strip_index % MythenDetector::strips_per_module());

            double solid_angle_correction =
                mythen_detector->average_solid_angle() / solid_angle_;
            flat_field(strip_index, 0) *= solid_angle_correction;

            flat_field(strip_index, 1) *=
                solid_angle_correction *
                solid_angle_correction; // error propagation

            flat_field(strip_index, 1) =
                std::sqrt(flat_field(strip_index, 1)); // store std

            sum_good_strips += flat_field(strip_index, 0);
            LOG(TLogLevel::logDEBUG1) << fmt::format(
                "strip {}: sum_good_strips  {}", strip_index, sum_good_strips);
            ++num_good_strips;
        }
    }

    // normalize flatfield to one
    double normalization_factor = sum_good_strips / num_good_strips;
    if (num_good_strips == 0) {
        throw std::runtime_error(
            "Only bad strips found for flatfield creation.");
    }
    LOG(TLogLevel::logDEBUG) << fmt::format(
        "normalization factor for flatfield: {}", normalization_factor);
    for (ssize_t strip_index = 0; strip_index < mythen_detector->num_strips();
         ++strip_index) {
        if (flat_field(strip_index, 0) != -1.0) {
            flat_field(strip_index, 0) /= normalization_factor;
            flat_field(strip_index, 1) /= normalization_factor;
        }
    }
}

void FlatField::read_normalized_flatfield_from_file(
    const std::string &filename,
    const std::shared_ptr<SimpleFileInterface> file_reader) {
    flat_field = NDArray<double, 2>(std::array<ssize_t, 2>{
        static_cast<ssize_t>(mythen_detector->num_strips()), 2});
    file_reader->open(filename);
    file_reader->read_into(reinterpret_cast<std::byte *>(flat_field.data()),
                           8); // TODO: should be int though !!!

    file_reader->close();
}

void FlatField::read_bad_channels_from_file(
    const std::filesystem::path &filename,
    const std::shared_ptr<SimpleFileInterface> file_reader) {

    file_reader->open(filename);
    file_reader->read_into(reinterpret_cast<std::byte *>(bad_channels.data()));
    file_reader->close();
}

void FlatField::set_bad_channels(const NDArray<bool, 1> &bad_channels_) {
    bad_channels = bad_channels_;
}

NDArray<bool, 1> FlatField::get_bad_channels() const { return bad_channels; }

void FlatField::set_normalized_flatfield(
    const NDArray<double, 2> &flat_field_) {
    flat_field = flat_field_;
}

void FlatField::set_normalized_flatfield(NDArray<double, 2> &&flat_field_) {
    flat_field = std::move(flat_field_);
}

// TODO: should it be view? or copy?
NDArray<double, 2> FlatField::get_normalized_flatfield() const {
    return flat_field;
}

NDView<const double, 2> FlatField::get_normalized_flatfield_view() const {
    return flat_field.view();
}

} // namespace angcal