#include "FlatField.hpp"
#include "logger.hpp"

using namespace aare;
namespace angcal {

FlatField::FlatField(
    std::shared_ptr<MythenDetectorSpecifications> mythen_detector_)
    : mythen_detector(mythen_detector_),
      DGparameters(mythen_detector->max_modules) {

    bad_channels =
        NDArray<bool, 1>(std::array<ssize_t, 1>{mythen_detector->num_strips()},
                         false); // default initialized to false
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

double FlatField::elastic_correction(const double detector_angle) const {
    return detector_angle +
           mythen_detector->elastic_correction_factor *
               (-std::sin(M_PI / 180.0 *
                          (detector_angle -
                           mythen_detector->detector_vertical_axis_offset)));
}

double FlatField::diffraction_angle_from_DG_parameters(
    const size_t module_index, const double detector_angle, size_t strip_index,
    const double distance_to_strip) const {

    double offset = DGparameters.offsets(module_index);
    double center = DGparameters.centers(module_index);
    double conversion = DGparameters.conversions(module_index);

    strip_index =
        std::signbit(conversion)
            ? MythenDetectorSpecifications::strips_per_module - 1 - strip_index
            : strip_index; // TODO: are the values sored in reserve?

    return offset +
           180.0 / M_PI *
               (center * std::abs(conversion) -
                std::atan((center - (strip_index + distance_to_strip)) *
                          std::abs(conversion))) +
           elastic_correction(detector_angle) + mythen_detector->offset +
           mythen_detector->sample_detector_offset;
}

double FlatField::solid_angle(const size_t module_index,
                              size_t strip_index) const {

    // convert to EE parameters
    const auto [normal_distance, module_center_distance, angle] =
        DGparameters.convert_to_EEParameters(module_index);

    strip_index =
        std::signbit(normal_distance)
            ? MythenDetectorSpecifications::strips_per_module - 1 - strip_index
            : strip_index;

    const double distance_sample_pixel = std::sqrt(
        std::pow(normal_distance, 2) +
        std::pow(module_center_distance - mythen_detector->pitch * strip_index,
                 2));

    // orthogonal projection of strip area onto plane normal to beam
    double projection_strip_area =
        mythen_detector->pitch * mythen_detector->transverse_width *
        std::abs(normal_distance) /
        distance_sample_pixel; // normal_distance/distance_sample_pixel is
                               // the cosine of the angle between the strip
                               // normal and the diffracted beam

    return projection_strip_area / std::pow(distance_sample_pixel, 2);
}

/*
double FlatField::transverse_width_correction_factor(const size_t module_index,
                                                     size_t strip_index) const {

    // convert to EE parameters
    const auto [normal_distance, module_center_distance, angle] =
        DGparameters.convert_to_EEParameters(module_index);

    strip_index =
        std::signbit(normal_distance)
            ? MythenDetectorSpecifications::strips_per_module - strip_index - 1
            : strip_index;

    const double distance_sample_pixel = std::sqrt(
        std::pow(normal_distance, 2) +
        std::pow(module_center_distance - mythen_detector->pitch * strip_index,
                 2));

    const double transverse_width_correction_factor =
        (2 * std::atan(mythen_detector->transverse_width /
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

void FlatField::create_normalized_flatfield_from_filelist(
    const std::vector<std::filesystem::path> &filelist,
    std::shared_ptr<MythenFileReader> file_reader) {

    flat_field = NDArray<double, 2>(
        std::array<ssize_t, 2>{
            static_cast<ssize_t>(mythen_detector->num_strips()), 2},
        -1.0);

    NDArray<double, 1> sum_statistical_weights(
        std::array<ssize_t, 1>{mythen_detector->num_strips()}, 0.0);

    constexpr uint32_t mighells_correction_constant = 1;

    size_t file_index = 0;

    for (const auto &file : filelist) {
        auto frame = file_reader->read_frame(file);

        uint64_t incident_intensity = frame.incident_intensity;

        LOG(TLogLevel::logDEBUG1)
            << fmt::format("incident intensity {}", incident_intensity);

        if (frame.incident_intensity == 0) {
            throw std::runtime_error(
                fmt::format("incident intensity is zero for file {}. "
                            "Cannot apply I0 correction",
                            file.string()));
        }

        // IncidentIntensityLogFile.append(
        // fmt::format("{}\n", incident_intensity));

        const double I0_correction_factor =
            scale_factor / static_cast<double>(incident_intensity);

        LOG(TLogLevel::logDEBUG1)
            << fmt::format("detector angle: {}", frame.detector_angle);

        for (ssize_t strip_index = 0;
             strip_index < mythen_detector->num_strips(); ++strip_index) {
            if (bad_channels(strip_index)) {
                continue; // skip bad channels
            }

            const size_t module_index =
                strip_index / MythenDetectorSpecifications::strips_per_module;
            const size_t local_strip_index =
                strip_index % MythenDetectorSpecifications::strips_per_module;

            const double left_strip_boundary_angle =
                diffraction_angle_from_DG_parameters(
                    module_index, frame.detector_angle, local_strip_index,
                    -0.5); // left strip boundary in angles [degrees]

            const double right_strip_boundary_angle =
                diffraction_angle_from_DG_parameters(
                    module_index, frame.detector_angle, local_strip_index,
                    +0.5); // right strip boundary in angles [degrees]

            // skip if channel is completely outside of soft window
            if (right_strip_boundary_angle < soft_window.first ||
                left_strip_boundary_angle > soft_window.second) {
                continue;
            }

            // mighells correction
            double photon_counts = frame.photon_counts(strip_index) != 0
                                       ? frame.photon_counts(strip_index) +
                                             mighells_correction_constant
                                       : 0.0;

            double variance_photon_counts = photon_counts; // poisson statistics

            if (photon_counts <= std::numeric_limits<double>::epsilon()) {
                // bad_channels(strip_index) = true; // mark as bad channel
                flat_field(strip_index, 0) = 0.0;
                flat_field(strip_index, 1) = 0.0;
                continue; // skip bad channels
            }

            // incident intensity correction
            photon_counts *= I0_correction_factor;
            variance_photon_counts *=
                I0_correction_factor * I0_correction_factor;

            // scale by soft_window coverage
            double soft_window_coverage =
                (std::min(right_strip_boundary_angle, soft_window.second) -
                 std::max(left_strip_boundary_angle, soft_window.first)) /
                (right_strip_boundary_angle - left_strip_boundary_angle);

            // StripWidthLogFile.append(
            //  fmt::format("{}\n", right_strip_boundary_angle -
            //  left_strip_boundary_angle));

            // CoverageLogFile.append(fmt::format("{}\n",
            // soft_window_coverage));

            LOG(TLogLevel::logDEBUG1)
                << fmt::format("strip {}, soft_window_coverage {}", strip_index,
                               soft_window_coverage);

            if (flat_field(strip_index, 0) == -1.0) {
                flat_field(strip_index, 0) =
                    0.0; // initialize to zero if not set yet
                flat_field(strip_index, 1) =
                    0.0; // initialize to zero if not set yet
            }
            flat_field(strip_index, 0) += photon_counts * soft_window_coverage;
            flat_field(strip_index, 1) += variance_photon_counts *
                                          soft_window_coverage *
                                          soft_window_coverage;
        }
        ++file_index;
    }

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
            double solid_angle_ = solid_angle(
                strip_index / MythenDetectorSpecifications::strips_per_module,
                strip_index % MythenDetectorSpecifications::strips_per_module);

            double solid_angle_correction =
                mythen_detector->average_solid_angle / solid_angle_;

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

void FlatField::read_module_parameters_from_file(
    const std::filesystem::path &filename,
    const std::shared_ptr<SimpleFileInterface> file_reader) {

    file_reader->open(filename);
    file_reader->read_into(DGparameters.parameters.buffer(), 8);

    file_reader->close();
}

void FlatField::read_bad_channels_from_file(
    const std::filesystem::path &filename,
    const std::shared_ptr<SimpleFileInterface> file_reader) {

    file_reader->open(filename);
    file_reader->read_into(reinterpret_cast<std::byte *>(bad_channels.data()));
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

NDView<double, 2> FlatField::get_normalized_flatfield_view() const {
    return flat_field.view();
}

} // namespace angcal