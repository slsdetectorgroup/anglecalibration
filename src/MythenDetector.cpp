#include "MythenDetector.hpp"
#include "logger.hpp"
#include <cmath>
#include <limits>

namespace angcal {

MythenDetector::MythenDetector(const size_t num_modules_)
    : m_num_modules(num_modules_) {}

const DGParameters &MythenDetector::get_DGparameters() const {
    return DGparameters;
}

DGParameters &MythenDetector::get_DGparameters() { return DGparameters; }

void MythenDetector::set_DGparameters(const DGParameters &dgparameters) {
    DGparameters = dgparameters;
}

const BCParameters &MythenDetector::get_BCparameters() const {
    return BCparameters;
}

BCParameters &MythenDetector::get_BCparameters() { return BCparameters; }

void MythenDetector::set_BCparameters(const BCParameters &bcparameters) {
    BCparameters = bcparameters;
}

double MythenDetector::diffraction_angle_from_BC_parameters(
    const size_t module_index, const double detector_angle, size_t strip_index,
    const double distance_to_strip) const {

    const double angle_module_center_normal =
        BCparameters.angle_center_module_normal(module_index);
    const double distance_center_sample =
        BCparameters.module_center_sample_distances(module_index);
    const double angle_module_center_beam =
        BCparameters.angle_center_beam(module_index);

    strip_index = std::signbit(distance_center_sample)
                      ? m_strips_per_module - strip_index - 1
                      : strip_index; // TODO: are the values sored in reserve?

    double diffraction_angle =
        angle_module_center_beam + angle_module_center_normal -
        180.0 / M_PI *
            std::atan(
                (std::abs(distance_center_sample) *
                     std::sin(M_PI / 180.0 * angle_module_center_normal) +
                 (m_strips_per_module * 0.5 - strip_index - distance_to_strip) *
                     m_pitch) /
                (std::abs(distance_center_sample) *
                 std::cos(M_PI / 180.0 * angle_module_center_normal))) +
        elastic_correction(detector_angle) + m_offset +
        m_sample_detector_offset;

    return sample_displacement_correction(diffraction_angle);
}

double MythenDetector::diffraction_angle_from_DG_parameters(
    const size_t module_index, const double detector_angle, size_t strip_index,
    const double distance_to_strip) const {

    double module_offset = DGparameters.offsets(module_index);
    double center = DGparameters.centers(module_index);
    double conversion = DGparameters.conversions(module_index);

    strip_index = std::signbit(conversion)
                      ? m_strips_per_module - strip_index - 1
                      : strip_index; // TODO: are the values sored in reserve?

    double diffraction_angle =
        module_offset +
        180.0 / M_PI *
            (center * std::abs(conversion) -
             std::atan((center - (strip_index + distance_to_strip)) *
                       std::abs(conversion))) +
        elastic_correction(detector_angle) + m_offset +
        m_sample_detector_offset;

    return sample_displacement_correction(diffraction_angle);
}

double MythenDetector::diffraction_angle_from_EE_parameters(
    const double module_center_distance, const double normal_distance,
    const double angle, const double detector_angle, size_t strip_index,
    const double distance_to_strip) const {

    strip_index = std::signbit(normal_distance)
                      ? m_strips_per_module - strip_index - 1
                      : strip_index;

    double diffraction_angle =
        angle -
        180.0 / M_PI *
            std::atan((module_center_distance -
                       m_pitch * (strip_index + distance_to_strip)) /
                      std::abs(normal_distance)) +
        elastic_correction(detector_angle) + m_offset +
        m_sample_detector_offset;
    return sample_displacement_correction(diffraction_angle);
}

double MythenDetector::sample_displacement_correction(
    const double diffraction_angle) const {

    /*
    double sample_displacement_correction =
        std::atan(std::tan(M_PI / 180.0 * diffraction_angle) *
                  std::sqrt(std::pow(sample_x_displacement, 2)
    + std::pow(sample_y_displacement, 2)) /
                  average_distance_sample_pixel);
    */

    double whatever1 =
        std::sin(M_PI / 180.0 * diffraction_angle) +
        m_sample_y_displacement / m_average_distance_sample_pixel;

    double whatever2 =
        std::cos(M_PI / 180.0 * diffraction_angle) -
        m_sample_x_displacement / m_average_distance_sample_pixel;

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

double MythenDetector::angular_strip_width_from_DG_parameters(
    const size_t module_index, const size_t local_strip_index) const {

    return std::abs(diffraction_angle_from_DG_parameters(
                        module_index, 0.0, local_strip_index, 0.0) -
                    diffraction_angle_from_DG_parameters(
                        module_index, 0.0, local_strip_index, 1.0));
}

double MythenDetector::angular_strip_width_from_BC_parameters(
    const size_t module_index, const size_t local_strip_index) const {

    return std::abs(diffraction_angle_from_BC_parameters(
                        module_index, 0.0, local_strip_index, 0.0) -
                    diffraction_angle_from_BC_parameters(
                        module_index, 0.0, local_strip_index, 1.0));
}

double MythenDetector::angular_strip_width_from_EE_parameters(
    const double module_center_distance, const double normal_distance,
    const double angle, const size_t local_strip_index) const {

    return std::abs(diffraction_angle_from_EE_parameters(
                        module_center_distance, normal_distance, angle, 0.0,
                        local_strip_index, 0.0) -
                    diffraction_angle_from_EE_parameters(
                        module_center_distance, normal_distance, angle, 0.0,
                        local_strip_index, 1.0));
}

double MythenDetector::elastic_correction(const double detector_angle) const {
    return detector_angle +
           m_elastic_correction_factor *
               (-std::sin(M_PI / 180.0 *
                          (detector_angle - m_detector_vertical_axis_offset)));
}

double MythenDetector::transverse_width(const size_t module_index,
                                        size_t strip_index) const {

    // convert to EE parameters
    const auto [normal_distance, module_center_distance, angle] =
        BCparameters.convert_to_EEParameters(module_index);

    strip_index = std::signbit(normal_distance)
                      ? m_strips_per_module - 1 - strip_index
                      : strip_index;

    const double distance_sample_pixel =
        std::sqrt(std::pow(normal_distance, 2) +
                  std::pow(module_center_distance - m_pitch * strip_index, 2));

    const double transverse_width =
        (2 *
         std::atan(m_strip_transverse_width / (2 * distance_sample_pixel))) /
        (2 * std::atan(m_strip_transverse_width /
                       (2 * m_average_distance_sample_pixel)));

    return transverse_width;
}

double MythenDetector::solid_angle(const size_t module_index,
                                   size_t strip_index) const {

    // convert to EE parameters
    const auto [normal_distance, module_center_distance, angle] =
        DGparameters.convert_to_EEParameters(module_index);

    strip_index = std::signbit(normal_distance)
                      ? MythenDetector::strips_per_module() - 1 - strip_index
                      : strip_index;

    const double distance_sample_pixel =
        std::sqrt(std::pow(normal_distance, 2) +
                  std::pow(module_center_distance - m_pitch * strip_index, 2));

    // orthogonal projection of strip area onto plane normal to beam
    double projection_strip_area =
        m_pitch * m_strip_transverse_width * std::abs(normal_distance) /
        distance_sample_pixel; // normal_distance/distance_sample_pixel is
                               // the cosine of the angle between the strip
                               // normal and the diffracted beam

    return projection_strip_area / std::pow(distance_sample_pixel, 2);
}

} // namespace angcal