
#include "Parameters.hpp"
#include "MythenDetectorSpecifications.hpp"
#include <cmath>
#include <tuple>

namespace angcal {
std::tuple<double, double, double>
BCParameters::convert_to_DGParameters(const size_t module_index) const {
    double center =
        std::abs(module_center_sample_distances(module_index)) *
            std::sin(M_PI / 180.0 * angle_center_module_normal(module_index)) /
            MythenDetectorSpecifications::pitch() +
        (MythenDetectorSpecifications::strips_per_module() * 0.5);

    double conversion =
        MythenDetectorSpecifications::pitch() /
        (std::abs(module_center_sample_distances(module_index)) *
         std::cos(M_PI / 180.0 * angle_center_module_normal(module_index)));

    if (std::signbit(module_center_sample_distances(module_index)))
        conversion = -conversion;

    double offset =
        angle_center_beam(module_index) +
        angle_center_module_normal(module_index) -
        180.0 / M_PI *
            (std::tan(M_PI / 180.0 * angle_center_module_normal(module_index)) +
             MythenDetectorSpecifications::strips_per_module() * 0.5 *
                 std::abs(conversion));

    return std::make_tuple(center, conversion, offset);
}

void BCParameters::convert_to_DGParameters(DGParameters &dgparameters) const {
    for (ssize_t module_index = 0;
         module_index < dgparameters.parameters.shape(0); ++module_index) {
        auto [center, conversion, offset] =
            convert_to_DGParameters(module_index);
        dgparameters.centers(module_index) = center;
        dgparameters.conversions(module_index) = conversion;
        dgparameters.offsets(module_index) = offset;
    }
}

std::tuple<double, double, double>
BCParameters::convert_to_EEParameters(const ssize_t module_index) const {
    double angle = angle_center_beam(module_index) +
                   angle_center_module_normal(module_index);

    double module_center_distance =
        std::abs(module_center_sample_distances(module_index)) *
            std::sin(M_PI / 180.0 * angle_center_module_normal(module_index)) +
        (MythenDetectorSpecifications::strips_per_module() * 0.5) *
            MythenDetectorSpecifications::pitch();

    double normal_distance =
        std::abs(module_center_sample_distances(module_index)) *
        std::cos(M_PI / 180.0 * angle_center_module_normal(module_index));

    if (std::signbit(module_center_sample_distances(module_index)))
        normal_distance = -normal_distance;

    return std::make_tuple(normal_distance, module_center_distance, angle);
}

void BCParameters::convert_to_EEParameters(EEParameters &eeparameters) const {
    for (ssize_t i = 0; i < parameters.shape(0); ++i) {
        auto [normal_distance, module_center_distance, angle] =
            convert_to_EEParameters(i);
        eeparameters.normal_distances(i) = normal_distance;
        eeparameters.module_center_distances(i) = module_center_distance;
        eeparameters.angles(i) = angle;
    }
}

std::tuple<double, double, double>
DGParameters::convert_to_BCParameters(const size_t module_index) const {

    double angle_center_module_normal =
        180.0 / M_PI *
        atan((centers(module_index) -
              0.5 * MythenDetectorSpecifications::strips_per_module()) *
             std::abs(conversions(
                 module_index))); // TODO in Antonios code it is minus?
    double distance_module_center_sample =
        MythenDetectorSpecifications::pitch() /
        std::abs(conversions(module_index)) *
        std::sqrt(
            1 +
            std::pow(
                std::abs(conversions(module_index)) *
                    (centers(module_index) -
                     0.5 * MythenDetectorSpecifications::strips_per_module()),
                2));

    distance_module_center_sample = std::signbit(conversions(module_index))
                                        ? -distance_module_center_sample
                                        : distance_module_center_sample;

    double angle_beam_module_center = offsets(module_index) +
                                      180.0 / M_PI *
                                          std::abs(conversions(module_index)) *
                                          centers(module_index) -
                                      angle_center_module_normal;

    return std::make_tuple(angle_center_module_normal,
                           distance_module_center_sample,
                           angle_beam_module_center);
}

void DGParameters::convert_to_BCParameters(BCParameters &bcparameters) const {
    for (ssize_t module_index = 0;
         module_index < bcparameters.parameters.shape(0); ++module_index) {
        auto [angle_center_module_normal, distance_module_center_sample,
              angle_beam_module_center] = convert_to_BCParameters(module_index);
        bcparameters.angle_center_module_normal(module_index) =
            angle_center_module_normal;
        bcparameters.module_center_sample_distances(module_index) =
            distance_module_center_sample;
        bcparameters.angle_center_beam(module_index) = angle_beam_module_center;
    }
}

std::tuple<double, double, double>
DGParameters::convert_to_EEParameters(const size_t module_index) const {

    const double module_center_distance =
        centers(module_index) * MythenDetectorSpecifications::pitch();
    double normal_distance = MythenDetectorSpecifications::pitch() /
                             std::abs(conversions(module_index));

    if (std::signbit(conversions(module_index))) {
        normal_distance = -normal_distance;
    }
    const double angle =
        offsets(module_index) + 180.0 / M_PI * centers(module_index) *
                                    std::abs(conversions(module_index));

    return std::make_tuple(module_center_distance, normal_distance, angle);
}

void DGParameters::convert_to_EEParameters(EEParameters &eeparameters) const {

    for (ssize_t i = 0; i < parameters.shape(0); ++i) {
        auto [module_center_distance, normal_distance, angle] =
            convert_to_EEParameters(i);
        eeparameters.normal_distances(i) = normal_distance;
        eeparameters.module_center_distances(i) = module_center_distance;
        eeparameters.angles(i) = angle;
    }
}
} // namespace angcal