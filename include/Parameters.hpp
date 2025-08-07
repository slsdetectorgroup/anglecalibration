
#include <cmath>
#include <cstdint>

#include "MythenDetectorSpecifications.hpp"
#include "aare/NDArray.hpp"

namespace angcal {

/**
 * geometric parameters
 */
struct EEParameters {

    EEParameters() = default;

    EEParameters(const ssize_t n_modules) {
        parameters =
            aare::NDArray<double, 2>(std::array<ssize_t, 2>{n_modules, 3});
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

    aare::NDArray<double, 2> parameters{};
};

/**
 * best computing parameters
 */
struct BCParameters {

    BCParameters() = default;

    BCParameters(const ssize_t n_modules) {
        parameters =
            aare::NDArray<double, 2>(std::array<ssize_t, 2>{n_modules, 3});
    }

    double &operator()(const size_t module_index,
                       const size_t parameter_index) {
        return parameters(module_index, parameter_index);
    }

    double operator()(const size_t module_index,
                      const size_t parameter_index) const {
        return parameters(module_index, parameter_index);
    }

    /** angle between center of module and module normal (from sample) (delta)
     */
    double &angle_center_module_normal(const size_t module_index) {
        return parameters(module_index, 0);
    }

    double angle_center_module_normal(const size_t module_index) const {
        return parameters(module_index, 0);
    }

    /**
     * euclidean distance between center of module and sample (L)
     */
    double &module_center_sample_distances(const size_t module_index) {
        return parameters(module_index, 1);
    }
    double module_center_sample_distances(const size_t module_index) const {
        return parameters(module_index, 1);
    }

    /**
     * diffraction angle betwen center module and beam  (phi)
     */
    double &angle_center_beam(const size_t module_index) {
        return parameters(module_index, 2);
    }

    double angle_center_beam(const size_t module_index) const {
        return parameters(module_index, 2);
    }

    aare::NDArray<double, 2> parameters{};
};

// TODO add units
/**
 * historical Detector Group parameters
 */
struct DGParameters {

    DGParameters() = default;

    DGParameters(const ssize_t n_modules) {
        parameters =
            aare::NDArray<double, 2>(std::array<ssize_t, 2>{n_modules, 3});
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

    std::tuple<double, double, double>
    convert_to_BCParameters(const size_t module_index) const {
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
                         0.5 *
                             MythenDetectorSpecifications::strips_per_module()),
                    2));
        double angle_beam_module_center =
            offsets(module_index) +
            180.0 / M_PI * std::abs(conversions(module_index)) *
                centers(module_index) -
            angle_center_module_normal;

        return std::make_tuple(angle_center_module_normal,
                               distance_module_center_sample,
                               angle_beam_module_center);
    }

    void convert_to_BCParameters(BCParameters &bcparameters) const {
        for (ssize_t module_index = 0;
             module_index < bcparameters.parameters.shape(0); ++module_index) {
            auto [angle_center_module_normal, distance_module_center_sample,
                  angle_beam_module_center] =
                convert_to_BCParameters(module_index);
            bcparameters.angle_center_module_normal(module_index) =
                angle_center_module_normal;
            bcparameters.module_center_sample_distances(module_index) =
                distance_module_center_sample;
            bcparameters.angle_center_beam(module_index) =
                angle_beam_module_center;
        }
    }

    std::tuple<double, double, double>
    convert_to_EEParameters(const size_t module_index) const {

        const double module_center_distance =
            centers(module_index) * MythenDetectorSpecifications::pitch();
        const double normal_distance = MythenDetectorSpecifications::pitch() /
                                       std::abs(conversions(module_index));
        const double angle =
            offsets(module_index) + 180.0 / M_PI * centers(module_index) *
                                        std::abs(conversions(module_index));

        return std::make_tuple(module_center_distance, normal_distance, angle);
    }

    EEParameters convert_to_EEParameters() const {

        EEParameters EEparameters(parameters.shape(0));

        for (ssize_t i = 0; i < parameters.shape(0); ++i) {
            auto [module_center_distance, normal_distance, angle] =
                convert_to_EEParameters(i);
            EEparameters.normal_distances(i) = normal_distance;
            EEparameters.module_center_distances(i) = module_center_distance;
            EEparameters.angles(i) = angle;
        }

        return EEparameters;
    }

    aare::NDArray<double, 2> parameters{};
};

} // namespace angcal