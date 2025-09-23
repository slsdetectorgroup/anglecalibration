
#include <cmath>
#include <cstdint>

#include "MythenDetectorSpecifications.hpp"
#include "aare/NDArray.hpp"

namespace angcal {

struct DGParameters; // forward declaration

/**
 * base class
 */
struct Parameters {

    ~Parameters() = default;

    Parameters() = default;

    Parameters(const ssize_t n_modules) {
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

    ssize_t num_modules() const { return parameters.shape(0); }

    aare::NDArray<double, 2> parameters{};
};

// TODO add abstract base class for parameters
/**
 * geometric parameters
 */
struct EEParameters : public Parameters {

    ~EEParameters() = default;

    EEParameters() = default;

    EEParameters(const ssize_t n_modules) : Parameters(n_modules) {};

    /**
     * normal distance between sample and detector (R)
     */
    double &normal_distances(const size_t module_index) {
        return Parameters::operator()(module_index, 0);
    }

    double normal_distances(const size_t module_index) const {
        return Parameters::operator()(module_index, 0);
    }

    /**
     * distances between intersection point of sample normal and module origin
     * (D)
     */
    double &module_center_distances(const size_t module_index) {
        return Parameters::operator()(module_index, 1);
    }
    double module_center_distances(const size_t module_index) const {
        return Parameters::operator()(module_index, 1);
    }

    /** angles between undiffracted beam and orthogonal sample projection on
     * detector (phi)
     */
    double &angles(const size_t module_index) {
        return Parameters::operator()(module_index, 2);
    }

    double angles(const size_t module_index) const {
        return Parameters::operator()(module_index, 0);
    }
};

/**
 * best computing parameters
 */
struct BCParameters : public Parameters {

    ~BCParameters() = default;

    BCParameters() = default;

    BCParameters(const ssize_t n_modules) : Parameters(n_modules) {};

    /** angle between center of module and module normal (from sample) [degrees]
     * (delta)
     */
    // TODO: check the order of filling - for optimization
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
     * diffraction angle betwen center module and beam  (phi) [degrees]
     */
    double &angle_center_beam(const size_t module_index) {
        return parameters(module_index, 2);
    }

    double angle_center_beam(const size_t module_index) const {
        return parameters(module_index, 2);
    }

    std::tuple<double, double, double>
    convert_to_DGParameters(const size_t module_index) const;

    void convert_to_DGParameters(DGParameters &dgparameters) const;

    std::tuple<double, double, double>
    convert_to_EEParameters(const ssize_t module_index) const;

    void convert_to_EEParameters(EEParameters &eeparameters) const;
};

// TODO add units
/**
 * historical Detector Group parameters
 */
struct DGParameters : public Parameters {

    ~DGParameters() = default;

    DGParameters() = default;

    DGParameters(const ssize_t n_modules) : Parameters(n_modules) {};

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
    convert_to_BCParameters(const size_t module_index) const;

    void convert_to_BCParameters(BCParameters &bcparameters) const;

    std::tuple<double, double, double>
    convert_to_EEParameters(const size_t module_index) const;

    void convert_to_EEParameters(EEParameters &eeparameters) const;
};

} // namespace angcal