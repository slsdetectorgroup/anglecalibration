#pragma once
#include <cmath>
#include <cstdint>
#include <vector>

#include "Parameters.hpp"

namespace angcal {

/**
 * @brief Mythen detector specifications and geometry (module parameters and
 * methods to calculate diffraction angles)
 */
class MythenDetector {

  public:
    /**
     * @brief Constructor for MythenDetector class
     * @param num_modules_ number of modules in detector
     */
    explicit MythenDetector(const size_t num_modules_ = 48);

    ~MythenDetector() = default;

    size_t num_modules() const { return m_num_modules; }

    ssize_t num_strips() const {
        return static_cast<ssize_t>(m_num_modules * m_strips_per_module);
    }

    static constexpr ssize_t strips_per_module() { return m_strips_per_module; }

    static constexpr double pitch() { return m_pitch; }

    static constexpr double strip_transverse_width() {
        return m_strip_transverse_width;
    }
    /**
     * @brief average distance from sample to pixel [mm]
     */
    double average_distance_sample_pixel() const {
        return m_average_distance_sample_pixel;
    }

    /**
     * @brief set average distance from sample to pixel [mm]
     * @param average_distance_sample_pixel_ average distance from sample to
     * pixel [mm]
     */
    void
    average_distance_sample_pixel(const double average_distance_sample_pixel_) {
        m_average_distance_sample_pixel = average_distance_sample_pixel_;
    }

    /**
     * @brief average solid angle of pixels in Mythen detector [sr]
     */
    double average_solid_angle() const { return m_average_solid_angle; }

    /**
     * @brief set average solid angle of pixels in Mythen detector [sr]
     * @param average_solid_angle_ average solid angle of pixels in Mythen
     * detector [sr]
     */
    void average_solid_angle(const double average_solid_angle_) {
        m_average_solid_angle = average_solid_angle_;
    }

    /**
     * @brief Offset between sample horizontal plane and detector [degrees]
     */
    double sample_detector_offset() const { return m_sample_detector_offset; }

    /**
     * @brief set offset between sample horizontal plane and detector [degrees]
     * @param sample_detector_offset_ offset between sample horizontal plane and
     * detector [degrees]
     */
    void sample_detector_offset(const double sample_detector_offset_) {
        m_sample_detector_offset = sample_detector_offset_;
    }

    /**
     * @brief additional offset to sample detector offset (can change in
     * experimental setup) [degrees]
     */
    double offset() const { return m_offset; }

    /**
     * @brief set additional offset to sample detector offset (can change in
     * experimental setup) [degrees]
     * @param offset_ additional offset to sample detector offset (can change in
     * experimental setup) [degrees]
     */
    void offset(const double offset_) { m_offset = offset_; }

    /**
     * @brief elastic correction factor
     */
    double elastic_correction_factor() const {
        return m_elastic_correction_factor;
    }

    /**
     * @brief set elastic correction factor
     * @param elastic_correction_factor_ elastic correction factor
     */
    void elastic_correction_factor(const double elastic_correction_factor_) {
        m_elastic_correction_factor = elastic_correction_factor_;
    }

    /**
     * @brief offset of detector to vertical axis (used for elastic correction)
     * [degrees]
     */
    double detector_vertical_axis_offset() const {
        return m_detector_vertical_axis_offset;
    }

    /**
     * @brief set offset of detector to vertical axis (used for elastic
     * correction) [degrees]
     * @param detector_vertical_axis_offset_ offset of detector to vertical
     * axis (used for elastic correction) [degrees]
     */
    void
    detector_vertical_axis_offset(const double detector_vertical_axis_offset_) {
        m_detector_vertical_axis_offset = detector_vertical_axis_offset_;
    }

    /**
     * @brief displacement of sample in x direction of horizontal sample plane
     * [mm]
     */
    double sample_x_displacement() const { return m_sample_x_displacement; }

    /**
     * @brief set displacement of sample in x direction of horizontal sample
     * plane [mm]
     * @param sample_x_displacement_ displacement of sample in x direction of
     * horizontal sample plane [mm]
     */
    void sample_x_displacement(const double sample_x_displacement_) {
        m_sample_x_displacement = sample_x_displacement_;
    }

    /**
     * @brief displacement of sample in y direction of horizontal sample plane
     * [mm]
     */
    double sample_y_displacement() const { return m_sample_y_displacement; }

    /**
     * @brief set displacement of sample in y direction of horizontal sample
     * plane [mm]
     * @param sample_y_displacement_ displacement of sample in y direction of
     * horizontal sample plane [mm]
     */
    void sample_y_displacement(const double sample_y_displacement_) {
        m_sample_y_displacement = sample_y_displacement_;
    }

    /**
     * @brief measured dead-time [s]
     */
    double dead_time() const { return m_dead_time; }

    /**
     * @brief set measured dead-time [s]
     * @param dead_time_ measured dead-time [s]
     */
    void dead_time(const double dead_time_) { m_dead_time = dead_time_; }

    /**
     * @brief get the unconnected modules
     */
    std::vector<ssize_t> &unconnected_modules() {
        return m_unconnected_modules;
    }

    /**
     * @brief get the histoic DG parameters
     */
    const DGParameters &get_DGparameters() const; // TODO do I need that?

    DGParameters &get_DGparameters();

    /**
     * @brief set the histoic DG parameters
     */
    void set_DGparameters(const DGParameters &dgparameters);

    /**
     * @brief get the current BC parameters
     */
    const BCParameters &get_BCparameters() const;

    BCParameters &get_BCparameters();

    /**
     * @brief set the current BC parameters
     */
    void set_BCparameters(const BCParameters &bcparameters);

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

    /**
     * @brief calculates the elastic correction factor for given detector angle
     * based on torsional compliance correction model
     * @param detector_angle detector axis motor position [degrees]
     */
    double elastic_correction(const double detector_angle) const;

    /**
     * @brief calculates the transverse width of the strip projected onto a
     * plane
     * @param module_index index of module
     * @param strip_index local strip index of module e.g. 0-1279
     * @return transverse width of strip projected onto a plane [mm] // TODO: is
     * it mm
     */
    double transverse_width(const size_t module_index,
                            size_t strip_index) const;

    /**
     * @brief calculates the solid angle of a strip
     * @param module_index index of module
     * @param strip_index local strip index of module e.g. 0-1279
     * @return solid angle of strip [sr]
     */
    double solid_angle(const size_t module_index, size_t strip_index) const;

  private:
    /**
     * @brief calculates the correction of the diffraction angle due to
     * sample displacement in the x and y direction of the horizontal sample
     * plane
     * @param diffraction_angle diffraction angle calculated at center of
     * strip without sample displacement correction [degrees]
     * @return correction of diffraction angle due to sample displacement
     * [degrees]
     */
    double sample_displacement_correction(const double diffraction_angle) const;

    // memebers

    /// @brief number of strips/channels per module
    static constexpr size_t m_strips_per_module = 1280;

    /// @brief Strip/channel width of Mythen detector [mm]
    static constexpr double m_pitch = 0.05;

    /// @brief Strip/channel height of Mythen detector [mm]
    static constexpr double m_strip_transverse_width = 8.0;

    /// @brief number of modules in detector
    const size_t m_num_modules{};

    /// @brief average euclidean distance between sample and pixel [mm]
    double m_average_distance_sample_pixel =
        2500.0 / M_PI; // TODO why are two values
                       // 4420.97064144153710469121564923651006_DP (R_std_H) in
                       // Antonios code

    /// @brief average solid angle of pixels in Mythen detector [sr]
    double m_average_solid_angle = 6.31654681669718951605407423992073673e-7;
    // 2.0465611686098894032015200537343187d - 8 TODO again two values

    /// @brief Offset between sample horizontal plane and detector
    /// [degrees]
    double m_sample_detector_offset = 1.4715;

    /// @brief additional offset to sample detector offset (can change in
    /// experimental setup) [degrees]
    // offset of detector to sample
    double m_offset = 0.0;

    /// @brief elastic correction factor
    double m_elastic_correction_factor = 0.0;

    /// @brief offste of detector to vertical axis (used for elastic correction)
    /// [degrees]
    double m_detector_vertical_axis_offset = 0.0;

    /// @brief displacement of sample in x direction of horizontal sample plane
    /// [mm]
    double m_sample_x_displacement = 0.0;

    /// @brief displacement of sample in y direction of horizontal sample plane
    /// [mm]
    double m_sample_y_displacement = 0.0;

    /// @brief measured dead-time [s]
    double m_dead_time =
        76.08e-9; // measured dead-time 31.8E-9_DP, 95.E-9_DP other beam
                  // settings ! TODO: should it be configurable?

    /// @brief list of module indices which are unconnected
    std::vector<ssize_t> m_unconnected_modules{};

    DGParameters DGparameters;

    /// @brief mythen detector BC parameters (best computing) for current
    /// detector geometry
    BCParameters BCparameters;
};

} // namespace angcal