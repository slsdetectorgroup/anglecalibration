#pragma once
#include "AngleCalibration.hpp"
#include "matplotlibcpp.h"
#include <numeric>

namespace angcal {

template <typename T>
void plot(const std::string &plot_title, const NDView<T, 1> y) {

    matplotlibcpp::backend("TkAgg");
    // matplotlibcpp::figure();
    matplotlibcpp::figure_size(800, 600);
    matplotlibcpp::xlabel("Diffraction Angle [degree]");
    matplotlibcpp::ylabel("Photon Counts");
    matplotlibcpp::title(plot_title);

    std::vector<double> bins(y.size());
    std::iota(bins.begin(), bins.end(), 0.0);

    std::vector<double> y_vec;
    y_vec.assign(y.begin(), y.end());

    matplotlibcpp::plot(bins, y_vec);

    matplotlibcpp::show(false);
}

class PlotHelper {

  public:
    /**
     * @brief constructor for PlotHelper class
     * @param anglecalibration
     */
    PlotHelper(std::shared_ptr<AngleCalibration> anglecalibration);

    ~PlotHelper();

    /// @brief pause the program until user input, e.g. before destroying all
    /// plots
    static void pause();

    /// @brief overwrite plot for plotting multiple frames
    /// in one figure
    void overwrite_plot();

    void create_plot() const;

    void initialize_axis(const std::string &plot_title) const;

    /*
    void plot_module_redistributed_to_fixed_angle_width_bins(
        const std::filesystem::path &file_path,
        std::optional<size_t> module_index = std::nullopt) const;
    */

    /**  @brief plot diffraction pattern for given photon counts
     * @param motor_position The position of the motor for resulting
     * diffraction pattern
     * @param photon_counts resulting photon counts of diffraction pattern
     * @param module_index optional module index for plotting only specific
     * module region
     */
    void plot_diffraction_pattern(
        const double motor_position, const NDView<double, 1> &photon_counts,
        std::optional<size_t> module_index = std::nullopt) const;

    /**
     * @brief plot base peak for given photon counts
     * @param photon_counts resulting photon counts of diffraction pattern
     * @param module_index optional module index if only one module contributed
     * to diffraction pattern
     */
    void plot_base_peak(
        const NDView<double, 1> photon_counts,
        const std::optional<size_t> module_index = std::nullopt) const;

  private:
    inline size_t left_module_boundary_as_fixed_angle_bin_index(
        const size_t module_index, const double detector_angle) const {

        return (m_anglecalibration->diffraction_angle_from_DG_parameters(
                    module_index, detector_angle, 0, -0.5) -
                m_anglecalibration->get_angular_range().first) /
               m_anglecalibration->get_histogram_bin_width();
    }

    inline size_t right_module_boundary_as_fixed_angle_bin_index(
        const size_t module_index, const double detector_angle) const {

        return (m_anglecalibration->diffraction_angle_from_DG_parameters(
                    module_index, detector_angle,
                    m_anglecalibration->get_detector_specifications()
                        ->strips_per_module,
                    +0.5) -
                m_anglecalibration->get_angular_range().first) /
               m_anglecalibration->get_histogram_bin_width();
    }

  private:
    /// @brief shared pointer to AngleCalibration class for accessing
    /// calibration specific parameters and functions
    const std::shared_ptr<AngleCalibration> m_anglecalibration{};

    /// @brief function to convert bin index to diffraction angle
    std::function<double(const size_t)> bin_to_diffraction_angle{};

    /// @brief function to convert bin index to diffraction angle for base peak
    /// ROI only
    std::function<double(const size_t)>
        bin_to_diffraction_angle_base_peak_ROI_only{};

    bool m_overwrite_plot = false;
};

} // namespace angcal