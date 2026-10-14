#pragma once
#include "AngleCalibration.hpp"
#include "matplotlibcpp.h"
#include <atomic>
#include <numeric>
#include <thread>

namespace angcal {

/*
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

    matplotlibcpp::show();

    // matplotlibcpp::detail::_interpreter::kill();
}
*/

class PlotHelper {

  public:
    /**
     * @brief constructor for PlotHelper class
     * @param mythen_detector shared pointer to MythenDetector holding the
     * detector geometry parameters
     * @param angular_range angular range {min_angle, max_angle} of the
     * diffraction pattern [degrees]
     * @param histogram_bin_width width of the fixed angle width bins [degrees]
     */
    PlotHelper(std::shared_ptr<MythenDetector> mythen_detector,
               const std::tuple<double, double> angular_range,
               const double histogram_bin_width)
        : m_mythen_detector(mythen_detector), m_angular_range(angular_range),
          m_histogram_bin_width(histogram_bin_width) {

        bin_to_diffraction_angle = [this](const size_t bin_index) {
            return bin_index * m_histogram_bin_width +
                   std::get<0>(m_angular_range); // m_min_angle;
        };

        // matplotlibcpp::backend("TkAgg");
        matplotlibcpp::ion();
    };

    ~PlotHelper() {
        // TODO: not good if used in local scope - should be a static function
        // !!
        matplotlibcpp::detail::_interpreter::kill();
    }

    /// @brief pause the program until user input, e.g. before destroying all
    /// plots
    static void pause() {
        std::cout << "Press Enter to continue..." << std::endl;
        std::cin.get();
    }

    /// @brief overwrite plot for plotting multiple frames
    /// in one figure
    void overwrite_plot() {
        m_overwrite_plot = true;
        create_plot();
    }

    /**
     * @brief reuse plot for plotting multiple frames in one figure without
     * overwriting the previous one
     */
    void reuse_plot() {
        m_overwrite_plot = false;
        m_reuse_plot = true;
        create_plot();
    }

    void create_plot() const {
        matplotlibcpp::figure();
        // matplotlibcpp::pause(0.05);
    }

    void initialize_axis(const std::string &plot_title) const {
        matplotlibcpp::xlabel("Diffraction Angle [degree]");
        matplotlibcpp::ylabel("Photon Counts");
        matplotlibcpp::title(plot_title);
    }

    void show_plot() const {
        matplotlibcpp::draw();

        std::atomic<bool> stop_interactive{false};
        std::thread input_thread([&stop_interactive]() {
            std::cout << "Place plot window to intended location and press "
                         "Enter to continue..."
                      << std::endl;
            std::cin.get();
            stop_interactive = true;
        });
        while (!stop_interactive) {
            matplotlibcpp::pause(0.05); // allow rendering
        }
        input_thread.detach();
    }

    /**  @brief plot diffraction pattern for given photon counts
     * @param motor_position The position of the motor for resulting
     * diffraction pattern
     * @param photon_counts resulting photon counts of diffraction pattern
     * @param module_index optional module index for plotting only specific
     * module region
     */
    void plot_diffraction_pattern(
        const NDView<double, 2> photon_counts,
        std::optional<size_t> module_index = std::nullopt,
        std::optional<double> motor_position = std::nullopt) const;

    /**
     * @brief plot base peak for given photon counts
     * @param photon_counts resulting photon counts of diffraction pattern
     * @param base_peak_angle angle of the base peak [degrees]
     * @param base_peak_roi_width width of the region of interest around the
     * base peak [degrees]
     * @param module_index optional module index if only one module contributed
     * to diffraction pattern
     */
    void plot_base_peak(
        const NDView<double, 2> photon_counts, const double base_peak_angle,
        const double base_peak_roi_width,
        const std::optional<size_t> module_index = std::nullopt) const;

  private:
    inline size_t left_module_boundary_as_fixed_angle_bin_index(
        const size_t module_index, const double detector_angle) const {

        return static_cast<size_t>(
            (m_mythen_detector->diffraction_angle_from_DG_parameters(
                 module_index, detector_angle, 0, 0.0) -
             std::get<0>(m_angular_range)) /
            m_histogram_bin_width);
    }

    inline size_t right_module_boundary_as_fixed_angle_bin_index(
        const size_t module_index, const double detector_angle) const {

        // right edge of last strip (strip index denotes left strip boundary)
        return static_cast<size_t>(
            (m_mythen_detector->diffraction_angle_from_DG_parameters(
                 module_index, detector_angle,
                 MythenDetector::strips_per_module() - 1, +1.0) -
             std::get<0>(m_angular_range)) /
            m_histogram_bin_width);
    }

  private:
    /// @brief shared pointer to detector geometry of m_anglecalibration for
    /// calculating diffraction angles
    const std::shared_ptr<MythenDetector> m_mythen_detector{};

    /// @brief angular range of diffraction pattern [degrees]
    const std::tuple<double, double> m_angular_range{};

    /// @brief histogram bin width [degrees]
    const double m_histogram_bin_width{};

    /// @brief function to convert bin index to diffraction angle
    std::function<double(const size_t)> bin_to_diffraction_angle{};

    bool m_overwrite_plot = false;

    bool m_reuse_plot = false;
};

inline void
PlotHelper::plot_base_peak(const NDView<double, 2> photon_counts,
                           const double base_peak_angle,
                           const double base_peak_roi_width,
                           const std::optional<size_t> module_index) const {

    std::string plot_title{};
    if (module_index.has_value()) {
        plot_title =
            fmt::format("Base Peak for module {}", module_index.value());
    } else {
        plot_title = "Base Peak for all modules";
    }

    if (!m_reuse_plot && !m_overwrite_plot) {
        create_plot();
    }

    if (m_overwrite_plot) {
        matplotlibcpp::cla();
    }

    initialize_axis(plot_title);

    /*
    auto bin_to_diffraction_angle_base_peak_ROI_only =
        [this, base_peak_angle, base_peak_roi_width](const ssize_t bin_index) {
            return (bin_index * m_histogram_bin_width - base_peak_roi_width +
                    base_peak_angle);
        };
    */

    ssize_t num_bins_base_peak_ROI =
        static_cast<ssize_t>((base_peak_angle + base_peak_roi_width) /
                             m_histogram_bin_width) -
        static_cast<ssize_t>((base_peak_angle - base_peak_roi_width) /
                             m_histogram_bin_width) +
        1;

    std::vector<double> bins(num_bins_base_peak_ROI);

    std::generate(
        bins.begin(), bins.end(),
        [this, n = size_t{0}, base_peak_angle, base_peak_roi_width]() mutable {
            return (n++) * m_histogram_bin_width - base_peak_roi_width +
                   base_peak_angle;
        });

    size_t left_boundary_bin =
        (base_peak_angle - base_peak_roi_width - std::get<0>(m_angular_range)) /
        m_histogram_bin_width;
    size_t right_boundary_bin =
        (base_peak_angle + base_peak_roi_width - std::get<0>(m_angular_range)) /
            m_histogram_bin_width +
        1;

    std::vector<double> photon_counts_vec(right_boundary_bin -
                                          left_boundary_bin);
    for (size_t i = left_boundary_bin; i < right_boundary_bin; ++i) {
        photon_counts_vec[i - left_boundary_bin] = photon_counts(i, 0);
    }

    matplotlibcpp::plot(bins, photon_counts_vec);

    if (!m_reuse_plot) {
        std::cout << "showing plot for base peak" << std::endl;
        show_plot();
    }
}

// defined in hpp file otherwise random seg faults as singleton _interpreter of
// matplotlibcpp is not properly initialized - not sure why
inline void PlotHelper::plot_diffraction_pattern(
    const NDView<double, 2> photon_counts, std::optional<size_t> module_index,
    std::optional<double> motor_position) const {

    std::string plot_title{};

    if (module_index.has_value() && !motor_position.has_value()) {
        throw std::invalid_argument(
            "If module_index is given, motor_position must also be given for "
            "plotting diffraction pattern of specific module region");
    }

    if (module_index.has_value()) {
        plot_title = fmt::format("Diffraction Pattern for module {}",
                                 module_index.value());
    } else {
        plot_title = "Diffraction Pattern";
    }

    if (!m_overwrite_plot && !m_reuse_plot) {
        create_plot();
    }

    if (m_overwrite_plot) {
        matplotlibcpp::cla();
    }

    initialize_axis(plot_title);

    size_t left_bin_boundary = 0;
    size_t right_bin_boundary = photon_counts.shape(0);

    const ssize_t num_fixed_angle_width_bins =
        std::ceil(std::get<1>(m_angular_range) / m_histogram_bin_width) -
        std::floor(std::get<0>(m_angular_range) / m_histogram_bin_width);

    if (module_index.has_value() &&
        (photon_counts.shape(0) == num_fixed_angle_width_bins)) {
        left_bin_boundary = left_module_boundary_as_fixed_angle_bin_index(
            module_index.value(), motor_position.value());
        right_bin_boundary = right_module_boundary_as_fixed_angle_bin_index(
            module_index.value(), motor_position.value());
    }

    std::vector<double> bins(right_bin_boundary - left_bin_boundary);
    std::generate(bins.begin(), bins.end(),
                  [this, n = left_bin_boundary]() mutable {
                      return bin_to_diffraction_angle(n++);
                  });

    std::vector<double> photon_counts_vec(right_bin_boundary -
                                          left_bin_boundary);

    for (size_t i = left_bin_boundary; i < right_bin_boundary; ++i) {
        photon_counts_vec[i - left_bin_boundary] =
            photon_counts(i, 0); // TODO: better utility for slicing in NDView
    }

    // photon_counts_vec.assign(photon_counts.begin() + left_bin_boundary,
    //  photon_counts.begin() + right_bin_boundary);

    matplotlibcpp::plot(bins, photon_counts_vec);

    if (!m_reuse_plot) {
        show_plot();
    }
}

} // namespace angcal