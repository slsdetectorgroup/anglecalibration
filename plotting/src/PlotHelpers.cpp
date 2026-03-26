#include "PlotHelpers.hpp"

namespace angcal {

PlotHelper::PlotHelper(std::shared_ptr<AngleCalibration> anglecalibration)
    : m_anglecalibration(anglecalibration) {

    bin_to_diffraction_angle = [this](const size_t bin_index) {
        return bin_index * m_anglecalibration->get_histogram_bin_width() +
               m_anglecalibration->get_angular_range().first; // m_min_angle;
    };

    bin_to_diffraction_angle_base_peak_ROI_only =
        [this](const ssize_t bin_index) {
            return (bin_index * m_anglecalibration->get_histogram_bin_width() -
                    m_anglecalibration->get_base_peak_ROI_width() +
                    m_anglecalibration->get_base_peak_angle());
        };

    matplotlibcpp::backend("TkAgg");
    matplotlibcpp::ion();
}

PlotHelper::~PlotHelper() {
    // before destroying all plots wait
    pause();
}

void PlotHelper::pause() {
    std::cout << "Press Enter to continue..." << std::endl;
    std::cin.get();
}

void PlotHelper::create_plot() const {
    matplotlibcpp::figure();
    matplotlibcpp::clf();
    matplotlibcpp::pause(0.05);
    // matplotlibcpp::figure_size(800, 600);
}

void PlotHelper::overwrite_plot() {
    m_overwrite_plot = true;
    create_plot();
}

void PlotHelper::initialize_axis(const std::string &plot_title) const {
    matplotlibcpp::xlabel("Diffraction Angle [degree]");
    matplotlibcpp::ylabel("Photon Counts");
    matplotlibcpp::title(plot_title);
}

/*
void plot_module_redistributed_to_fixed_angle_width_bins(
    const std::filesystem::path &file_path,
    std::optional<size_t> module_index = std::nullopt) const {
    const double frame_angle =
        m_anglecalibration->get_file_reader()->read_detector_angle(
            file_path);

    NDArray<double, 1> module_redistributed_to_fixed_angle_bins{};

    if (module_index.has_value()) {
        module_redistributed_to_fixed_angle_bins =
            m_anglecalibration->convert({file_path}, module_index.value());
    } else {
        module_redistributed_to_fixed_angle_bins =
            m_anglecalibration->convert({file_path});
    }

    plot_module_redistributed_to_fixed_angle_width_bins(
        frame_angle, module_redistributed_to_fixed_angle_bins.view(),
        module_index);
}
*/

void PlotHelper::plot_diffraction_pattern(
    const double motor_position, const NDView<double, 1> &photon_counts,
    std::optional<size_t> module_index) const {

    std::string plot_title{};
    if (module_index.has_value()) {
        plot_title = fmt::format("Diffraction Pattern for module {}",
                                 module_index.value());
    } else {
        plot_title = "Diffraction Pattern";
    }

    if (!m_overwrite_plot) {
        create_plot();
    } else {
        matplotlibcpp::cla();
    }

    initialize_axis(plot_title);

    size_t left_bin_boundary = 0;
    size_t right_bin_boundary = photon_counts.size();
    if (module_index.has_value()) {
        left_bin_boundary = left_module_boundary_as_fixed_angle_bin_index(
            module_index.value(), motor_position);
        right_bin_boundary = right_module_boundary_as_fixed_angle_bin_index(
            module_index.value(), motor_position);
    }

    std::vector<double> bins(right_bin_boundary - left_bin_boundary);
    std::generate(bins.begin(), bins.end(),
                  [this, n = left_bin_boundary]() mutable {
                      return bin_to_diffraction_angle(n++);
                  });

    std::vector<double> photon_counts_vec;
    photon_counts_vec.assign(photon_counts.begin() + left_bin_boundary,
                             photon_counts.begin() + right_bin_boundary);

    matplotlibcpp::plot(bins, photon_counts_vec);

    matplotlibcpp::draw();
    // matplotlibcpp::pause(2.0);
    matplotlibcpp::pause(0.05);
}

void PlotHelper::plot_base_peak(
    const NDView<double, 1> photon_counts,
    const std::optional<size_t> module_index) const {

    std::string plot_title{};
    if (module_index.has_value()) {
        plot_title =
            fmt::format("Base Peak for module {}", module_index.value());
    } else {
        plot_title = "Base Peak for all modules";
    }

    if (!m_overwrite_plot) {
        create_plot();
    } else {
        matplotlibcpp::cla();
    }

    // matplotlibcpp::figure_size(800, 600);
    initialize_axis(plot_title);

    std::vector<double> bins(m_anglecalibration->get_base_peak_ROI_num_bins());

    std::generate(bins.begin(), bins.end(), [this, n = size_t{0}]() mutable {
        return bin_to_diffraction_angle_base_peak_ROI_only(n++);
    });

    size_t left_boundary_bin = (m_anglecalibration->get_base_peak_angle() -
                                m_anglecalibration->get_base_peak_ROI_width() -
                                m_anglecalibration->get_angular_range().first) /
                               m_anglecalibration->get_histogram_bin_width();
    size_t right_boundary_bin =
        (m_anglecalibration->get_base_peak_angle() +
         m_anglecalibration->get_base_peak_ROI_width() -
         m_anglecalibration->get_angular_range().first) /
            m_anglecalibration->get_histogram_bin_width() +
        1;

    std::vector<double> photon_counts_vec;
    photon_counts_vec.assign(photon_counts.begin() + left_boundary_bin,
                             photon_counts.begin() + right_boundary_bin);

    matplotlibcpp::plot(bins, photon_counts_vec);
    matplotlibcpp::draw();
    matplotlibcpp::pause(0.05);
}

} // namespace angcal