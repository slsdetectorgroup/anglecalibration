#include "AngleCalibration.hpp"
#include "plot_histogram.hpp"

namespace angcal {
class PlotHelper {

  public:
    PlotHelper(const std::shared_ptr<AngleCalibration> anglecalibration)
        : m_anglecalibration(anglecalibration) {

        bin_to_diffraction_angle = [this](const size_t bin_index) {
            return bin_index * m_anglecalibration->get_histogram_bin_width() +
                   m_anglecalibration->get_detector_specifications()
                       ->min_angle();
        };

        bin_to_diffraction_angle_base_peak_ROI_only =
            [this](const size_t bin_index) {
                return (static_cast<ssize_t>(bin_index) -
                        static_cast<ssize_t>(
                            m_anglecalibration->get_base_peak_ROI())) *
                           m_anglecalibration->get_histogram_bin_width() +
                       m_anglecalibration->get_base_peak_angle();
            };
    }

    void plot_module_redistributed_to_fixed_angle_width_bins(
        const size_t module_index,
        const NDView<double, 1> redistributed_photon_counts,
        const double detector_angle,
        std::shared_ptr<Gnuplot> gp = nullptr) const {

        std::string plot_title =
            fmt::format("Fixed angle bin-width photon count "
                        "histogram for module {}",
                        module_index);

        if (gp == nullptr) {
            gp = std::make_shared<Gnuplot>();
            *gp << "set terminal qt persist\n";
        }

        *gp << "set xlabel \"Diffraction Angle [degree]\"\n";
        *gp << "set ylabel \"Photon Counts\"\n";

        plot_photon_counts_for_fixed_angle_width_bins(
            *gp, redistributed_photon_counts,
            {left_module_boundary_as_fixed_angle_bin_index(module_index,
                                                           detector_angle),
             right_module_boundary_as_fixed_angle_bin_index(module_index,
                                                            detector_angle)},
            plot_title, bin_to_diffraction_angle);
    }

    void plot_base_peak_region_of_interest(
        const size_t module_index,
        const NDView<double, 1> base_peak_ROI_photon_counts,
        std::shared_ptr<Gnuplot> gp = nullptr) const {

        std::string plot_title =
            fmt::format("Base Peak for module {}", module_index);

        if (gp == nullptr) {
            gp = std::make_shared<Gnuplot>();
            *gp << "set terminal qt persist\n";
        }

        *gp << "set xlabel \"Diffraction Angle [degree]\"\n";
        *gp << "set ylabel \"Photon Counts\"\n";
        plot_photon_counts_for_fixed_angle_width_bins(
            *gp, base_peak_ROI_photon_counts,
            {0, m_anglecalibration->get_base_peak_ROI_num_bins()}, plot_title,
            bin_to_diffraction_angle_base_peak_ROI_only);
    }

    void plot_redistributed_photon_counts(
        const NDView<double, 1> redistributed_photon_counts,
        std::shared_ptr<Gnuplot> gp = nullptr) {
        std::string plot_title =
            "Redistributed Photon Counts to fixed angle with bins";

        if (gp == nullptr) {
            gp = std::make_shared<Gnuplot>();
            *gp << "set terminal qt persist\n";
        }

        *gp << "set xlabel \"Diffraction Angle [degree]\"\n";
        *gp << "set ylabel \"Photon Counts\"\n";

        plot_photon_counts_for_fixed_angle_width_bins(
            *gp, redistributed_photon_counts,
            {0, m_anglecalibration->new_number_of_bins()}, plot_title,
            bin_to_diffraction_angle);
    }

    void static pause() {
        std::cout << "Press Enter to continue...\n";
        std::cin.get();
    }

  private:
    inline size_t left_module_boundary_as_fixed_angle_bin_index(
        const size_t module_index, const double detector_angle) const {

        return (m_anglecalibration->diffraction_angle_from_DG_parameters(
                    module_index, detector_angle, 0, -0.5) -
                m_anglecalibration->get_detector_specifications()
                    ->min_angle()) /
               m_anglecalibration->get_histogram_bin_width();
    }

    inline size_t right_module_boundary_as_fixed_angle_bin_index(
        const size_t module_index, const double detector_angle) const {

        return (m_anglecalibration->diffraction_angle_from_DG_parameters(
                    module_index, detector_angle,
                    m_anglecalibration->get_detector_specifications()
                        ->strips_per_module(),
                    +0.5) -
                m_anglecalibration->get_detector_specifications()
                    ->min_angle()) /
               m_anglecalibration->get_histogram_bin_width();
    }

  private:
    const std::shared_ptr<AngleCalibration> m_anglecalibration{};

    std::function<double(const size_t)> bin_to_diffraction_angle{};

    std::function<double(const size_t)>
        bin_to_diffraction_angle_base_peak_ROI_only{};
};

/**
 * @brief plots the module for all frames
 * @tparam true only plot frames where base_peak ROI overlaps with module
 * boundary
 * @param base_peak_boundary [degrees]: artifically enlarge base_peak ROI (for
 * debugging) [base_peak-base_peak_boundary, base_peak+base_peak_boundary]
 */
template <bool only_frames_in_basepeakregion = false>
void plot_module_in_angle_for_all_frames(
    const AngleCalibration &anglecalibration, const PlotHelper &plotter,
    MythenFileReader &mythenfilereader, std::vector<std::string> &filelist,
    const size_t module_index,
    std::optional<double> base_peak_boundary = std::nullopt) {

    std::shared_ptr<Gnuplot> gp = std::make_shared<Gnuplot>();
    size_t frame_number = 0;
    const ssize_t new_num_bins = anglecalibration.new_number_of_bins();

    for (const auto &file : filelist) {

        MythenFrame frame = mythenfilereader.read_frame(file);

        if constexpr (only_frames_in_basepeakregion) {
            if (anglecalibration.base_peak_is_in_module(
                    module_index, frame.detector_angle, base_peak_boundary)) {

                NDArray<double, 1> fixed_angle_width_bins_photon_counts(
                    std::array<ssize_t, 1>{new_num_bins}, 0.0);
                NDArray<double, 1>
                    fixed_angle_width_bins_photon_counts_variance(
                        std::array<ssize_t, 1>{new_num_bins}, 0.0);

                anglecalibration
                    .redistribute_photon_counts_to_fixed_angle_width_bins<
                        false>(
                        module_index, frame,
                        fixed_angle_width_bins_photon_counts.view(),
                        fixed_angle_width_bins_photon_counts_variance.view());

                plotter.plot_module_redistributed_to_fixed_angle_width_bins(
                    module_index, fixed_angle_width_bins_photon_counts.view(),
                    frame.detector_angle, gp);

                std::this_thread::sleep_for(std::chrono::milliseconds(500));
                *gp << "clear \n"; // plot in the same window
            }
        } else {

            NDArray<double, 1> fixed_angle_width_bins_photon_counts(
                std::array<ssize_t, 1>{new_num_bins}, 0.0);
            NDArray<double, 1> fixed_angle_width_bins_photon_counts_variance(
                std::array<ssize_t, 1>{new_num_bins}, 0.0);

            anglecalibration
                .redistribute_photon_counts_to_fixed_angle_width_bins<false>(
                    module_index, frame,
                    fixed_angle_width_bins_photon_counts.view(),
                    fixed_angle_width_bins_photon_counts_variance.view());

            plotter.plot_module_redistributed_to_fixed_angle_width_bins(
                module_index, fixed_angle_width_bins_photon_counts.view(),
                frame.detector_angle, gp);

            plotter.plot_module_redistributed_to_fixed_angle_width_bins(
                module_index, fixed_angle_width_bins_photon_counts.view(),
                frame.detector_angle, gp);

            std::this_thread::sleep_for(std::chrono::milliseconds(500));
            *gp << "clear \n"; // plot in the same window
        }

        LOG(TLogLevel::logDEBUG1) << fmt::format("frame {} of {} frames",
                                                 frame_number, filelist.size());
        ++frame_number;
    }
    plotter.pause();
}

} // namespace angcal