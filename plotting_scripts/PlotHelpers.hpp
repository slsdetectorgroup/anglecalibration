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

    void plot_base_peak_of_module(
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

} // namespace angcal