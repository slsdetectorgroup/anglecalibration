#include "PlotCalibrationProcess.hpp"
#include "AngleCalibration.hpp"

namespace angcal {

PlotCalibrationProcess::PlotCalibrationProcess(
    const AngleCalibration *anglecalibration, const std::string &plot_title)
    : m_anglecalibration(anglecalibration), m_plot_title(plot_title) {

    bin_to_diffraction_angle_base_peak_ROI_only =
        [this](const ssize_t bin_index) {
            return (bin_index * m_anglecalibration->get_histogram_bin_width() -
                    m_anglecalibration->get_base_peak_ROI_width() +
                    m_anglecalibration->get_base_peak_angle());
        };

    bins.resize(m_anglecalibration->get_base_peak_ROI_num_bins());
    std::generate(bins.begin(), bins.end(), [this, n = size_t{0}]() mutable {
        return bin_to_diffraction_angle_base_peak_ROI_only(n++);
    });
}

} // namespace angcal