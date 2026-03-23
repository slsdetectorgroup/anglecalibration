#include "PlotHelper.hpp"
#include "AngleCalibration.hpp"

namespace angcal {

/*
PlotHelper::PlotHelper(const std::shared_ptr<AngleCalibration> anglecalibration)
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
}
*/

PlotCalibrationProcess::PlotCalibrationProcess(
    const AngleCalibration *anglecalibration, const std::string &plot_title)
    : m_anglecalibration(anglecalibration) {

    bin_to_diffraction_angle_base_peak_ROI_only =
        [this](const ssize_t bin_index) {
            return (bin_index * m_anglecalibration->get_histogram_bin_width() -
                    m_anglecalibration->get_base_peak_ROI_width() +
                    m_anglecalibration->get_base_peak_angle());
        };

    m_fig = matplot::figure(true);
    m_ax = m_fig->current_axes();

    initializeplot(plot_title);
}

void PlotCalibrationProcess::initializeplot(const std::string &plot_title) {
    m_ax->xlabel("Diffraction Angle [degree]");
    m_ax->ylabel("Photon Counts");
    m_ax->title(plot_title);
    m_ax->hold(matplot::on); // keep previous curves
    m_ax->line_width(2);
    m_fig->size(800, 600);
}

void PlotCalibrationProcess::add_curve(
    const NDView<double, 1> &photon_counts_base_peak_ROI) {
    std::vector<double> bins(photon_counts_base_peak_ROI.size());
    std::generate(bins.begin(), bins.end(), [this, n = size_t{0}]() mutable {
        return bin_to_diffraction_angle_base_peak_ROI_only(n++);
    });

    std::vector<double> counts(photon_counts_base_peak_ROI.size());
    counts.assign(photon_counts_base_peak_ROI.begin(),
                  photon_counts_base_peak_ROI.end()); // uff need std::span
    m_ax->plot(bins, counts); // append line to the existing axes
}

void PlotCalibrationProcess::show() const {
    m_ax->draw();
    std::this_thread::sleep_for(
        std::chrono::seconds(1)); // give the plot some time to render
    m_ax->clear();                // clear the axes for the next plot
}

/*
void PlotHelper::plotCalibrationStep() {

    std::vector<double> bins(redistributed_base_peaks[0].size());
    std::transform(bins.begin(), bins.end(), bins.begin(),
                   [this](size_t i) { return bin_to_diffraction_angle(i); });

    for (const auto &redistributed_base_peak : redistributed_base_peaks) {

        matplot::plot(bins, redistributed_base_peak);
        matplot::hold(on);
    }

    matplot::xlabel("Diffraction Angle [degree]");
    matplot::ylabel("Photon Counts");
    matplot::title("Base Peak Region of Interest");
    matplot::show();

    redistributed_base_peaks.clear();

    redistributed_base_peaks.reserve(
        m_anglecalibration->get_file_list().size() *
        2); // reserve space for all frames and max 2 modules
}
*/

} // namespace angcal