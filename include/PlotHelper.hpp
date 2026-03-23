#pragma once
#include "aare/NDView.hpp"
#include "matplot/matplot.h"

namespace angcal {

class AngleCalibration; // forward declaration to avoid circular dependency

/*
class PlotHelper {

  public:
    PlotHelper(const std::shared_ptr<AngleCalibration> anglecalibration);

    ~PlotHelper() = default;

    void plotCalibrationStep();

    void plotFrameRedistributedToFixedAngleWidthBins(
        const size_t module_index,
        const NDView<double, 1> redistributed_photon_counts,
        const double detector_angle) const;

    void plotOptimizationError();

  private:
    std::shared_ptr<AngleCalibration> m_anglecalibration;

    std::function<double(const size_t)> bin_to_diffraction_angle;
    std::function<double(const size_t)>
        bin_to_diffraction_angle_base_peak_ROI_only;

    std::vector<NDArray<double, 1>> redistributed_base_peaks;
};
*/

class PlotCalibrationProcess {

  public:
    PlotCalibrationProcess(const AngleCalibration *anglecalibration,
                           const std::string &plot_title);

    ~PlotCalibrationProcess() = default;

    void initializeplot(const std::string &plot_title);

    void add_curve(const aare::NDView<double, 1> &photon_counts_base_peak_ROI);

    void show() const;

  private:
    const AngleCalibration *m_anglecalibration;
    std::function<double(const size_t)>
        bin_to_diffraction_angle_base_peak_ROI_only;
    matplot::figure_handle m_fig;
    matplot::axes_handle m_ax;
};

} // namespace angcal