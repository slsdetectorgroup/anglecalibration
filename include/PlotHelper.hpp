#pragma once
#include "aare/NDView.hpp"
// #include "matplot/matplot.h"
#include "matplotlibcpp.h"

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

    inline void initializeplot();

    void add_curve(const aare::NDView<double, 1> &photon_counts_base_peak_ROI);

    void close() { matplotlibcpp::close(); }

    void show();

  private:
    const AngleCalibration *m_anglecalibration;
    std::function<double(const size_t)>
        bin_to_diffraction_angle_base_peak_ROI_only;
    const std::string m_plot_title{};
    std::vector<double> bins{};
};

} // namespace angcal