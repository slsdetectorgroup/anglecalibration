#pragma once
#include "aare/NDView.hpp"
// #include "matplot/matplot.h"
#include "matplotlibcpp.h"

namespace angcal {

class AngleCalibration; // forward declaration to avoid circular dependency

class PlotCalibrationProcess {

  public:
    PlotCalibrationProcess(const AngleCalibration *anglecalibration,
                           const std::string &plot_title);

    ~PlotCalibrationProcess() = default;

    static void kill_python_interpreter() {
        // matplotlibcpp::detail::_interpreter::kill();
    }

    inline void initializeplot() {
        matplotlibcpp::xlabel("Diffraction Angle [degree]");
        matplotlibcpp::ylabel("Photon Counts");
        matplotlibcpp::title(m_plot_title);
    }

    inline void initializematplotlib() {
        if (!Py_IsInitialized()) {
            Py_Initialize();
        }
        // matplotlibcpp::backend("TkAgg");
        matplotlibcpp::figure_size(800, 600);
        matplotlibcpp::ion(); // interactive mode on to update the plot after
                              // each curve is added
        initializeplot();
    }

    void add_curve(const aare::NDView<double, 1> &photon_counts_base_peak_ROI) {
        std::vector<double> counts(photon_counts_base_peak_ROI.size());
        counts.assign(photon_counts_base_peak_ROI.begin(),
                      photon_counts_base_peak_ROI.end()); // uff need std::span

        matplotlibcpp::plot(bins,
                            counts); // append line to the existing figure
    }

    void close() { matplotlibcpp::close(); }

    void show() {
        matplotlibcpp::draw();
        matplotlibcpp::pause(0.5); // allow rendering
        matplotlibcpp::cla();      // clear the figure for the next plot
        initializeplot();
    }

  private:
    const AngleCalibration *m_anglecalibration;
    std::function<double(const size_t)>
        bin_to_diffraction_angle_base_peak_ROI_only;
    const std::string m_plot_title{};
    std::vector<double> bins{};
};

} // namespace angcal