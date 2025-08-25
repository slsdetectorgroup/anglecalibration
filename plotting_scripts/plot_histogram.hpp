#pragma once

#include "MythenDetectorSpecifications.hpp"
#include "MythenFileReader.hpp"
#include "aare/NDView.hpp"
#include <functional>
#include <gnuplot-iostream.h>
#include <map>

// size_t window_id = 0;

namespace angcal {

inline size_t plot_window = 0;

// TODO should be more generic
template <typename T>
void plot_photon_counts(
    const aare::NDView<T, 1> photon_counts,
    const std::pair<size_t, size_t> bin_range, const std::string &plot_title,
    std::optional<std::shared_ptr<MythenDetectorSpecifications>>
        mythen_detector_specifications = std::nullopt) {

    // TODO: change to size_t
    std::vector<std::pair<int, T>> plot_data;
    const size_t num_bins = bin_range.second - bin_range.first;

    LOG(TLogLevel::logDEBUG) << "plotting " << num_bins << " bins from "
                             << bin_range.first << " to " << bin_range.second;

    plot_data.reserve(num_bins);

    for (size_t bin = bin_range.first; bin < bin_range.second; ++bin) {
        if (!mythen_detector_specifications.has_value()) {
            plot_data.emplace_back(bin, photon_counts(bin));

        } else if (mythen_detector_specifications.has_value()) {

            if (!mythen_detector_specifications.value()->get_bad_channels()(
                    bin)) {
                plot_data.emplace_back(bin, photon_counts(bin));
            } else {
                plot_data.emplace_back(bin, 0);
            }
        }
    }

    Gnuplot gp;

    gp << "set terminal qt persist\n";

    gp << "set title \" " << plot_title << " \" \n";

    gp << "set xlabel 'Strip Index'\n";
    gp << "set ylabel 'Photon Count'\n";
    gp << "plot '-' using 1:2 with lines lc rgb 'orange' notitle \n";

    gp.send1d(plot_data);

    std::cout << "Press Enter to continue...\n";
    std::cin.get();
}

// maybe add a range
// TODO add st::enable_if for FUnc
template <typename Func>
inline void plot_photon_counts_for_fixed_angle_width_bins(
    Gnuplot &gp, const aare::NDView<double, 1> photon_counts,
    const std::pair<size_t, size_t> bin_range, const std::string &plot_title,
    Func &&bin_index_to_angle_conversion) {

    gp << "set title \"" << plot_title << "\"\n";

    std::vector<std::pair<double, double>> plot_data;
    const size_t num_bins = bin_range.second - bin_range.first;
    plot_data.reserve(num_bins);

    for (size_t bin = bin_range.first; bin < bin_range.second; ++bin) {
        plot_data.emplace_back(bin_index_to_angle_conversion(bin),
                               photon_counts(bin));
    }

    gp << "plot '-' using 1:2 with lines lc rgb 'orange' notitle\n";

    gp.send1d(plot_data);

    gp.flush();
}

class PlotCalibrationProcess {

  public:
    PlotCalibrationProcess(
        const std::string &plot_title,
        const std::string &xlabel = "Diffraction Angle [degree]",
        const std::string &ylabel = "Photon Counts") {
        gp << "set terminal qt " << plot_window++ << " persist\n";

        gp << "set title \"" << plot_title << "\"\n";

        gp << "set xlabel \"" << xlabel << "\"\n";
        gp << "set ylabel \"" << ylabel << "\"\n";

        gp << "plot sin(x) with lines lc rgb 'white' notitle \n"; // dummy
        // gp << "plot NaN title ''\n";

        // gp << "set multiplot \n";
    }

    ~PlotCalibrationProcess() {
        gp << "quit\n";
        gp.flush();
    }

    inline void clear() {
        gp << "clear\n";
        gp.flush();
        gp << "plot sin(x) with lines lc rgb 'white' lw 0 notitle \n"; // dummy
                                                                       // value
    }

    void static pause() {
        std::cout << "Press Enter to continue...\n";
        std::cin.get();
    }

    inline void flush() { gp.flush(); }

    template <typename Func>
    inline void append_to_plot(const aare::NDView<double, 1> photon_counts,
                               const std::pair<size_t, size_t> bin_range,
                               Func &&bin_index_to_angle_conversion,
                               const std::string &datasetname) {

        std::ofstream out(datasetname);
        for (size_t bin = bin_range.first; bin < bin_range.second; ++bin) {
            out << bin_index_to_angle_conversion(bin) << " "
                << photon_counts(bin) << "\n";
        }

        gp << "replot '" << datasetname << "' using 1:2 with lines notitle\n";
        gp.flush();
    }

  private:
    Gnuplot gp;
    static size_t plot_window;
};

inline size_t PlotCalibrationProcess::plot_window = 0;

} // namespace angcal