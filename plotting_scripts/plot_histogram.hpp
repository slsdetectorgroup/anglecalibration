#pragma once

#include "MythenDetectorSpecifications.hpp"
#include "MythenFileReader.hpp"
#include "aare/NDView.hpp"
#include <functional>
#include <gnuplot-iostream.h>
#include <map>

// size_t window_id = 0;

namespace angcal {

inline void
plot_photon_counts(const aare::NDView<uint32_t, 1> photon_counts,
                   const std::pair<size_t, size_t> bin_range,
                   std::optional<const size_t> module_index,
                   std::optional<std::shared_ptr<MythenDetectorSpecifications>>
                       mythen_detector_specifications = std::nullopt) {

    std::vector<std::pair<int, int>> plot_data;
    const size_t num_bins = bin_range.second - bin_range.first;
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
    // gp << "set terminal qt" << window_id << "\n";
    if (module_index.has_value()) {
        gp << "set title 'Photon Counts for module " << module_index.value()
           << "'\n";
    } else {
        gp << "set title 'Photon Counts for all modules'\n";
    }
    gp << "set xlabel 'Strip Index'\n";
    gp << "set ylabel 'Photon Count'\n";
    gp << "plot '-' using 1:2 with lines lc rgb 'orange' title 'Photon "
          "Counts'\n";
    // gp << "plot '-' using 1:2 with points pointtype 7 pointsize 1 lc rgb "
    //"'orange' title 'Photon Counts'\n";

    gp.send1d(plot_data);
    gp << "bind all 'alt-End' 'exit gnuplot'\n";
    // gp << "pause mouse close\n";
    gp << "pause -1 'Close window or press Enter to continue'\n";

    std::cout << "Press Enter to continue...\n";
    std::cin.get();

    //++window_id;
}

// maybe add a range
// TODO add st::enable_if for FUnc
template <typename Func>
inline void plot_photon_counts_for_fixed_angle_width_bins(
    const aare::NDView<double, 1> photon_counts,
    const std::pair<size_t, size_t> bin_range,
    Func &&bin_index_to_angle_conversion,
    std::optional<const size_t> module_index = std::nullopt) {

    std::vector<std::pair<double, double>> plot_data;
    const size_t num_bins = bin_range.second - bin_range.first;
    plot_data.reserve(num_bins);

    for (size_t bin = bin_range.first; bin < bin_range.second; ++bin) {
        plot_data.emplace_back(bin_index_to_angle_conversion(bin),
                               photon_counts(bin));
    }

    Gnuplot gp;
    gp << "set terminal qt persist\n";
    if (module_index.has_value()) {
        gp << "set title 'Photon Counts for module " << module_index.value()
           << "'\n";
    } else {
        gp << "set title 'Photon Counts for all modules'\n";
    }
    gp << "set xlabel 'Diffraction Angle'\n";
    gp << "set ylabel 'Photon Count'\n";
    gp << "plot '-' using 1:2 with lines lc rgb 'orange' title 'Photon "
          "Counts'\n";
    // gp << "plot '-' using 1:2 with points pointtype 7 pointsize 1 lc rgb "
    //"'orange' title 'Photon Counts'\n";

    gp.send1d(plot_data);
    gp << "bind all 'alt-End' 'exit gnuplot'\n";
    // gp << "pause mouse close\n";
    gp << "pause -1 'Close window or press Enter to continue'\n";

    std::cout << "Press Enter to continue...\n";
    std::cin.get();
}

} // namespace angcal