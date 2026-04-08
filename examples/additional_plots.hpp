#include "AngleCalibration.hpp"
#include "PlotHelpers.hpp"

namespace angcal {

inline void
plot_all_base_peaks(const std::shared_ptr<AngleCalibration> anglecalibration,
                    const std::shared_ptr<MythenFileReader> mythen_file_reader,
                    const std::vector<std::string> &file_list,
                    std::optional<size_t> module_index = std::nullopt) {

    PlotHelper plotter(anglecalibration);

    plotter.reuse_plot();

    std::vector<size_t> modules_to_plot{};

    if (module_index.has_value()) {
        modules_to_plot.push_back(module_index.value());
    } else {
        modules_to_plot.resize(
            anglecalibration->get_detector_specifications()->max_modules);
        std::iota(modules_to_plot.begin(), modules_to_plot.end(), 0);
    }

    for (const auto &file : file_list) {
        const double detector_angle =
            mythen_file_reader->read_detector_angle(file);

        for (const size_t module_idx : modules_to_plot) {
            if (anglecalibration->base_peak_is_in_module(module_idx,
                                                         detector_angle)) {
                auto diffraction_pattern =
                    anglecalibration->convert({file}, module_idx);

                plotter.plot_base_peak(diffraction_pattern.view());
            }
        }
    }

    plotter.show_plot();
}

} // namespace angcal