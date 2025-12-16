#pragma once
#include "AngleCalibration.hpp"
#include "logger.hpp"

#ifdef ANGCAL_PLOT
#include "PlotHelpers.hpp"
#endif

using namespace angcal;

std::pair<double, double>
get_detector_range(std::shared_ptr<MythenFileReader> mythen_file_reader,
                   const std::vector<std::string> &file_list) {

    auto frame = mythen_file_reader->read_frame(file_list[0]);

    double min_detector_angle =
        mythen_file_reader->read_frame(file_list[0]).detector_angle;
    double max_detector_angle =
        mythen_file_reader->read_frame(file_list[file_list.size() - 1])
            .detector_angle; // assuming frames are in order

    // Files are ordered by angle already
    /*
    std::vector<double> detector_angles(file_list.size());
    int n = 0;
    std::generate(detector_angles.begin(), detector_angles.end(),
                  [&mythen_file_reader, &file_list, &n]() mutable {
                      return mythen_file_reader->read_frame(file_list[n++])
                          .detector_angle;
                  });

    auto [min, max] =
        std::minmax_element(detector_angles.begin(), detector_angles.end());
    */

    return std::make_pair(max_detector_angle, min_detector_angle);
}

void get_angle_range_of_module(
    std::shared_ptr<MythenFileReader> mythen_file_reader,
    const std::vector<std::string> &file_list,
    AngleCalibration &anglecalibration, const size_t module_index) {

    double detector_angle =
        mythen_file_reader->read_frame(file_list[0]).detector_angle;
    double max_left_module_angle =
        anglecalibration.diffraction_angle_from_DG_parameters(
            module_index, detector_angle, 0, -0.5);
    double max_right_module_angle =
        anglecalibration.diffraction_angle_from_DG_parameters(
            module_index, detector_angle, 1279, +0.5);

    double min_left_module_angle = max_left_module_angle;
    double min_right_module_angle = max_right_module_angle;

    for (auto &frame : file_list) {
        detector_angle = mythen_file_reader->read_frame(frame).detector_angle;
        double left_module_angle =
            anglecalibration.diffraction_angle_from_DG_parameters(
                module_index, detector_angle, 0, -0.5);
        double right_module_angle =
            anglecalibration.diffraction_angle_from_DG_parameters(
                module_index, detector_angle, 1279, +0.5);

        if (left_module_angle < min_left_module_angle &&
            right_module_angle < min_right_module_angle) {
            min_left_module_angle = left_module_angle;
            min_right_module_angle = right_module_angle;
        }
        if (left_module_angle > max_left_module_angle &&
            right_module_angle > max_right_module_angle) {
            max_left_module_angle = left_module_angle;
            max_right_module_angle = right_module_angle;
        }
    }

    LOG(TLogLevel::logINFO) << fmt::format(
        "minimum module ROI [degrees]: {},{}, maximum module ROI [degrees]: "
        "{},{} for module "
        "{}",
        min_left_module_angle, min_right_module_angle, max_left_module_angle,
        max_right_module_angle, module_index);
}

void get_module_angle_ranges(
    std::shared_ptr<MythenFileReader> mythen_file_reader,
    const std::vector<std::string> &file_list,
    AngleCalibration &anglecalibration, const ssize_t num_modules) {

    for (ssize_t module_idx = 0; module_idx < num_modules; ++module_idx) {
        get_angle_range_of_module(mythen_file_reader, file_list,
                                  anglecalibration, module_idx);
    }
}

#ifdef ANGCAL_PLOT
/**
 * @brief select base peak by plotting module 0 redistributed to fixed
 * angle-width bins for each frame in filelist with a detector angle between
 * 6° and 33°
 */
void select_base_peak(std::shared_ptr<AngleCalibration> anglecalibration,
                      std::shared_ptr<MythenFileReader> mythen_file_reader,
                      const std::vector<std::string> &filelist,
                      const size_t module_index = 0) {
    PlotHelper plotter(anglecalibration);

    std::shared_ptr<Gnuplot> gp = std::make_shared<Gnuplot>();

    auto detector_angle_range = std::make_pair(6.0, 33.0);
    double strip_angle = anglecalibration->diffraction_angle_from_DG_parameters(
        module_index, 0.0, 0, -0.5);
    detector_angle_range.first -= strip_angle - 5.0;
    detector_angle_range.second -= strip_angle + 5.0;

    LOG(TLogLevel::logINFO) << fmt::format(
        "adjusted detector angle range for module {}: {} to {}", module_index,
        detector_angle_range.first, detector_angle_range.second);

    for (const auto &file : filelist) {

        LOG(TLogLevel::logDEBUG) << fmt::format("reading file {}", file);
        auto frame = mythen_file_reader->read_frame(file);

        if (frame.detector_angle > detector_angle_range.first &&
            frame.detector_angle < detector_angle_range.second) {

            LOG(TLogLevel::logINFO)
                << fmt::format("plotting file {} with detector angle {}", file,
                               frame.detector_angle);

            // plot module
            auto module_redistributed_to_fixed_angle_bins =
                anglecalibration
                    ->redistribute_photon_counts_to_fixed_angle_width_bins(
                        frame, module_index);

            plotter.plot_module_redistributed_to_fixed_angle_width_bins(
                module_index, module_redistributed_to_fixed_angle_bins.view(),
                frame.detector_angle, gp);

            plotter.pause();
        }
    }
}
#endif
