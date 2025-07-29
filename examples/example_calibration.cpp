#include "AngleCalibration.hpp"
#include "logger.hpp"
#include <filesystem>

#ifdef ANGCAL_PLOT
#include "plot_histogram.hpp"
#endif

inline auto data_path() {
    if (const char *env_p = std::getenv("ANGCAL_TEST_DATA")) {
        return std::filesystem::path(env_p);
    } else {
        throw std::runtime_error("Path to test data: $AARE_TEST_DATA not set");
    }
}

using namespace angcal;

int main() {

    auto file_path =
        data_path() / "VariaMay2025/Antonio20250512/angcal_M3_Mar21_2/";

    // auto acquisition_path = data_path() / "VariaMay2025/Antonio20250512";

    assert(std::filesystem::exists(file_path));

    std::shared_ptr<MythenDetectorSpecifications> mythen_detector_ptr =
        std::make_shared<MythenDetectorSpecifications>();

    std::string bad_channels_filename = file_path / "bcX.txt";

    assert(std::filesystem::exists(bad_channels_filename));

    mythen_detector_ptr->read_bad_channels_from_file(bad_channels_filename);

    LOG(TLogLevel::logDEBUG) << "read bad channels";

    std::string unconnected_modules_filename = file_path / "ModOut.txt";

    /*
    auto unconnected_modules =
        read_unconnected_modules(unconnected_modules_filename);

    mythen_detector_ptr->set_unconnected_modules(unconnected_modules);
    */

    std::shared_ptr<FlatField> flat_field_ptr =
        std::make_shared<FlatField>(mythen_detector_ptr);

    std::string flatfield_filename =
        file_path / "Flatfield_E17p5keV_T8751eV_MIX_Mar2021_open_WS.raw";

    flat_field_ptr->read_flatfield_from_file(flatfield_filename);

    LOG(TLogLevel::logDEBUG) << "read flatfield from file";

    AngleCalibration anglecalibration(mythen_detector_ptr, flat_field_ptr,
                                      file_path);

    std::string initial_angles_filename = file_path / "angcal_Mar2021_P10.off";

    anglecalibration.read_initial_calibration_from_file(
        initial_angles_filename);

    LOG(TLogLevel::logDEBUG) << "read initial parameters from file";

    /*
    std::vector<std::string> filelist(20);
    int i = 20;
    std::generate(filelist.begin(), filelist.end(), [&i]() {
        return "ang1up_22keV_LaB60p3mm_48M_a_03" + std::to_string(i++) + ".h5";
    });
    */

    std::vector<std::string> filelist{"ang1dnSi0p3mm_0160.h5"};

    const double base_peak_angle = 36.0568;
    // 35.550; // angpeak=69.225 - maybe select
    // one yourself !!! -could this be off?

    LOG(TLogLevel::logDEBUG) << "starting calibration";
    // anglecalibration.calculate_similarity_of_peaks(0);
    // anglecalibration.calibrate(filelist, base_peak_angle);

    anglecalibration.set_base_peak_angle(base_peak_angle);

    MythenFileReader mythen_file_reader(file_path);

    NDArray<double, 1> S2(
        std::array<ssize_t, 1>{anglecalibration.get_base_peak_ROI_num_bins()},
        0.0);
    NDArray<double, 1> S1(
        std::array<ssize_t, 1>{anglecalibration.get_base_peak_ROI_num_bins()},
        0.0);
    NDArray<double, 1> S0(
        std::array<ssize_t, 1>{anglecalibration.get_base_peak_ROI_num_bins()},
        0.0);

    for (const auto &file : filelist) {
        MythenFrame frame = mythen_file_reader.read_frame(file);

        std::cout << "frame_detector_angle: " << frame.detector_angle
                  << std::endl;

#ifdef ANGCAL_PLOT
        // plot everything
        plot_photon_counts(frame.photon_counts.view(),
                           {0, mythen_detector_ptr->num_strips()}, std::nullopt,
                           *mythen_detector_ptr);

        auto new_fixed_angle_width_bins_photon_counts =
            anglecalibration
                .redistribute_photon_counts_to_fixed_angle_width_bins(frame);

        const double bin_width = anglecalibration.get_histogram_bin_width();

        // plot all
        double most_left_strip_boundary_angle =
            (anglecalibration.diffraction_angle_from_DG_parameters(0, 0, -0.5) +
             frame.detector_angle + 180.0) /
            bin_width;

        double most_right_strip_boundary_angle =
            (anglecalibration.diffraction_angle_from_DG_parameters(
                 mythen_detector_ptr->max_modules() / 2, 1279, 0.5) +
             frame.detector_angle + 180.0) /
            bin_width;

        auto bin_to_diffraction_angle = [&bin_width](const size_t bin_index) {
            return bin_index * bin_width - 180.0;
        };

        plot_photon_counts_for_fixed_angle_width_bins(
            new_fixed_angle_width_bins_photon_counts.view(),
            {most_left_strip_boundary_angle, most_right_strip_boundary_angle},
            bin_to_diffraction_angle);

#endif

        for (size_t module_index = 0;
             module_index < mythen_detector_ptr->max_modules();
             ++module_index) {
            // base peak angle is in module
            if (anglecalibration.base_peak_is_in_module(module_index,
                                                        frame.detector_angle)) {

                LOG(TLogLevel::logINFO)
                    << "base_peak is in module: " << module_index << "\n";
                /*
                auto [fixed_angle_width_bins_photon_counts,
                      fixed_angle_width_bins_photon_count_variance] =
                    anglecalibration
                        .redistribute_photon_counts_to_fixed_angle_bins(
                            module_index, frame, S0.view(), S1.view(),
                            S2.view());
                */

#ifdef ANGCAL_PLOT

                // plot module
                double left_strip_boundary_angle =
                    (anglecalibration.diffraction_angle_from_DG_parameters(
                         module_index, 0, -0.5) +
                     frame.detector_angle + 180.0) /
                    anglecalibration.get_histogram_bin_width();

                double right_strip_boundary_angle =
                    (anglecalibration.diffraction_angle_from_DG_parameters(
                         module_index, 1279, 0.5) +
                     frame.detector_angle + 180.0) /
                    anglecalibration.get_histogram_bin_width();

                plot_photon_counts_for_fixed_angle_width_bins(
                    new_fixed_angle_width_bins_photon_counts.view(),
                    {left_strip_boundary_angle, right_strip_boundary_angle},
                    bin_to_diffraction_angle, module_index);

                // plot roi_of_peak
                NDArray<double, 1> base_peak_roi_photon_counts(
                    std::array<ssize_t, 1>{
                        anglecalibration.get_base_peak_ROI_num_bins()},
                    0.0);
                NDArray<double, 1> base_peak_roi_photon_counts_variance(
                    std::array<ssize_t, 1>{
                        anglecalibration.get_base_peak_ROI_num_bins()},
                    0.0);

                const size_t base_peak_roi =
                    anglecalibration.get_base_peak_ROI();
                auto bin_to_diffraction_angle_base_peak_ROI_only =
                    [&bin_width, &base_peak_angle,
                     &base_peak_roi](const size_t bin_index) {
                        return (static_cast<ssize_t>(bin_index) -
                                static_cast<ssize_t>(base_peak_roi)) *
                                   bin_width +
                               base_peak_angle;
                    };

                anglecalibration
                    .redistribute_photon_counts_to_fixed_angle_width_bins<true>(
                        module_index, frame, base_peak_roi_photon_counts.view(),
                        base_peak_roi_photon_counts_variance.view());

                plot_photon_counts_for_fixed_angle_width_bins(
                    base_peak_roi_photon_counts.view(), {0, 101},
                    bin_to_diffraction_angle_base_peak_ROI_only, module_index);

                /*
                plot_photon_counts(
                    frame.photon_counts.view(),
                    {module_index * mythen_detector_ptr->strips_per_module(),
                     (module_index + 1) *
                         mythen_detector_ptr->strips_per_module()},
                    module_index, *mythen_detector_ptr);

                size_t left_module_boundary_angle = static_cast<size_t>(
                    (anglecalibration.diffraction_angle_from_DG_parameters(
                         module_index, 0, 0.5) +
                     frame.detector_angle) /
                    anglecalibration.get_histogram_bin_width());
                size_t right_module_boundary_angle = static_cast<size_t>(
                    (anglecalibration.diffraction_angle_from_DG_parameters(
                         module_index, mythen_detector_ptr->strips_per_module(),
                         -0.5) +
                     frame.detector_angle) /
                    anglecalibration.get_histogram_bin_width());

                plot_photon_counts_for_fixed_angle_width_bins(
                    fixed_angle_width_bins_photon_counts.view(),
                    {left_module_boundary_angle, right_module_boundary_angle},
                    module_index, anglecalibration.get_histogram_bin_width());
                */
#endif
            }
        }
    }
}
