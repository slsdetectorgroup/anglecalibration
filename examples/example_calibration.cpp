#include "AngleCalibration.hpp"
#include "aare/File.hpp"
#include "logger.hpp"
#include <filesystem>

#ifdef ANGCAL_PLOT
#include "PlotHelpers.hpp"
#endif

/**
 * How to run: All data for this example can be found in the github directory
 * https://gitea.psi.ch/angcal/VariaMay2025 before running execute: export
 * ANGCAL_TEST_PATH=Path_to_git_repository
 *
 */
inline auto data_path() {
    if (const char *env_p = std::getenv("ANGCAL_TEST_DATA")) {
        return std::filesystem::path(env_p);
    } else {
        throw std::runtime_error(
            "Path to test data: $ANGCAL_TEST_DATA not set");
    }
}

using namespace angcal;

// helper functions
NDArray<double, 1>
get_redistributed_photon_counts(const AngleCalibration &anglecalibration,
                                const size_t module_index,
                                const MythenFrame &frame) {

    const ssize_t new_num_bins = anglecalibration.new_number_of_bins();

    NDArray<double, 1> fixed_angle_width_bins_photon_counts(
        std::array<ssize_t, 1>{new_num_bins}, 0.0);
    NDArray<double, 1> fixed_angle_width_bins_photon_counts_variance(
        std::array<ssize_t, 1>{new_num_bins}, 0.0);

    anglecalibration
        .redistribute_photon_counts_to_fixed_angle_width_bins<false>(
            module_index, frame, fixed_angle_width_bins_photon_counts.view(),
            fixed_angle_width_bins_photon_counts_variance.view());

    return fixed_angle_width_bins_photon_counts;
}

NDArray<double, 1>
get_base_peak_region_of_interest(const AngleCalibration &anglecalibration,
                                 const size_t module_index,
                                 const MythenFrame &frame) {

    NDArray<double, 1> base_peak_roi_photon_counts(
        std::array<ssize_t, 1>{anglecalibration.get_base_peak_ROI_num_bins()},
        0.0);
    NDArray<double, 1> base_peak_roi_photon_counts_variance(
        std::array<ssize_t, 1>{anglecalibration.get_base_peak_ROI_num_bins()},
        0.0);

    anglecalibration.redistribute_photon_counts_to_fixed_angle_width_bins<true>(
        module_index, frame, base_peak_roi_photon_counts.view(),
        base_peak_roi_photon_counts_variance.view());

    return base_peak_roi_photon_counts;
}

int main() {

    auto file_path = data_path() / "Antonio20250512" / "angcal_M3_Mar21_2";

    assert(std::filesystem::exists(file_path));

    std::shared_ptr<MythenDetectorSpecifications> mythen_detector_ptr =
        std::make_shared<MythenDetectorSpecifications>();

    std::string bad_channels_filename = file_path / "bcX.txt";

    assert(std::filesystem::exists(bad_channels_filename));

    auto bad_channels = mythen_detector_ptr->get_bad_channels();

    mythen_detector_ptr->read_bad_channels_from_file(bad_channels_filename);

    LOG(TLogLevel::logDEBUG) << "read bad channels";

    /*
    std::string unconnected_modules_filename = file_path / "ModOut.txt";

    auto unconnected_modules =
        read_unconnected_modules(unconnected_modules_filename);

    mythen_detector_ptr->set_unconnected_modules(unconnected_modules);
    */

    std::shared_ptr<FlatField> flat_field_ptr =
        std::make_shared<FlatField>(mythen_detector_ptr);

    std::string flatfield_filename =
        file_path / "Flatfield_E17p5keV_T8751eV_MIX_Mar2021_open_WS.raw";

    // TODO: do I want to pass an inverse flatfield or only read flatfield and
    // calculate in Flatfield class
    //  inverse_flatfield is only properly allocated when calling
    //  inverse_flatfield
    NDArray<double, 1> inverse_normalized_flatfield(
        std::array<ssize_t, 1>{mythen_detector_ptr->num_strips()});
    CustomFlatFieldFile flatfieldfilereader;
    flatfieldfilereader.open(flatfield_filename);
    flatfieldfilereader.read_into(inverse_normalized_flatfield.buffer(), 8);
    flat_field_ptr->set_inverse_normalized_flatfield(
        inverse_normalized_flatfield);

#ifdef ANGCAL_PLOT
    plot_photon_counts(inverse_normalized_flatfield.view(),
                       {0, mythen_detector_ptr->num_strips()},
                       "Inverse normalized flatfield", mythen_detector_ptr);
    PlotHelper::pause();
#endif

    LOG(TLogLevel::logDEBUG) << "read flatfield from file";

    // TODO change MythenFileReader to not take path
    MythenFileReader mythen_file_reader(file_path);

    AngleCalibration anglecalibration(mythen_detector_ptr, flat_field_ptr,
                                      file_path);

    std::string initial_angles_filename = file_path / "angcal_Mar2021_P10.off";

    anglecalibration.read_initial_calibration_from_file(
        initial_angles_filename);

    LOG(TLogLevel::logDEBUG) << "read initial parameters from file";

    std::vector<std::string> filelist(1001); // 1001

    size_t i = 0;
    std::generate(filelist.begin(), filelist.end(), [&i]() {
        return "ang1dnSi0p3mm_" + fmt::format("{:04}", i++) + ".h5";
    });

#ifdef ANGCAL_PLOT
    // uncomment if you want to see all frames and select a base peak
    /*
    PlotHelper temp_plotter(
        std::make_shared<AngleCalibration>(anglecalibration));
    plot_module_in_angle_for_all_frames(anglecalibration, temp_plotter,
                                        mythen_file_reader, filelist, 3);
    */
#endif

    const double base_peak_angle = -31.73; // 36.0568;
    // 35.550; // angpeak=69.225 - maybe select
    // one yourself !!! -could this be off?

    anglecalibration.set_base_peak_angle(base_peak_angle);

    // plot some stuff
    MythenFrame frame = mythen_file_reader.read_frame("ang1dnSi0p3mm_0170.h5");

#ifdef ANGCAL_PLOT
    plot_photon_counts(frame.photon_counts.view(),
                       {0, mythen_detector_ptr->num_strips()}, "Photon Counts",
                       mythen_detector_ptr);

    PlotHelper::pause();

    // normalize photon counts
    NDArray<double, 1> normalized_photon_counts{frame.photon_counts.shape()};
    for (ssize_t index = 0; index < frame.photon_counts.size(); ++index) {
        normalized_photon_counts(index) =
            frame.photon_counts(index) * inverse_normalized_flatfield(index);
    }

    plot_photon_counts(normalized_photon_counts.view(),
                       {0, mythen_detector_ptr->num_strips()},
                       "Normalized Photon Counts", mythen_detector_ptr);

    PlotHelper::pause();
#endif

// plot everything redistributed to fixed angle width bins
#ifdef ANGCAL_PLOT
    PlotHelper plotter(std::make_shared<AngleCalibration>(anglecalibration));

    auto new_fixed_angle_width_bins_photon_counts =
        anglecalibration.redistribute_photon_counts_to_fixed_angle_width_bins(
            frame);

    plotter.plot_redistributed_photon_counts(
        new_fixed_angle_width_bins_photon_counts.view());

    plotter.pause();
#endif

    size_t module_index = 3;
// plot one module for fixed angle width bins
#ifdef ANGCAL_PLOT
    // plot module
    module_index = 3;

    auto module_redistributed_to_fixed_angle_bins =
        get_redistributed_photon_counts(anglecalibration, module_index, frame);
    plotter.plot_module_redistributed_to_fixed_angle_width_bins(
        module_index, module_redistributed_to_fixed_angle_bins.view(),
        frame.detector_angle);

    plotter.pause();
#endif

    // plot base peak for one module

#ifdef ANGCAL_PLOT
    module_index = 3;
    auto base_peak_for_module =
        get_base_peak_region_of_interest(anglecalibration, module_index, frame);
    plotter.plot_base_peak_region_of_interest(module_index,
                                              base_peak_for_module.view());
#endif

    // anglecalibration.calibrate(filelist, base_peak_angle);
}
