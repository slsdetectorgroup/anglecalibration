#include "AngleCalibration.hpp"
#include "PlotHelpers.hpp"
#include "aare/File.hpp"
#include "logger.hpp"
#include "utils.hpp"
#include <filesystem>

/**
 * How to run: All data for this example can be found in the github directory
 * https://gitea.psi.ch/angcal/VariaMay2025 before running execute: export
 * ANGCAL_TEST_PATH=Path_to_git_repository
 *
 */

using namespace angcal;

int main() {

    auto file_path = data_path();

    assert(std::filesystem::exists(file_path));

    // export ANGCAL_TEST_DATA=/home/mazzol_a/ANGCALDATA/angcal22_Jul25/

    auto bad_channels_filename = file_path / "bc2025_001_RING.chans";

    auto flatfield_filename =
        file_path /
        "Flatfield_EkeV22p0_T11000eV_up_TESTFF1_clean_Jun2025_open_WS.raw";

    auto initial_angles_filename = file_path / "angcal_Jul2025_P12_0p0105.off";

    std::string acquisition_fileprefix = "ang1up_22keV_MIX_0p5mm_48M_a_";

    const size_t num_files = 1501;

    auto output_filename =
        file_path / "angcal_Jul2025_P12_0p0105_new_calibrated.off";

    // export
    // ANGCAL_TEST_DATA=~/Documents/VariaMay2025/Antonio20250512/angcal_M3_Mar21_2

    /*
    auto bad_channels_filename = file_path / "bcX.txt";

    auto flatfield_filename =
        file_path / "Flatfield_E17p5keV_T8751eV_MIX_Mar2021_open_WS.raw";

    auto initial_angles_filename = file_path / "angcal_Mar2021_P10.off";

    std::string acquisition_fileprefix = "ang1upSi0p3mm_";

    const size_t num_files = 1001;

    auto output_filename = file_path / "angcal_Mar2021_P10_calibrated.off";

    // base peaks: 31.588, 23.5126, 28.771
    */

    std::shared_ptr<MythenDetectorSpecifications> mythen_detector_ptr =
        std::make_shared<MythenDetectorSpecifications>();

    assert(std::filesystem::exists(bad_channels_filename));

    /*
    std::string unconnected_modules_filename = file_path / "ModOut.txt";

    auto unconnected_modules =
        read_unconnected_modules(unconnected_modules_filename);

    mythen_detector_ptr->set_unconnected_modules(unconnected_modules);
    */

    std::shared_ptr<FlatField> flat_field_ptr =
        std::make_shared<FlatField>(mythen_detector_ptr);

    LOG(TLogLevel::logINFO) << "read normalized flatfield from file";

    flat_field_ptr->read_normalized_flatfield_from_file(flatfield_filename);

    auto mythen_file_reader = std::make_shared<EpicsMythenFileReader>();

    AngleCalibration anglecalibration(mythen_detector_ptr, flat_field_ptr,
                                      mythen_file_reader);

    anglecalibration.read_initial_calibration_from_file(
        initial_angles_filename);

    LOG(TLogLevel::logINFO) << "read initial parameters from file";

    anglecalibration.read_bad_channels_from_file(bad_channels_filename);

    // plot("Bad Channels", anglecalibration.get_bad_channels());

    LOG(TLogLevel::logINFO) << "read bad channels";

    std::vector<std::string> filelist(num_files); // 1501

    size_t i = 0;
    std::generate(filelist.begin(), filelist.end(),
                  [&i, &file_path, &acquisition_fileprefix]() {
                      return file_path / (acquisition_fileprefix +
                                          fmt::format("{:04}", i++) + ".h5");
                  });

    auto [left_angle, right_angle] =
        get_detector_range(mythen_file_reader, filelist);

    LOG(TLogLevel::logINFO) << fmt::format("detector range [degrees]: [{},{}]",
                                           left_angle, right_angle);

    LOG(TLogLevel::logINFO) << fmt::format(
        "range of detector [degrees]: {}",
        anglecalibration.diffraction_angle_from_DG_parameters(0, 0.0, 0, -0.5) -
            anglecalibration.diffraction_angle_from_DG_parameters(23, 0.0, 1279,
                                                                  +0.5));

    LOG(TLogLevel::logDEBUG) << "selecting base peak";

    select_base_peak(std::make_shared<AngleCalibration>(anglecalibration),
                     mythen_file_reader, filelist, 0);

    // take a tabulated peak as base peak
    // or take a base peak for module 0 that is well inside the detector range
    // and not at the module boundaries

    // TODO: choosing the base peak is very hard - good average
    const double base_peak_angle = 19.0678; // 19.0648;

    anglecalibration.set_base_peak_angle(base_peak_angle);

    anglecalibration.set_base_peak_ROI_width(0.18);

    auto first_frame = mythen_file_reader->read_frame(filelist[0]);
    LOG(TLogLevel::logINFO)
        << fmt::format("incident intensity of first frame: {}",
                       first_frame.incident_intensity);
    anglecalibration.set_scale_factor(first_frame.incident_intensity);

    // anglecalibration.set_histogram_bin_width(0.01);

    // plot everything redistributed to fixed angle width bins
    PlotHelper plotter(std::make_shared<AngleCalibration>(anglecalibration));

    auto test_frame_name = file_path / (acquisition_fileprefix + "0209.h5");

    // plot some stuff
    double test_frame_angle =
        mythen_file_reader->read_detector_angle(test_frame_name);

    // plot one module for fixed angle width bins
    auto module_redistributed_to_fixed_angle_bins =
        anglecalibration.convert({test_frame_name}, 1);

    plotter.plot_diffraction_pattern(
        test_frame_angle, module_redistributed_to_fixed_angle_bins.view(), 1);

    // plot base peak for one module
    plotter.plot_base_peak(module_redistributed_to_fixed_angle_bins.view());

    // anglecalibration.calibrate<true>(filelist, base_peak_angle, 1);

    /*
    #ifdef ANGCAL_PLOT
        anglecalibration.set_base_peak_angle(base_peak_angle);
        anglecalibration.set_calibration_files(filelist);
        anglecalibration.read_initial_calibration_from_file(output_filename);
        std::string plot_title = fmt::format("Overlaping Base Peaks");
        auto plot = std::make_shared<PlotCalibrationProcess>(plot_title);

        anglecalibration.plot_all_base_peaks(plot);
    #endif
    */

    /*
    anglecalibration.calibrate(filelist, base_peak_angle);

    auto bcparameters = anglecalibration.get_BCparameters();
    auto dgparameters = anglecalibration.get_DGparameters();
    bcparameters.convert_to_DGParameters(dgparameters);

    anglecalibration.write_DG_parameters_to_file(output_filename, dgparameters);
    */
}
