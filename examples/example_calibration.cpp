#include "AngleCalibration.hpp"
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
        file_path / "angcal_Jul2025_P12_0p0105_calibrated.off";

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

    mythen_detector_ptr->read_bad_channels_from_file(bad_channels_filename);

    LOG(TLogLevel::logINFO) << "read bad channels";

#ifdef ANGCAL_PLOT

    plot_photon_counts(mythen_detector_ptr->get_bad_channels(),
                       {0, mythen_detector_ptr->num_strips()}, "Bad channels",
                       std::nullopt);
#endif

    /*
    std::string unconnected_modules_filename = file_path / "ModOut.txt";

    auto unconnected_modules =
        read_unconnected_modules(unconnected_modules_filename);

    mythen_detector_ptr->set_unconnected_modules(unconnected_modules);
    */

    std::shared_ptr<FlatField> flat_field_ptr =
        std::make_shared<FlatField>(mythen_detector_ptr);

    LOG(TLogLevel::logINFO) << "read normalized flatfield from file";

    // TODO: probably normalized flatfield is stored or inverse?
    // TODO: maybe allocate normalized flatfield somewhere else and directly
    // read into it- decouble filesystem from FlatField class completely
    // Flatfield and normalized flatfield use different file formats e.g. no
    // error stored in flatfield - flatfield is a simple mythenfile!!!! - not a
    // txt file
    flat_field_ptr->read_normalized_flatfield_from_file(flatfield_filename);

    flat_field_ptr->calculate_inverse_normalized_flatfield<true>();

    auto mythen_file_reader = std::make_shared<EpicsMythenFileReader>();

    AngleCalibration anglecalibration(mythen_detector_ptr, flat_field_ptr,
                                      mythen_file_reader);

    anglecalibration.read_initial_calibration_from_file(
        initial_angles_filename);

    LOG(TLogLevel::logINFO) << "read initial parameters from file";

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

    // get_module_angle_ranges(mythen_file_reader, filelist, anglecalibration,
    // mythen_detector_ptr->max_modules());

#ifdef ANGCAL_PLOT
    // select_base_peak(std::make_shared<AngleCalibration>(anglecalibration),
    // mythen_file_reader, filelist, 18);
#endif

    // take a tabulated peak as base peak
    // or take a base peak for module 0 that is well inside the detector range
    // and not at the module boundaries
    const double base_peak_angle = 17.599; // 26.7731, 20.5902, 14.0686

    anglecalibration.set_base_peak_angle(base_peak_angle);

    anglecalibration.set_base_peak_ROI_width(0.18);

    // anglecalibration.set_histogram_bin_width(0.01);

    // plot some stuff
    MythenFrame frame = mythen_file_reader->read_frame(
        file_path / acquisition_fileprefix.append("1014.h5")); // 0165.h5

    size_t module_index = 18; // 0

#ifdef ANGCAL_PLOT

    // normalize photon counts
    NDArray<double, 1> normalized_photon_counts{
        std::array<ssize_t, 1>{frame.size()}, 0.0};

    for (ssize_t index = 0; index < frame.size(); ++index) {
        normalized_photon_counts(index) =
            frame.photon_counts(index) *
            flat_field_ptr->get_inverse_normalized_flatfield()(index, 0);
    }

    plot_photon_counts(normalized_photon_counts.view(),
                       {0, mythen_detector_ptr->num_strips()},
                       "Normalized Photon Counts", mythen_detector_ptr);

    plot_photon_counts(
        normalized_photon_counts.view(),
        {module_index * mythen_detector_ptr->strips_per_module(),
         (module_index + 1) * mythen_detector_ptr->strips_per_module()},
        fmt::format("Normalized Photon Counts for module {}", module_index),
        mythen_detector_ptr);

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

// plot one module for fixed angle width bins
#ifdef ANGCAL_PLOT
    // plot module

    auto module_redistributed_to_fixed_angle_bins =
        anglecalibration.redistribute_photon_counts_to_fixed_angle_width_bins(
            frame, module_index);

    plotter.plot_module_redistributed_to_fixed_angle_width_bins(
        module_index, module_redistributed_to_fixed_angle_bins.view(),
        frame.detector_angle);

    plotter.pause();
#endif

// plot base peak for one module
#ifdef ANGCAL_PLOT

    auto base_peak_for_module =
        anglecalibration.redistributed_photon_counts_in_base_peak_ROI(
            frame, module_index);

    plotter.plot_base_peak_region_of_interest(module_index,
                                              base_peak_for_module.view());
    plotter.pause();

#endif

    // anglecalibration.calibrate(filelist, base_peak_angle, module_index);
    anglecalibration.calibrate(filelist, base_peak_angle); //, output_filename);

    // anglecalibration.write_to_file(output_filename);
}
