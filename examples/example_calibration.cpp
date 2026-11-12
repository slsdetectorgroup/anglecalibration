#include "AngleCalibration.hpp"
#include "AngularConversion.hpp"
#include "PlotHelpers.hpp"
#include "aare/File.hpp"
#include "logger.hpp"
#include "utils.hpp"
#include <filesystem>

/**
 * How to run: All data for this example can be found on
 * https://zenodo.org/records/20645666. Before running execute: export
 * ANGCAL_CALIBRATION_DATA=Path/to/mythenanglecalibrationtestdata
 */

using namespace angcal;

inline auto data_path() {
    if (const char *env_p = std::getenv("ANGCAL_CALIBRATION_DATA")) {
        return std::filesystem::path(env_p);
    } else {
        throw std::runtime_error(
            "Path to test data: $ANGCAL_CALIBRATION_DATA not set");
    }
}

int main() {

    auto file_path = data_path();

    assert(std::filesystem::exists(file_path));

    auto bad_channels_filename = file_path / "bc2025_001_RING.chans";

    auto flatfield_filename =
        file_path /
        "Flatfield_EkeV22p0_T11000eV_up_TESTFF1_clean_Jun2025_open_WS.raw";

    auto initial_angles_filename = file_path / "angcal_Jul2025_P12_0p0105.csv";

    std::string acquisition_fileprefix = "ang1up_22keV_MIX_0p5mm_48M_a_";

    const size_t num_files = 1501;

    // auto output_filename =
    // file_path / "angcal_Jul2025_P12_0p0105_new_calibrated.csv";

    const size_t num_modules = 48;

    DGParameters dgparameters(num_modules);

    read_csv(initial_angles_filename, dgparameters.parameters.view(), true);

    std::shared_ptr<MythenDetector> mythen_detector_ptr =
        std::make_shared<MythenDetector>(dgparameters);

    /*
    std::string unconnected_modules_filename = file_path / "ModOut.txt";

    auto unconnected_modules =
        read_unconnected_modules(unconnected_modules_filename);

    mythen_detector_ptr->set_unconnected_modules(unconnected_modules);
    */

    std::shared_ptr<FlatField> flat_field_ptr =
        std::make_shared<FlatField>(mythen_detector_ptr);

    LOG(TLogLevel::logINFO) << "read normalized flatfield from file";

    flat_field_ptr->set_normalized_flatfield(
        read_flatfield(flatfield_filename, num_modules));

    AngleCalibration anglecalibration(mythen_detector_ptr, flat_field_ptr);

    LOG(TLogLevel::logINFO) << "read bad channels";

    anglecalibration.set_bad_channels(
        read_bad_channels(bad_channels_filename,
                          mythen_detector_ptr->num_strips())
            .view());

    EpicsMythenFileReader mythen_file_reader(
        file_path / (acquisition_fileprefix + "0000.h5"));

    auto [left_angle, right_angle] = get_detector_range(mythen_file_reader);

    LOG(TLogLevel::logINFO) << fmt::format("detector range [degrees]: [{},{}]",
                                           left_angle, right_angle);

    LOG(TLogLevel::logINFO) << fmt::format(
        "range of detector [degrees]: {}",
        mythen_detector_ptr->diffraction_angle_from_DG_parameters(0, 0.0, 0,
                                                                  0.0) -
            mythen_detector_ptr->diffraction_angle_from_DG_parameters(
                23, 0.0, 1279, +1.0));

    auto [first_header, first_frame] = mythen_file_reader.read_frame();
    LOG(TLogLevel::logINFO)
        << fmt::format("incident intensity of first frame: {}",
                       first_header.incident_intensity);

    LOG(TLogLevel::logDEBUG) << "selecting base peak";

    // select_base_peak(std::make_shared<AngleCalibration>(anglecalibration),
    // mythen_file_reader, filelist, 0);

    // take a tabulated peak as base peak
    // or take a base peak for module 0 that is well inside the detector range
    // and not at the module boundaries

    const double base_peak_angle = 19.0678; // 19.0648;

    anglecalibration.set_base_peak_angle(base_peak_angle);

    anglecalibration.set_base_peak_ROI_width(0.18);

    // plot everything redistributed to fixed angle width bins
    PlotHelper plotter(mythen_detector_ptr,
                       anglecalibration.get_angular_range(),
                       anglecalibration.get_histogram_bin_width());

    auto [test_header, test_frame] = mythen_file_reader.read_frame(209);

    // plot some stuff
    double test_frame_angle =
        test_header
            .detector_angle; // mythen_file_reader->read_detector_angle(test_frame_name);

    HistogramMetaData histogram_meta_data{
        anglecalibration.get_angular_range(),
        anglecalibration.get_histogram_bin_width()};

    const double scale_factor = first_header.incident_intensity;

    // plot one module for fixed angle width bins
    auto module_redistributed_to_fixed_angle_bins = convert(
        *mythen_detector_ptr, histogram_meta_data, test_frame.photon_counts(),
        test_header.detector_angle, test_header.incident_intensity,
        test_header.exposure_time, anglecalibration.get_bad_channels(),
        flat_field_ptr->get_normalized_flatfield_view(), scale_factor);

    plotter.plot_diffraction_pattern(
        module_redistributed_to_fixed_angle_bins.view(), 1, test_frame_angle);

    // plot base peak for one module
    plotter.plot_base_peak(module_redistributed_to_fixed_angle_bins.view(),
                           base_peak_angle,
                           anglecalibration.get_base_peak_ROI_width());

    mythen_file_reader.seek(0); // reset to first frame

    auto [headers, frames] = mythen_file_reader.read_n(
        num_files); // TODO: maybe read_all would be great
    NDArray<double, 1> detector_angles(
        std::array<ssize_t, 1>{static_cast<ssize_t>(num_files)});
    std::transform(
        headers.begin(), headers.end(), detector_angles.begin(),
        [](const EpicsMythenHeader &header) { return header.detector_angle; });

    anglecalibration.calibrate<true>(frames.view(), detector_angles.view(),
                                     base_peak_angle, 1);

    /*
    anglecalibration.calibrate<false>(frames.view(), detector_angles.view(),
                                      base_peak_angle);

    auto bcparameters = anglecalibration.get_BCparameters();
    auto dgparameters = anglecalibration.get_DGparameters();
    bcparameters.convert_to_DGParameters(dgparameters);
    */
}
