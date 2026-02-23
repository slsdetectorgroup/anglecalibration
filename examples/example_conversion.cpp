#include "AngleCalibration.hpp"
#include "aare/File.hpp"
#include "logger.hpp"
#include "utils.hpp"
#include <filesystem>

/** How To use
 * set environmnet variable to where youre test data is located
 *  export ANGCAL_TEST_DATA=/home/mazzol_a/Documents/AngularConversionTestData
 */

int main() {

    // define test files names
    auto file_path = data_path();
    assert(std::filesystem::exists(file_path));

    auto bad_channels_filename = file_path / "bc2023_003_RING.chans";

    assert(std::filesystem::exists(bad_channels_filename));

    auto flatfield_filename =
        file_path /
        "Flatfield_E17p5keV_T12500eV_up_AUGCAL2_Sep2023_open_WS_C_X_X.raw";

    assert(std::filesystem::exists(flatfield_filename));

    auto initial_angles_filename = file_path / "Angcal_2E_Feb2023_P29.off";

    assert(std::filesystem::exists(initial_angles_filename));

    // cretae MythenDetectorSpecifications with all detector specifications
    std::shared_ptr<MythenDetectorSpecifications> mythen_detector_ptr =
        std::make_shared<MythenDetectorSpecifications>();

    mythen_detector_ptr->elastic_correction_factor = 0.0045;
    mythen_detector_ptr->detector_vertical_axis_offset = 24.8;

    // FileaReader to read Mythen Files
    auto mythen_file_reader = std::make_shared<EpicsMythenFileReader>();

    // create FlatField
    std::shared_ptr<FlatField> flat_field_ptr =
        std::make_shared<FlatField>(mythen_detector_ptr);

    flat_field_ptr->read_normalized_flatfield_from_file(flatfield_filename);

    LOG(TLogLevel::logINFO) << "normalized flatfield calculated";

    // create AngleCalibration class
    AngleCalibration anglecalibration(mythen_detector_ptr, flat_field_ptr,
                                      mythen_file_reader);

    // read initial parameters from file
    anglecalibration.read_initial_calibration_from_file(
        initial_angles_filename);

    anglecalibration.read_bad_channels_from_file(bad_channels_filename);

    LOG(TLogLevel::logINFO) << "read bad channels";

#ifdef ANGCAL_PLOT
    plot_photon_counts(anglecalibration.get_bad_channels(),
                       {0, mythen_detector_ptr->num_strips()}, "Bad channels",
                       std::nullopt);
#endif

    std::vector<std::string> file_list(4);
    std::generate(file_list.begin(), file_list.end(),
                  [n = 60, &file_path]() mutable {
                      std::string filename =
                          fmt::format("Fructose_0p2_60_{:04d}.h5", n++);
                      return (file_path / filename).string();
                  });

    anglecalibration.set_scale_factor(5050880.0);

    anglecalibration.set_angular_range(0.0, 90.5); // 130

    auto redistributed_photon_counts = anglecalibration.convert(
        file_list); // redistributes and applies pixel wise correction

    LOG(TLogLevel::logDEBUG) << "redistributed photon counts shape: "
                             << redistributed_photon_counts.shape()[0];

    LOG(TLogLevel::logDEBUG) << "num fixed angle bins: "
                             << anglecalibration.num_fixed_angle_width_bins();

#ifdef ANGCAL_PLOT
    MythenFrame frame =
        mythen_file_reader->read_frame(file_path / "Fructose_0p2_60_0060.h5");

    auto photon_counts = frame.photon_counts();

    plot_photon_counts(frame.photon_counts(),
                       {0, mythen_detector_ptr->num_strips()},
                       "Raw Photon Counts", mythen_detector_ptr);
#endif

    // plot
#ifdef ANGCAL_PLOT

    PlotHelper plotter(std::make_shared<AngleCalibration>(anglecalibration));

    plotter.plot_redistributed_photon_counts(
        redistributed_photon_counts.view());

    plotter.pause();
#endif
}