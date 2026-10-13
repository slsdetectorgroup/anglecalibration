#include "AngularConversion.hpp"
#include "aare/File.hpp"
#include "logger.hpp"
#include "utils.hpp"
#include <filesystem>

/** How To use
 * set environmnet variable to where youre test data is located
 *  export ANGCAL_TEST_DATA=/path/to/local/mythenangleconversiontestdata
 */

using namespace angcal;

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

    auto initial_angles_filename = file_path / "Angcal_2E_Feb2023_P29.csv";

    assert(std::filesystem::exists(initial_angles_filename));

    // create MythenDetector with all detector specifications and geometry
    std::shared_ptr<MythenDetector> mythen_detector_ptr =
        std::make_shared<MythenDetector>();

    mythen_detector_ptr->elastic_correction_factor(0.0045);
    mythen_detector_ptr->detector_vertical_axis_offset(24.8);

    DGParameters dgparameters(mythen_detector_ptr->num_modules());
    read_csv(initial_angles_filename, dgparameters.parameters.view(), true);
    mythen_detector_ptr->set_DGparameters(
        dgparameters); // set detector parameters

    // FileaReader to read Mythen Files
    auto mythen_file_reader = std::make_shared<EpicsMythenFileReader>();

    NDArray<double, 2> flatfield =
        read_flatfield(flatfield_filename, mythen_detector_ptr->num_modules());

    LOG(TLogLevel::logINFO) << "read normalized flatfield from file";

    // create AngleCalibration class
    AngularConversion angularconversion(mythen_detector_ptr, 0.0036,
                                        std::make_tuple(0.0, 90.5)); // 130

    NDArray<bool, 1> bad_channels = read_bad_channels(
        bad_channels_filename, mythen_detector_ptr->num_strips());

    angularconversion.set_bad_channels(bad_channels.view());

    angularconversion.set_scale_factor(5050880.0);

    LOG(TLogLevel::logINFO) << "read bad channels";

    std::vector<std::string> file_list(4);
    std::generate(file_list.begin(), file_list.end(),
                  [n = 60, &file_path]() mutable {
                      std::string filename =
                          fmt::format("Fructose_0p2_60_{:04d}.h5", n++);
                      return (file_path / filename).string();
                  });

    const ssize_t num_frames = file_list.size();

    NDArray<uint32_t, 2> photon_counts(std::array<ssize_t, 2>{
        num_frames, static_cast<ssize_t>(mythen_detector_ptr->num_strips())});
    NDArray<double, 1> detector_angles(std::array<ssize_t, 1>{num_frames});
    NDArray<uint64_t, 1> incident_intensities(
        std::array<ssize_t, 1>{num_frames});

    double exposure_time =
        mythen_file_reader->read_frame(file_list[0]).exposure_time;

    for (ssize_t frame_index = 0; frame_index < num_frames; ++frame_index) {

        mythen_file_reader->read_into(photon_counts.buffer() +
                                          frame_index *
                                              photon_counts.byte_strides()[0],
                                      file_list[frame_index]);

        detector_angles[frame_index] =
            mythen_file_reader->read_detector_angle(file_list[frame_index]);

        incident_intensities(frame_index) =
            mythen_file_reader->read_incident_intensity(file_list[frame_index]);
    }

    angularconversion.set_scale_factor(incident_intensities[0]);

    auto redistributed_photon_counts =
        angularconversion.convert(photon_counts.view(), detector_angles.view(),
                                  incident_intensities.view(), exposure_time);

    // plot
    PlotHelper plotter(mythen_detector_ptr,
                       angularconversion.get_angular_range(),
                       angularconversion.get_histogram_bin_width());

    plotter.plot_diffraction_pattern(redistributed_photon_counts.view());
}