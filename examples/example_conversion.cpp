#include "AngularConversion.hpp"
#include "aare/File.hpp"
#include "logger.hpp"
#include "utils.hpp"
#include <filesystem>

/** How To use
 * set environmnet variable to where youre test data is located
 *  export ANGCAL_CONVERSION_DATA=/path/to/local/mythenangleconversiontestdata
 */

using namespace angcal;

inline auto data_path() {
    if (const char *env_p = std::getenv("ANGCAL_CONVERSION_DATA")) {
        return std::filesystem::path(env_p);
    } else {
        throw std::runtime_error(
            "Path to test data: $ANGCAL_CONVERSION_DATA not set");
    }
}

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

    const size_t num_modules = 48;

    DGParameters dgparameters(num_modules);
    read_csv(initial_angles_filename, dgparameters.parameters.view(), true);

    // create MythenDetector with all detector specifications and geometry
    std::shared_ptr<MythenDetector> mythen_detector_ptr =
        std::make_shared<MythenDetector>(dgparameters);

    mythen_detector_ptr->elastic_correction_factor(0.0045);
    mythen_detector_ptr->detector_vertical_axis_offset(24.8);

    NDArray<double, 2> flatfield =
        read_flatfield(flatfield_filename, mythen_detector_ptr->num_modules());

    LOG(TLogLevel::logINFO) << "read normalized flatfield from file";

    HistogramMetaData histogram_meta_data{std::make_tuple(0.0, 90.5), 0.0036};

    NDArray<bool, 1> bad_channels = read_bad_channels(
        bad_channels_filename, mythen_detector_ptr->num_strips());

    LOG(TLogLevel::logINFO) << "read bad channels";

    EpicsMythenFileReader mythen_file_reader(file_path /
                                             "Fructose_0p2_60_0060.h5");

    auto [headers, frames] = mythen_file_reader.read_n(4);

    double exposure_time = headers[0].exposure_time;

    NDArray<double, 1> detector_angles(
        std::array<ssize_t, 1>{static_cast<ssize_t>(headers.shape(0))});
    NDArray<uint64_t, 1> incident_intensities(
        std::array<ssize_t, 1>{static_cast<ssize_t>(headers.shape(0))});
    for (ssize_t frame_index = 0; frame_index < headers.shape(0);
         ++frame_index) {

        detector_angles[frame_index] = headers[frame_index].detector_angle;

        incident_intensities[frame_index] =
            headers[frame_index].incident_intensity;
    }

    const double scale_factor = incident_intensities[0];

    auto redistributed_photon_counts = convert(
        *mythen_detector_ptr, histogram_meta_data, frames.view(),
        detector_angles.view(), incident_intensities.view(), exposure_time,
        bad_channels.view(), flatfield.view(), scale_factor);

    // plot
    PlotHelper plotter(mythen_detector_ptr, histogram_meta_data.angular_range,
                       histogram_meta_data.bin_width);

    plotter.plot_diffraction_pattern(redistributed_photon_counts.view());
}