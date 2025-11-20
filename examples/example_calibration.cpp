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
        "minimum module angle: {},{}, maximum module angle: {},{} for module "
        "{}",
        min_left_module_angle, min_right_module_angle, max_left_module_angle,
        max_right_module_angle, module_index);
}

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

    auto output_filename =
        file_path / "angcal_Jul2025_P12_0p0105_calibrated.off";

    // export
    // ANGCAL_TEST_DATA=~/Documents/VariaMay2025/Antonio20250512/angcal_M3_Mar21_2

    /*
    auto bad_channels_filename = file_path / "bcX.txt";

    auto flatfield_filename =
        file_path / "Flatfield_E17p5keV_T8751eV_MIX_Mar2021_open_WS.raw";

    auto initial_angles_filename = file_path / "angcal_Mar2021_P10.off";

    std::string acquisition_fileprefix = "ang1dnSi0p3mm_";

    auto output_filename = file_path / "angcal_Mar2021_P10_calibrated.off";
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

    // TODO: do I want to pass an inverse flatfield or only read flatfield and
    // calculate in Flatfield class
    //  inverse_flatfield is only properly allocated when calling
    //  inverse_flatfield

    LOG(TLogLevel::logINFO) << "read flatfield from file";

    NDArray<double, 2> normalized_flatfield(
        std::array<ssize_t, 2>{mythen_detector_ptr->num_strips(), 2});
    CustomFlatFieldFile flatfieldfilereader;
    flatfieldfilereader.open(flatfield_filename);
    flatfieldfilereader.read_into(normalized_flatfield.buffer(), 8);

    NDArray<double, 2> inverse_normalized_flatfield(
        std::array<ssize_t, 2>{mythen_detector_ptr->num_strips(), 2});
    for (ssize_t index = 0; index < mythen_detector_ptr->num_strips();
         ++index) {
        inverse_normalized_flatfield(index, 0) =
            1.0 / normalized_flatfield(index, 0);
        inverse_normalized_flatfield(index, 1) = normalized_flatfield(
            index, 1); // TODO: dont know how to handle variance
    }
    flat_field_ptr->set_inverse_normalized_flatfield(
        inverse_normalized_flatfield);

#ifdef ANGCAL_PLOT
    /*
    plot_photon_counts(flat_field_ptr->get_inverse_normalized_flatfield(),
                       {0, mythen_detector_ptr->num_strips()},
                       "Inverse normalized flatfield", mythen_detector_ptr);
    */
    /*
    plot_photon_counts(flat_field_ptr->get_inverse_normalized_flatfield(),
                       {3 * mythen_detector_ptr->strips_per_module(),
                        4 * mythen_detector_ptr->strips_per_module()},
                       "Inverse normalized flatfield for module 3",
                       mythen_detector_ptr);
    */
#endif

    auto mythen_file_reader = std::make_shared<CustomMythenFileReader>();

    AngleCalibration anglecalibration(mythen_detector_ptr, flat_field_ptr,
                                      mythen_file_reader);

    anglecalibration.read_initial_calibration_from_file(
        initial_angles_filename);

    LOG(TLogLevel::logINFO) << "read initial parameters from file";

    std::vector<std::string> filelist(1501); // 1501

    size_t i = 0;
    std::generate(filelist.begin(), filelist.end(),
                  [&i, &file_path, &acquisition_fileprefix]() {
                      return file_path / (acquisition_fileprefix +
                                          fmt::format("{:04}", i++) + ".h5");
                  });

    auto [left_angle, right_angle] =
        get_detector_range(mythen_file_reader, filelist);

    LOG(TLogLevel::logINFO)
        << fmt::format("detector range: [{},{}]", left_angle, right_angle);

#ifdef ANGCAL_PLOT
    PlotHelper plotter(std::make_shared<AngleCalibration>(anglecalibration));

    std::shared_ptr<Gnuplot> gp = std::make_shared<Gnuplot>();

    for (const auto &file : filelist) {
        auto frame = mythen_file_reader->read_frame(file);

        if (frame.detector_angle > 6.0 || frame.detector_angle > 33.0) {

            LOG(TLogLevel::logINFO) << "file: " << file;

            // plot module
            auto module_redistributed_to_fixed_angle_bins =
                anglecalibration
                    .redistribute_photon_counts_to_fixed_angle_width_bins(frame,
                                                                          0);

            plotter.plot_module_redistributed_to_fixed_angle_width_bins(
                0, module_redistributed_to_fixed_angle_bins.view(),
                frame.detector_angle, gp);

            plotter.pause();
        }
    }

#endif

    // take a tabulated peak as base peak
    // or take a base peak for module 0 that is well inside the detector range
    // and not at the module boundaries
    const double base_peak_angle = -49.4786;

    anglecalibration.set_base_peak_angle(base_peak_angle);

    anglecalibration.set_base_peak_ROI(0.18);

    // anglecalibration.set_histogram_bin_width(0.01);

    /*
    for (size_t module_index = 0; module_index < 48; ++module_index)
        get_angle_range_of_module(mythen_file_reader, filelist,
                                  anglecalibration, module_index);
    */

    // plot some stuff
    MythenFrame frame = mythen_file_reader->read_frame(
        file_path / acquisition_fileprefix.append("0850.h5")); // 0170.h5

    LOG(TLogLevel::logINFO) << fmt::format(
        "angle range of detector: {}",
        anglecalibration.diffraction_angle_from_DG_parameters(0, 0.0, 0, -0.5) -
            anglecalibration.diffraction_angle_from_DG_parameters(23, 0.0, 1279,
                                                                  +0.5));

    size_t module_index = 0;
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
        frame.detector_angle, gp);

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
    // anglecalibration.calibrate(filelist, base_peak_angle, output_filename);

    // anglecalibration.write_to_file(output_filename);
}
