/************************************************
 * @file AngleCalibration.test.cpp
 * @short test case for angle calibration class
 ***********************************************/

#include "AngleCalibration.hpp"
#include "CustomFiles.hpp"
#include "FlatField.hpp"
#include "logger.hpp"

#include <filesystem>

#include "test_config.hpp"

#include <iomanip>
#include <type_traits>

#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace aare;
using namespace angcal;

TEST_CASE("read initial angle calibration file",
          "[anglecalibration] [.files]") {

    std::shared_ptr<MythenDetectorSpecifications> mythen_detector_ptr =
        std::make_shared<MythenDetectorSpecifications>();

    auto fpath = test_data_path() / "AngleCalibration_Test_Data";

    AngleCalibration anglecalibration(mythen_detector_ptr,
                                      std::shared_ptr<FlatField>{});

    std::string filename = fpath / "Angcal_2E_Feb2023_P29.off";

    REQUIRE(std::filesystem::exists(filename));

    anglecalibration.read_initial_calibration_from_file(filename);

    auto parameters = anglecalibration.get_DGparameters();

    std::cout.precision(17);

    CHECK(parameters.parameters.size() == 48 * 3);

    CHECK(parameters.centers(9) == Catch::Approx(660.342326));
    CHECK(parameters.offsets(47) == Catch::Approx(5.8053312));
    CHECK(parameters.conversions(27) == Catch::Approx(-0.6581179125e-4));
}

TEST_CASE("read bad channels",
          "[anglecalibration][mythenspecifications][.files]") {

    MythenDetectorSpecifications mythen_detector(
        std::make_shared<CustomBadChannelsFile>());

    std::string bad_channels_filename = test_data_path() /
                                        "AngleCalibration_Test_Data" /
                                        "bc2023_003_RING.chans";

    REQUIRE(std::filesystem::exists(bad_channels_filename));

    mythen_detector.read_bad_channels_from_file(bad_channels_filename);

    CHECK(mythen_detector.get_bad_channels().size() == 61440);

    CHECK(mythen_detector.get_bad_channels()[61437] == true);
    CHECK(std::all_of(mythen_detector.get_bad_channels().begin() + 30720,
                      mythen_detector.get_bad_channels().begin() + 61439,
                      [](const bool element) { return element; }));
}

TEST_CASE("read flatfield", "[anglecalibration][flatfield][.files]") {

    std::shared_ptr<MythenDetectorSpecifications> mythen_detector_ptr =
        std::make_shared<MythenDetectorSpecifications>();

    FlatField flatfield(mythen_detector_ptr);

    std::string flatfield_filename =
        test_data_path() / "AngleCalibration_Test_Data" /
        "Flatfield_E22p0keV_T11000eV_up_48M_a_LONG_Feb2023_open_WS_SUMC.raw";

    REQUIRE(std::filesystem::exists(flatfield_filename));

    flatfield.read_flatfield_from_file(flatfield_filename);

    auto flatfield_data = flatfield.get_flatfield();

    CHECK(flatfield_data.size() == 61440);

    CHECK(flatfield_data[0] == 0);
    CHECK(flatfield_data[21] == 4234186);
}

TEST_CASE("create flatfield", "[anglecalibration], [flatfield], [.files]") {

    ssize_t n_modules = 2;
    double exposure_time = 5.0;
    size_t n_counters = 1;
    std::shared_ptr<MythenDetectorSpecifications> mythen_detector_ptr =
        std::make_shared<MythenDetectorSpecifications>(n_modules, exposure_time,
                                                       n_counters);

    FlatField flatfield(mythen_detector_ptr);

    auto data_path = test_data_path() / "AngleCalibration_Test_Data" /
                     "Flatfieldacquisition";

    REQUIRE(std::filesystem::exists(data_path));

    CHECK_NOTHROW(flatfield.create_flatfield_from_rawfilesystem(data_path));

    auto flatfield_data = flatfield.get_flatfield();

    CHECK(flatfield_data.size() == n_modules * 1280);

    CHECK(flatfield_data[0] == 0);
    CHECK(flatfield_data[21] ==
          2 * 3 *
              21); // virtual data 2 angles, 3 frames - 21 as increasing numbers
}

TEST_CASE("check parameter conversion", "[anglecalibration]") {

    std::shared_ptr<MythenDetectorSpecifications> mythen_detector_ptr =
        std::make_shared<MythenDetectorSpecifications>();

    auto file_path = test_data_path();

    AngleCalibration anglecalibration(mythen_detector_ptr,
                                      std::shared_ptr<FlatField>{});

    std::string initial_angles_filename =
        test_data_path() / "AngleCalibration_Test_Data" / "AntoniosTestData" /
        "angcal_Mar2021_P10.off";

    REQUIRE(std::filesystem::exists(initial_angles_filename));

    anglecalibration.read_initial_calibration_from_file(
        initial_angles_filename);

    SECTION("clockwise module") {
        const ssize_t local_strip_index = 1;
        const ssize_t module_index = 0;

        double diffraction_angle_DG_param =
            anglecalibration.diffraction_angle_from_DG_parameters(
                module_index, 0.0, local_strip_index);

        double diffraction_angle_BC_param =
            anglecalibration.diffraction_angle_from_BC_parameters(
                module_index, 0.0, local_strip_index);

        auto [module_center_distance, normal_distance, angle] =
            anglecalibration.get_DGparameters().convert_to_EEParameters(
                module_index);

        double diffraction_angle_EE_param =
            anglecalibration.diffraction_angle_from_EE_parameters(
                module_center_distance, normal_distance, angle, 0.0,
                local_strip_index);

        CHECK(diffraction_angle_EE_param ==
              Catch::Approx(diffraction_angle_DG_param));

        CHECK(diffraction_angle_BC_param ==
              Catch::Approx(diffraction_angle_DG_param));
    }
    SECTION("counter clockwise module") {
        const ssize_t local_strip_index = 1278;
        const ssize_t module_index = 47;

        double diffraction_angle_DG_param =
            anglecalibration.diffraction_angle_from_DG_parameters(
                module_index, 0.0, local_strip_index);

        double diffraction_angle_BC_param =
            anglecalibration.diffraction_angle_from_BC_parameters(
                module_index, 0.0, local_strip_index);

        auto [module_center_distance, normal_distance, angle] =
            anglecalibration.get_DGparameters().convert_to_EEParameters(
                module_index);

        double diffraction_angle_EE_param =
            anglecalibration.diffraction_angle_from_EE_parameters(
                module_center_distance, normal_distance, angle, 0.0,
                local_strip_index);

        CHECK(diffraction_angle_EE_param ==
              Catch::Approx(diffraction_angle_DG_param));

        CHECK(diffraction_angle_BC_param ==
              Catch::Approx(diffraction_angle_DG_param));
    }
    SECTION("clockwise and counterclockwise module have the same diffraction "
            "angle") {
        double diffraction_angle_BC_param_clockwise_module =
            anglecalibration.diffraction_angle_from_BC_parameters(0, 0.0, 1);

        double diffraction_angle_BC_param_counterclockwise_module =
            anglecalibration.diffraction_angle_from_BC_parameters(47, 0.0,
                                                                  1278);

        // TODO: not sure if this should hold
        CHECK(
            diffraction_angle_BC_param_clockwise_module ==
            Catch::Approx(diffraction_angle_BC_param_counterclockwise_module));
    }
}
