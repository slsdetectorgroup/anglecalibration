/************************************************
 * @file MythenFileReader.test.cpp
 * @short test case for angle calibration class
 ***********************************************/

#include "MythenFileReader.hpp"

#include <filesystem>

#include "test_config.hpp"

#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace aare;
using namespace angcal;

TEST_CASE("MythenFileReader using raw file", "[.mythenfilereader][.files]") {

    SECTION("read mythen frame from raw file") {

        MythenFileReader file_reader;

        auto mythen_data = test_data_path() / "AngleCalibration_Test_Data" /
                           "DummyRawFiles" / "run_2_master_2.json";

        REQUIRE(std::filesystem::exists(mythen_data));

        MythenFrame frame = file_reader.read_frame(mythen_data);

        REQUIRE(frame.detector_angle == 0);
        REQUIRE(frame.incident_intensity == 0);
        REQUIRE(frame.photon_counts().size() == 1280 * 2);
        REQUIRE(frame.channel_mask == std::array<uint8_t, 3>{0, 0, 1});
    }
    SECTION(
        "read mythen frame from raw file including detector and I0 values") {

        auto mythen_data = test_data_path() / "AngleCalibration_Test_Data" /
                           "DummyRawFiles" / "run_2_master_2.json";

        REQUIRE(std::filesystem::exists(mythen_data));

        auto I0_data = test_data_path() / "AngleCalibration_Test_Data" /
                       "DummyRawFiles" / "I0.bin";

        auto detector_angles = test_data_path() / "AngleCalibration_Test_Data" /
                               "DummyRawFiles" / "detector_angles.bin";

        REQUIRE(std::filesystem::exists(I0_data));
        REQUIRE(std::filesystem::exists(detector_angles));

        MythenFileReader file_reader(detector_angles, I0_data);

        MythenFrame frame = file_reader.read_frame(mythen_data);

        REQUIRE(frame.detector_angle == 1.11);
        REQUIRE(frame.incident_intensity == 1002);
        REQUIRE(frame.photon_counts().size() == 1280 * 2);
        REQUIRE(frame.channel_mask == std::array<uint8_t, 3>{0, 0, 1});
    }
}

TEST_CASE("test custom mythenfilereader", "[.mythenfilereader][.files]") {

    auto fpath = test_data_path() / "AngleCalibration_Test_Data";

    REQUIRE(std::filesystem::exists(fpath));

    CustomMythenFileReader file_reader;

    auto frame =
        file_reader.read_frame(fpath / "ang1up_22keV_LaB60p3mm_48M_a_0320.h5");

    CHECK(frame.detector_angle == 0.99955);

    CHECK(frame.channel_mask == std::array<uint8_t, 3>{0, 0, 1});

    CHECK(frame.photon_counts().size() == 61440);

    CHECK(frame.incident_intensity == 0.0);
}
