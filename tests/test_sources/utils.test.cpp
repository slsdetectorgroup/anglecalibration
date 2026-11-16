/************************************************
 * @file utils.test.cpp
 * @short test cases for the file reading functions in utils.hpp
 ***********************************************/

#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "TempFileDirectory.hpp"
#include "test_config.hpp"
#include "utils.hpp"

using namespace aare;
using namespace angcal;

// ---------------------------------------------------------------------------
// bad channels (.chans): one strip index or an inclusive range "a-b" per line
// ---------------------------------------------------------------------------

TEST_CASE("read bad channels from temporary .chans file", "[utils]") {

    TemporaryDirectory tmp_dir("angcal-bad_channels");

    constexpr ssize_t num_strips = 2 * MythenDetector::strips_per_module();

    SECTION("single indices and inclusive ranges") {
        // same layout as bc2023_003_RING.chans
        // second module completely bad, first module has some bad strips
        const auto filename =
            tmp_dir.write_file("bad_channels.chans", "0\n"
                                                     "1\n"
                                                     "5\n"
                                                     "1274\n"
                                                     "1279\n"
                                                     "1280-2559\n");

        auto bad_channels = read_bad_channels(filename, num_strips);

        REQUIRE(bad_channels.size() == num_strips);

        std::set<ssize_t> expected_bad{0, 1, 5, 1274, 1279};
        for (ssize_t i = 1280; i < 2560; ++i) {
            expected_bad.insert(i);
        }

        for (ssize_t strip = 0; strip < num_strips; ++strip) {
            const bool is_expected_bad =
                expected_bad.find(strip) != expected_bad.end();
            INFO("strip " << strip);
            CHECK(bad_channels(strip) == is_expected_bad);
        }

        const auto num_bad =
            std::count(bad_channels.begin(), bad_channels.end(), true);
        CHECK(num_bad == static_cast<ssize_t>(expected_bad.size()));
    }

    SECTION("out of bound range throws") {
        const auto filename =
            tmp_dir.write_file("out_of_bounds.chans", "1\n"
                                                      "1280-2559\n"
                                                      "2560\n");

        CHECK_THROWS_AS(read_bad_channels(filename, num_strips),
                        std::runtime_error);
    }

    SECTION("empty file yields no bad channels") {
        const auto filename = tmp_dir.write_file("empty.chans", "");

        auto bad_channels = read_bad_channels(filename, num_strips);

        CHECK(bad_channels.size() == num_strips);
        CHECK(std::none_of(bad_channels.begin(), bad_channels.end(),
                           [](bool element) { return element; }));
    }

    SECTION("non existing file throws") {
        CHECK_THROWS_AS(
            read_bad_channels(tmp_dir.path() / "does_not_exist.chans",
                              num_strips),
            std::runtime_error);
    }
}

// ---------------------------------------------------------------------------
// flatfield (.raw): header line, then whitespace separated
// strip_index flatfield flatfield_error val_1 val_2
// ---------------------------------------------------------------------------

TEST_CASE("read flatfield from temporary .raw file", "[utils]") {

    TemporaryDirectory tmp_dir("angcal-flatfield");

    constexpr ssize_t num_modules = 2;
    constexpr ssize_t num_strips =
        num_modules * MythenDetector::strips_per_module();

    SECTION("values in same format as Flatfield_*.raw") {
        // mixes the formatting found in the real file: left aligned indices
        // for bad strips (-1), right aligned indices and Fortran style
        // exponents
        const auto filename = tmp_dir.write_file(
            "flatfield.raw",
            "# chan eff err_eff\n"
            "0  -1.0000000000000000       -1.0000000000000000        "
            "6.9028838384535474E-007   3.7617109252629684E-003\n"
            "5  -1.0000000000000000       -1.0000000000000000        "
            "6.9031705898936257E-007   3.7618151015469843E-003\n"
            "           6  0.76053714144724005        "
            "3.9644285189937315E-003   6.9032276742712450E-007   "
            "3.7618358400269827E-003\n"
            "        1280  1.2500000000000000        "
            "1.0000000000000000E-002   6.9032276742712450E-007   "
            "3.7618358400269827E-003\n"
            "2559   -1.0000000000000000       -1.0000000000000000       "
            "-6.8617006362804965E-007   3.7505050460968589E-003\n");

        auto flatfield = read_flatfield(filename, num_modules);

        REQUIRE(flatfield.shape(0) == num_strips);
        REQUIRE(flatfield.shape(1) == 2);

        CHECK(flatfield(0, 0) == Catch::Approx(-1.0));
        CHECK(flatfield(0, 1) == Catch::Approx(-1.0));
        CHECK(flatfield(5, 0) == Catch::Approx(-1.0));
        CHECK(flatfield(6, 0) == Catch::Approx(0.76053714144724005));
        CHECK(flatfield(6, 1) == Catch::Approx(3.9644285189937315E-003));
        CHECK(flatfield(1280, 0) == Catch::Approx(1.25));
        CHECK(flatfield(1280, 1) == Catch::Approx(1.0E-002));
        CHECK(flatfield(2559, 0) == Catch::Approx(-1.0));
        CHECK(flatfield(2559, 1) == Catch::Approx(-1.0));

        // strips not listed in the file are zero initialized
        CHECK(flatfield(1, 0) == 0.0);
        CHECK(flatfield(1, 1) == 0.0);
        CHECK(flatfield(7, 0) == 0.0);
        CHECK(flatfield(2558, 1) == 0.0);
    }

    SECTION("only header yields zero flatfield") {
        const auto filename =
            tmp_dir.write_file("header_only.raw", "# chan eff err_eff\n");

        auto flatfield = read_flatfield(filename, num_modules);

        CHECK(flatfield.shape(0) == num_strips);
        CHECK(std::all_of(flatfield.begin(), flatfield.end(),
                          [](double element) { return element == 0.0; }));
    }

    SECTION("non existing file throws") {
        CHECK_THROWS_AS(
            read_flatfield(tmp_dir.path() / "does_not_exist.raw", num_modules),
            std::runtime_error);
    }
}

// ---------------------------------------------------------------------------
// csv: header "center,conversion,offset" and one row per module
// ---------------------------------------------------------------------------

namespace {
const std::string dg_parameters_csv =
    "center,conversion,offset\n"
    "643.809720708243,6.57719636736681e-05,0.0\n"
    "642.197591224993,6.57694036246975e-05,5.00489288125167\n"
    "639.5,-6.581179125986344e-05,5.805331278507036\n";
} // namespace

TEST_CASE("read csv into preallocated buffer", "[utils][csv]") {

    TemporaryDirectory tmp_dir("angcal-csv");

    SECTION("values in same format as Angcal_*.csv") {
        const auto filename =
            tmp_dir.write_file("parameters.csv", dg_parameters_csv);

        NDArray<double, 2> parameters({3, 3}, 0.0);
        read_csv(filename, parameters.view());

        CHECK(parameters(0, 0) == Catch::Approx(643.809720708243));
        CHECK(parameters(0, 1) == Catch::Approx(6.57719636736681e-05));
        CHECK(parameters(0, 2) == Catch::Approx(0.0));
        CHECK(parameters(1, 2) == Catch::Approx(5.00489288125167));
        CHECK(parameters(2, 0) == Catch::Approx(639.5));
        CHECK(parameters(2, 1) == Catch::Approx(-6.581179125986344e-05));
        CHECK(parameters(2, 2) == Catch::Approx(5.805331278507036));
    }

    SECTION("empty lines and windows line endings are ignored") {
        const auto filename =
            tmp_dir.write_file("crlf.csv", "center,conversion,offset\r\n"
                                           "1.0,2.0,3.0\r\n"
                                           "\r\n"
                                           "\n"
                                           "4.0,5.0,6.0\r\n");

        NDArray<double, 2> parameters({2, 3}, 0.0);
        read_csv(filename, parameters.view());

        CHECK(parameters(0, 0) == 1.0);
        CHECK(parameters(0, 2) == 3.0);
        CHECK(parameters(1, 0) == 4.0);
        CHECK(parameters(1, 2) == 6.0);
    }

    SECTION("without header") {
        const auto filename =
            tmp_dir.write_file("no_header.csv", "1.5,2.5\n3.5,4.5\n");

        NDArray<double, 2> buffer({2, 2}, 0.0);
        read_csv(filename, buffer.view(), false);

        CHECK(buffer(0, 0) == 1.5);
        CHECK(buffer(0, 1) == 2.5);
        CHECK(buffer(1, 0) == 3.5);
        CHECK(buffer(1, 1) == 4.5);
    }

    SECTION("too many rows throws") {
        const auto filename =
            tmp_dir.write_file("parameters.csv", dg_parameters_csv);
        NDArray<double, 2> parameters({2, 3}, 0.0);
        CHECK_THROWS_AS(read_csv(filename, parameters.view()),
                        std::runtime_error);
    }

    SECTION("too few rows throws") {
        const auto filename =
            tmp_dir.write_file("parameters.csv", dg_parameters_csv);
        NDArray<double, 2> parameters({4, 3}, 0.0);
        CHECK_THROWS_AS(read_csv(filename, parameters.view()),
                        std::runtime_error);
    }

    SECTION("too many columns throws") {
        const auto filename =
            tmp_dir.write_file("parameters.csv", dg_parameters_csv);
        NDArray<double, 2> parameters({3, 2}, 0.0);
        CHECK_THROWS_AS(read_csv(filename, parameters.view()),
                        std::runtime_error);
    }

    SECTION("too few columns throws") {
        const auto filename =
            tmp_dir.write_file("parameters.csv", dg_parameters_csv);
        NDArray<double, 2> parameters({3, 4}, 0.0);
        CHECK_THROWS_AS(read_csv(filename, parameters.view()),
                        std::runtime_error);
    }

    SECTION("non numeric value throws") {
        const auto filename = tmp_dir.write_file(
            "invalid.csv", "center,conversion,offset\n1.0,abc,3.0\n");
        NDArray<double, 2> parameters({1, 3}, 0.0);
        CHECK_THROWS(read_csv(filename, parameters.view()));
    }

    SECTION("non existing file throws") {
        NDArray<double, 2> parameters({3, 3}, 0.0);
        CHECK_THROWS_AS(
            read_csv(tmp_dir.path() / "does_not_exist.csv", parameters.view()),
            std::runtime_error);
    }
}

TEST_CASE("read csv into newly allocated array", "[utils][csv]") {

    TemporaryDirectory tmp_dir("angcal-csv");

    SECTION("number of rows is determined from file") {
        const auto filename =
            tmp_dir.write_file("parameters.csv", dg_parameters_csv);

        auto parameters = read_csv<double, 3>(filename);

        REQUIRE(parameters.shape(0) == 3);
        REQUIRE(parameters.shape(1) == 3);
        CHECK(parameters(0, 0) == Catch::Approx(643.809720708243));
        CHECK(parameters(1, 1) == Catch::Approx(6.57694036246975e-05));
        CHECK(parameters(2, 2) == Catch::Approx(5.805331278507036));
    }

    SECTION("integral values") {
        const auto filename =
            tmp_dir.write_file("ints.csv", "a,b\n1,-2\n30,40\n");

        auto values = read_csv<int, 2>(filename);

        REQUIRE(values.shape(0) == 2);
        CHECK(values(0, 0) == 1);
        CHECK(values(0, 1) == -2);
        CHECK(values(1, 0) == 30);
        CHECK(values(1, 1) == 40);
    }

    SECTION("header only yields zero rows") {
        const auto filename =
            tmp_dir.write_file("header.csv", "center,conversion,offset\n");

        auto parameters = read_csv<double, 3>(filename);

        CHECK(parameters.shape(0) == 0);
        CHECK(parameters.shape(1) == 3);
    }

    SECTION("column mismatch throws") {
        const auto filename =
            tmp_dir.write_file("parameters.csv", dg_parameters_csv);
        CHECK_THROWS_AS((read_csv<double, 2>(filename)), std::runtime_error);
        CHECK_THROWS_AS((read_csv<double, 4>(filename)), std::runtime_error);
    }

    SECTION("non existing file throws") {
        CHECK_THROWS_AS(
            (read_csv<double, 3>(tmp_dir.path() / "does_not_exist.csv")),
            std::runtime_error);
    }
}

// ---------------------------------------------------------------------------
// real data files, requires $ANGCAL_TEST_DATA
// ---------------------------------------------------------------------------

TEST_CASE("read bad channels from real .chans file", "[utils][.files]") {

    const auto filename = test_data_path() / "bc2023_003_RING.chans";
    REQUIRE(std::filesystem::exists(filename));

    constexpr ssize_t num_strips = 48 * MythenDetector::strips_per_module();

    auto bad_channels = read_bad_channels(filename, num_strips);

    CHECK(bad_channels.size() == 61440);
    CHECK(bad_channels(0));
    CHECK(bad_channels(5));
    CHECK_FALSE(bad_channels(6));
    CHECK(bad_channels(8560));
    CHECK(bad_channels(8576));
    CHECK(bad_channels(61437));
    CHECK(std::all_of(bad_channels.begin() + 30720, bad_channels.end(),
                      [](bool element) { return element; }));
}

TEST_CASE("read flatfield from real .raw file", "[utils][.files]") {

    const auto filename =
        test_data_path() /
        "Flatfield_E17p5keV_T12500eV_up_AUGCAL2_Sep2023_open_WS_C_X_X.raw";
    REQUIRE(std::filesystem::exists(filename));

    auto flatfield = read_flatfield(filename, 48);

    CHECK(flatfield.shape(0) == 48 * MythenDetector::strips_per_module());
    CHECK(flatfield(0, 0) == Catch::Approx(-1.0));
    CHECK(flatfield(6, 0) == Catch::Approx(0.76053714144724005));
    CHECK(flatfield(6, 1) == Catch::Approx(3.9644285189937315E-003));
    CHECK(flatfield(21, 0) == Catch::Approx(0.74560040579886455));
    CHECK(flatfield(21, 1) == Catch::Approx(3.9179309274806973E-003));
    CHECK(flatfield(61439, 0) == Catch::Approx(-1.0));
}

TEST_CASE("read DG parameters from real .csv file", "[utils][.files]") {

    const auto filename = test_data_path() / "Angcal_2E_Feb2023_P29.csv";
    REQUIRE(std::filesystem::exists(filename));

    NDArray<double, 2> parameters({48, 3}, 0.0);
    read_csv(filename, parameters.view());

    CHECK(parameters(9, 0) == Catch::Approx(660.342326085818));
    CHECK(parameters(27, 1) == Catch::Approx(-6.581179125986344e-05));
    CHECK(parameters(47, 2) == Catch::Approx(5.805331278507036));

    auto allocated_parameters = read_csv<double, 3>(filename);
    CHECK(allocated_parameters.shape(0) == 48);
    CHECK(allocated_parameters(9, 0) == Catch::Approx(660.342326085818));
}
