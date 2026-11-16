
#include "MythenDetector.hpp"
#include "Parameters.hpp"
#include "logger.hpp"

#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace angcal;

TEST_CASE("check diffraction angle", "[mythen_detector]") {

    DGParameters dg_parameters(48);

    dg_parameters(0, 0) = 643.809720708243;
    dg_parameters(0, 1) = 6.57719636736681e-05;
    dg_parameters(0, 2) = 0.0;
    dg_parameters(5, 0) = 648.772192208828;
    dg_parameters(5, 1) = 6.57918991867043e-05;
    dg_parameters(5, 2) = 25.0071467102665;
    dg_parameters(15, 0) = 623.903939323805;
    dg_parameters(15, 1) = 6.58283077031797e-05;
    dg_parameters(15, 2) = 75.0695128118781;
    dg_parameters(28, 0) = 639.5;
    dg_parameters(28, 1) = -6.581179125986344e-05;
    dg_parameters(28, 2) = 100.80533127850704;
    dg_parameters(47, 0) = 639.5;
    dg_parameters(47, 1) = -6.581179125986344e-05;
    dg_parameters(47, 2) = 5.805331278507036;

    MythenDetector mythen_detector(dg_parameters);

    auto [module_index, local_strip_index] = GENERATE(
        std::make_pair<ssize_t>(0, 0), std::make_pair<ssize_t>(5, 20),
        std::make_pair<ssize_t>(15, 1250), std::make_pair<ssize_t>(47, 1278),
        std::make_pair<ssize_t>(28, 100));

    double diffraction_angle_DG_param =
        mythen_detector.diffraction_angle_from_DG_parameters(
            module_index, 0.0, local_strip_index, 0.0);

    double diffraction_angle_BC_param =
        mythen_detector.diffraction_angle_from_BC_parameters(
            module_index, 0.0, local_strip_index, 0.0);

    // Get the EE parameters for the given module index
    auto [normal_distance, module_center_distance, angle] =
        mythen_detector.get_DGparameters().convert_to_EEParameters(
            module_index);

    double diffraction_angle_EE_param =
        mythen_detector.diffraction_angle_from_EE_parameters(
            module_center_distance, normal_distance, angle, 0.0,
            local_strip_index, 0.0);

    CHECK(diffraction_angle_EE_param ==
          Catch::Approx(diffraction_angle_DG_param));

    CHECK(diffraction_angle_BC_param ==
          Catch::Approx(diffraction_angle_DG_param));
}

/*
TEST_CASE("diffraction angle of clockwise and counterclockwise modules",
          "[mythen_detector]") {

    DGParameters dg_parameters(48);

    dg_parameters(0, 0) = 643.809720708243;
    dg_parameters(0, 1) = 6.57719636736681e-05;
    dg_parameters(0, 2) = 0.0;
    dg_parameters(47, 0) = 639.5;
    dg_parameters(47, 1) = -6.581179125986344e-05;
    dg_parameters(47, 2) = 5.805331278507036;

    MythenDetector mythen_detector(dg_parameters); // dummy dg parameters

    double diffraction_angle_BC_param_clockwise_module =
        mythen_detector.diffraction_angle_from_BC_parameters(0, 0.0, 1);

    double diffraction_angle_BC_param_counterclockwise_module =
        mythen_detector.diffraction_angle_from_BC_parameters(47, 0.0, 1278);

    // TODO: not sure if this should hold
    CHECK(diffraction_angle_BC_param_clockwise_module ==
          Catch::Approx(diffraction_angle_BC_param_counterclockwise_module));
}
*/

TEST_CASE("check displacement correction", "[mythen_detector]") {

    auto [diffraction_angle, corrected_diffraction_angle] =
        GENERATE(std::make_pair(0.0, 0.0), std::make_pair(90.0, 90.0),
                 std::make_pair(-90.0, -90.0), std::make_pair(150.0, 150.0),
                 std::make_pair(-150.0, -150.0), std::make_pair(-75.0, -75.0),
                 std::make_pair(75.0, 75.0), std::make_pair(-180.0, -180.0),
                 std::make_pair(180.0, 180.0));

    MythenDetector mythen_detector(DGParameters(10)); // dummy dg parameters

    mythen_detector.sample_y_displacement(0.0); // no displacement
    mythen_detector.sample_x_displacement(0.0); // no displacement

    double sample_displacement_correction =
        mythen_detector.sample_displacement_correction(diffraction_angle);

    CHECK(sample_displacement_correction ==
          Catch::Approx(corrected_diffraction_angle));
}
