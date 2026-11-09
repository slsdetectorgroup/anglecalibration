#pragma once
#include <array>
#include <cstddef>
#include <functional>
#include <vector>

#include "aare/NDArray.hpp"

namespace angcal {

/** @brief part of the calibration a step belongs to */
enum class CalibrationPhase {
    /// optimization of module center sample distance and angle between module
    /// center and module normal
    CoupledParameters,
    /// optimization of angle between module center and beam
    Offset
};

/**
 * @brief state of the calibration of one module after an optimization step
 */
struct CalibrationStep {

    size_t module_index{};

    CalibrationPhase phase{};

    /// index of optimization step, '0' denotes the state before the first step
    size_t iteration{};

    /// similarity criterion of the base peaks
    double similarity{};

    /// BC parameters of the module: module center sample distance, angle
    /// between module center and module normal, angle between module center
    /// and beam
    std::array<double, 3> parameters{};

    /// diffraction angle of each bin in base peak ROI [degrees]
    std::vector<double> bin_angles{};

    /// detector angle of the acquisition of each curve [degrees]
    std::vector<double> detector_angles{};

    /// index of the module each curve was measured with
    std::vector<size_t> curve_modules{};

    /// photon counts redistributed to fixed angle width bins in base peak ROI,
    /// shape (curves, bins), NaN where the base peak ROI is not covered
    aare::NDArray<double, 2> base_peak_curves{};
};

/** @brief function called with the state of the calibration after each
 * optimization step */
using CalibrationObserver = std::function<void(const CalibrationStep &)>;

} // namespace angcal
