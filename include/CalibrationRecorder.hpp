#pragma once
#include <H5Cpp.h>
#include <filesystem>
#include <string>

#include "CalibrationStep.hpp"

namespace angcal {

/**
 * @brief records calibration steps to a hdf5 file
 *
 * Each module and phase gets its own group e.g. /module_07/coupled storing:
 * curves [step, curve, bin], parameters [step, 3], similarity [step],
 * bin_angles [bin], detector_angles [curve], curve_modules [curve]
 */
class CalibrationRecorder {

  public:
    /**
     * @param filename path to hdf5 file, an existing file is overwritten
     */
    explicit CalibrationRecorder(const std::filesystem::path &filename);

    /** @brief appends the step to the file */
    void record(const CalibrationStep &step);

    void operator()(const CalibrationStep &step) { record(step); }

    void close();

    /** @brief name of the group the steps of a module and phase are stored
     * in */
    static std::string group_name(const size_t module_index,
                                  const CalibrationPhase phase);

  private:
    void append_step(const CalibrationStep &step);

    H5::H5File m_file;
};

} // namespace angcal
