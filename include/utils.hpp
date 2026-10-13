#pragma once
#include "AngleCalibration.hpp"
#include "MythenDetector.hpp"
#include "PlotHelpers.hpp"
#include "logger.hpp"

inline auto data_path() {
    if (const char *env_p = std::getenv("ANGCAL_TEST_DATA")) {
        return std::filesystem::path(env_p);
    } else {
        throw std::runtime_error(
            "Path to test data: $ANGCAL_TEST_DATA not set");
    }
}

namespace angcal {

inline std::pair<double, double>
get_detector_range(std::shared_ptr<MythenFileReader> mythen_file_reader,
                   const std::vector<std::string> &file_list) {

    // Files are ordered by angle already
    double min_detector_angle =
        mythen_file_reader->read_detector_angle(file_list.front());
    double max_detector_angle =
        mythen_file_reader->read_detector_angle(file_list.back());

    return std::make_pair(max_detector_angle, min_detector_angle);
}

inline void
get_angle_range_of_module(std::shared_ptr<MythenFileReader> mythen_file_reader,
                          const std::vector<std::string> &file_list,
                          const MythenDetector &mythen_detector,
                          const size_t module_index) {

    constexpr size_t last_strip = MythenDetector::strips_per_module() - 1;

    double detector_angle =
        mythen_file_reader->read_detector_angle(file_list[0]);
    double max_left_module_angle =
        mythen_detector.diffraction_angle_from_DG_parameters(
            module_index, detector_angle, 0, 0.0);
    double max_right_module_angle =
        mythen_detector.diffraction_angle_from_DG_parameters(
            module_index, detector_angle, last_strip, +1.0);

    double min_left_module_angle = max_left_module_angle;
    double min_right_module_angle = max_right_module_angle;

    for (auto &frame : file_list) {
        detector_angle = mythen_file_reader->read_detector_angle(frame);
        double left_module_angle =
            mythen_detector.diffraction_angle_from_DG_parameters(
                module_index, detector_angle, 0, 0.0);
        double right_module_angle =
            mythen_detector.diffraction_angle_from_DG_parameters(
                module_index, detector_angle, last_strip, +1.0);

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
        "minimum module ROI [degrees]: {},{}, maximum module ROI [degrees]: "
        "{},{} for module "
        "{}",
        min_left_module_angle, min_right_module_angle, max_left_module_angle,
        max_right_module_angle, module_index);
}

inline void
get_module_angle_ranges(std::shared_ptr<MythenFileReader> mythen_file_reader,
                        const std::vector<std::string> &file_list,
                        const MythenDetector &mythen_detector) {

    for (size_t module_idx = 0; module_idx < mythen_detector.num_modules();
         ++module_idx) {
        get_angle_range_of_module(mythen_file_reader, file_list,
                                  mythen_detector, module_idx);
    }
}

/**
 * @brief select base peak by performing angular conversion for a specific
 * module and clicking through frames where module is within range of [6°, 33°]
 * @param anglecalibration shared pointer to AngleCalibration class for
 * accessing calibration specific parameters and functions
 * @param mythen_file_reader shared pointer to MythenFileReader for reading
 * mythen acquisition files
 * @param filelist vector of acquisition files to select base peak from
 * @param module_index optional module index if only one module contributed to
 * diffraction pattern, default is 0.
 */
inline void
select_base_peak(std::shared_ptr<MythenDetector> mythen_detector_ptr,
                 std::shared_ptr<AngularConversion> angularconversion,
                 std::shared_ptr<MythenFileReader> mythen_file_reader,
                 const std::vector<std::string> &filelist,
                 const size_t module_index = 0) {

    PlotHelper plotter(mythen_detector_ptr,
                       angularconversion->get_angular_range(),
                       angularconversion->get_histogram_bin_width());
    plotter.overwrite_plot();

    auto detector_angle_range = std::make_pair(
        6.0, 33.0); // TODO: is this absolute to beam source or relative?

    double left_module_strip_angle =
        mythen_detector_ptr->diffraction_angle_from_DG_parameters(module_index,
                                                                  0.0, 0, 0.0);
    LOG(TLogLevel::logDEBUG)
        << fmt::format("strip angle: {}", left_module_strip_angle);

    LOG(TLogLevel::logINFO) << fmt::format(
        "adjusted detector angle range for module {}: {} to {}", module_index,
        detector_angle_range.first, detector_angle_range.second);

    for (const auto &file : filelist) {

        LOG(TLogLevel::logDEBUG) << fmt::format("reading file {}", file);
        double detector_angle = mythen_file_reader->read_detector_angle(file);

        // -5.0, 5.0 adjust for movement
        if (detector_angle + left_module_strip_angle - 5.0 >
                detector_angle_range.first &&
            detector_angle + left_module_strip_angle + 5.0 <
                detector_angle_range.second) {

            LOG(TLogLevel::logINFO)
                << fmt::format("plotting file {} with detector angle {}", file,
                               detector_angle);

            auto frame = mythen_file_reader->read_frame(
                file); // read frame to get photon counts for
            // plot module
            auto module_redistributed_to_fixed_angle_bins =
                angularconversion->convert(
                    frame.photon_counts(), detector_angle,
                    frame.incident_intensity,
                    frame.exposure_time); // always plots diffraction pattern
                                          // independant of detector location.

            plotter.plot_diffraction_pattern(
                module_redistributed_to_fixed_angle_bins.view(), module_index,
                detector_angle);
        }
    }
}

/**
 * @brief reads a comma separated file into a 2D buffer (rows x cols)
 * e.g. DG parameters with header "center,conversion,offset" and one row per
 * module
 * @param filename path to csv file
 * @param buffer view of buffer to fill, shape must match the csv data
 * @param skip_header if true the first line is skipped
 */
inline void read_csv(const std::filesystem::path &filename,
                     aare::NDView<double, 2> buffer,
                     const bool skip_header = true) {

    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename.string());
    }

    std::string line;
    if (skip_header) {
        std::getline(file, line);
    }

    ssize_t row = 0;
    while (std::getline(file, line)) {
        if (line.empty() || line == "\r")
            continue;

        if (row >= buffer.shape(0)) {
            throw std::runtime_error("More rows in " + filename.string() +
                                     " than expected (" +
                                     std::to_string(buffer.shape(0)) + ")");
        }

        std::stringstream line_stream(line);
        std::string value;
        ssize_t col = 0;

        while (std::getline(line_stream, value, ',')) {
            if (col >= buffer.shape(1)) {
                throw std::runtime_error(
                    "More columns in row " + std::to_string(row) + " of " +
                    filename.string() + " than expected (" +
                    std::to_string(buffer.shape(1)) + ")");
            }
            buffer(row, col) = std::stod(value);
            ++col;
        }

        if (col != buffer.shape(1)) {
            throw std::runtime_error(
                "Expected " + std::to_string(buffer.shape(1)) +
                " columns in row " + std::to_string(row) + " of " +
                filename.string() + " but got " + std::to_string(col));
        }
        ++row;
    }

    if (row != buffer.shape(0)) {
        throw std::runtime_error("Expected " + std::to_string(buffer.shape(0)) +
                                 " rows in " + filename.string() + " but got " +
                                 std::to_string(row));
    }
}

/**
 * @brief reads flatfield file (columns: strip_index, flatfield, error, ...)
 * @return array of shape (num_strips, 2) with flatfield value and error per
 * strip
 */
inline NDArray<double, 2> read_flatfield(const std::filesystem::path &filename,
                                         const ssize_t num_modules) {

    if (!std::filesystem::exists(filename)) {
        throw std::runtime_error(
            fmt::format("File {} does not exist\n", filename.string()));
    }

    const ssize_t num_strips =
        num_modules * MythenDetector::strips_per_module();

    NDArray<double, 2> flatfield({num_strips, 2}, 0.0);

    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename.string());
    }

    ssize_t strip_index;
    double flatfield_value, flatfield_error, val_1, val_2;
    std::string line;

    std::getline(file, line); // skip header
    while (file >> strip_index >> flatfield_value >> flatfield_error >> val_1 >>
           val_2) {
        if (strip_index < 0 || strip_index >= num_strips) {
            throw std::runtime_error(
                fmt::format("strip index {} in {} out of range (num strips {})",
                            strip_index, filename.string(), num_strips));
        }
        flatfield(strip_index, 0) = flatfield_value;
        flatfield(strip_index, 1) = flatfield_error;
    }

    return flatfield;
}

inline NDArray<bool, 1> read_bad_channels(const std::filesystem::path &filename,
                                          const ssize_t num_strips) {

    if (!std::filesystem::exists(filename)) {
        throw std::runtime_error(
            fmt::format("File {} does not exist\n", filename.string()));
    }

    NDArray<bool, 1> bad_channels({num_strips});

    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename.string());
    }

    const size_t data_type_bytes = sizeof(bool);

    auto buffer = bad_channels.buffer();

    std::string line;

    size_t index = 0;
    try {
        while (std::getline(file, line)) {

            std::size_t pos = line.find("-");

            if (pos == std::string::npos) {
                index = std::stoi(line) * data_type_bytes;
                buffer[index] = static_cast<std::byte>(true);
            } else {
                size_t line_size = line.size();
                for (int i = std::stoi(line.substr(0, pos));
                     i <= std::stoi(line.substr(pos + 1, line_size - pos));
                     ++i) {
                    index = i * data_type_bytes;
                    buffer[index] = static_cast<std::byte>(true);
                }
            }
        }
    } catch (const std::exception &e) {
        LOG(TLogLevel::logERROR) << e.what();
    }

    file.close();

    return bad_channels;
}

} // namespace angcal