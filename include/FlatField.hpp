
/**
 * stores flatfield for angle calibration
 */

#pragma once
#include <cmath>
#include <cstdint>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>

#include "CustomFiles.hpp"
#include "MythenDetectorSpecifications.hpp"
#include "aare/File.hpp"
#include "aare/NDArray.hpp"
#include "helpers/FileInterface.hpp"

using namespace aare;
namespace angcal {

/*
template <class CustomFile> struct custom_file_compatibility {
    static constexpr bool value =
        std::is_base_of<DetectorFileInterface, CustomFile>::value &&
        std::is_constructible<CustomFile, std::string, ssize_t, ssize_t,
                              std::string>::value;
};

template <class CustomFile>
constexpr bool custom_file_compatibility_v =
    custom_file_compatibility<CustomFile>::value;
*/

// TODO change to uint32_t
class FlatField {

  public:
    /**
     * @tparam custom_file_ptr: Fileclass that inherits from
     * DetectorFileInterface class needs to overload read_frame and constructor
     * needs to have the signature Fileclass(std::string filename, ssize_t rows,
     * ssize_t cols, std::string reading_mode) rows, and cols define the
     * intended image size to read
     */
    FlatField(std::shared_ptr<MythenDetectorSpecifications> mythen_detector_,
              std::optional<std::shared_ptr<SimpleFileInterface>>
                  custom_file_ptr = std::nullopt)
        : mythen_detector(mythen_detector_) {

        if (custom_file_ptr.has_value()) {
            m_custom_detector_file_ptr = custom_file_ptr.value();
        }

        // TODO is it double or uint32_t
        flat_field = NDArray<uint32_t, 1>(
            std::array<ssize_t, 1>{mythen_detector->num_strips()}, 0);
    }

    /**
     * @brief sums up the photon counts for multiple acquisitions
     * @param file_path: path to filesystem - the filesystem should contain
     * multiple acqisitions for different detector angles acquired using
     * slsReceiver and mythen3Detector
     */
    // TODO: adjust for different counters!!!!
    void create_flatfield_from_rawfilesystem(
        const std::filesystem::path &file_path) {

        try {
            for (const auto &file_in_path :
                 std::filesystem::recursive_directory_iterator(file_path)) {
                if (file_in_path.is_regular_file()) {
                    std::string filename =
                        file_in_path.path().filename().string();
                    if (filename.find("master") != std::string::npos) {
                        File f(file_in_path);
                        auto frames = f.read_n(f.total_frames());
                        for (const auto &frame : frames) {
                            if (frame.rows() * frame.cols() !=
                                mythen_detector->num_strips() *
                                    mythen_detector->num_counters()) {
                                throw std::runtime_error(fmt::format(
                                    "sizes mismatch. Expect a size of "
                                    "{} - frame has a size of {}",
                                    mythen_detector->num_strips() *
                                        mythen_detector->num_counters(),
                                    frame.rows() * frame.cols()));
                            }
                            for (ssize_t row = 0; row < frame.rows(); ++row)
                                for (ssize_t col = 0; col < frame.cols();
                                     ++col) {
                                    flat_field(row * frame.cols() + col) +=
                                        *reinterpret_cast<uint32_t *>(
                                            frame.pixel_ptr(
                                                row,
                                                col)); // TODO inefficient as
                                                       // one has to copy twice
                                                       // into frame and into
                                                       // flat_field
                                }
                        }
                    }
                }
            }
        } catch (const std::filesystem::filesystem_error &e) {
            LOG(TLogLevel::logERROR)
                << "Filesystem error: " << e.what() << '\n';
        } catch (const std::exception &e) {
            LOG(TLogLevel::logERROR) << "Runtime error: " << e.what() << '\n';
        }
    }

    /**
     * @brief sums up the photon counts for multiple acquisitions
     * @param filelist: path to file that stores the file paths to the aquired
     * data the list should contain multiple acquisitions for different detector
     * angles
     * @warning only works if member m_custom_detectot_file_ptr supports reading
     * the format
     */
    void create_flatfield_from_filelist(const std::filesystem::path &filelist) {
        std::ifstream file_filelist(filelist);

        try {
            std::string filename;
            while (std::getline(file_filelist, filename)) {
                m_custom_detector_file_ptr->open(filename);
                Frame frame(mythen_detector->num_strips(), 1,
                            aare::Dtype::TypeIndex::UINT32);

                m_custom_detector_file_ptr->read_into(frame.data(),
                                                      frame.dtype().bytes());

                if (frame.rows() * frame.cols() !=
                    mythen_detector->num_strips()) {
                    throw std::runtime_error(
                        fmt::format("sizes mismatch. Expect a size of "
                                    "{} - frame has a size of {}",
                                    mythen_detector->num_strips(),
                                    frame.rows() * frame.cols()));
                }
                for (ssize_t row = 0; row < frame.rows(); ++row)
                    for (ssize_t col = 0; col < frame.cols(); ++col) {
                        flat_field(row * frame.cols() + col) +=
                            *reinterpret_cast<uint32_t *>(frame.pixel_ptr(
                                row,
                                col)); // TODO inefficient as one has to
                                       // copy twice into frame and into
                                       // flat_field
                    }
            }
            file_filelist.close();
        } catch (const std::exception &e) {
            LOG(TLogLevel::logERROR) << "Error: " << e.what() << '\n';
        }
    }

    /**
     * read flatfield from file
     * @warning only works if member m_custom_detectot_file_ptr supports reading
     * the format
     */
    void read_flatfield_from_file(const std::string &filename) {
        m_custom_detector_file_ptr->open(filename);
        std::cout << "ok im reading into" << std::endl;
        m_custom_detector_file_ptr->read_into(
            reinterpret_cast<std::byte *>(flat_field.data()), 4);
    }

    /**
     * set flatfield from NDArray
     */
    void set_flatfield(const NDArray<uint32_t, 1> &flat_field_) {
        flat_field = flat_field_;
    }

    NDArray<uint32_t, 1> get_flatfield() const { return flat_field; }

    void set_inverse_normalized_flatfield(
        const NDArray<double, 1> &inverse_normalized_flatfield_) {
        inverse_normalized_flat_field = inverse_normalized_flatfield_;
    }

    NDView<double, 1> get_inverse_normalized_flatfield() const {
        return inverse_normalized_flat_field.view();
    }

    NDView<double, 1> get_normalized_flatfield() const {
        return normalized_flat_field.view();
    }

    void inverse_normalized_flatfield() {
        double mean = calculate_mean();

        inverse_normalized_flat_field = NDArray<double, 1>(flat_field.shape());

        for (ssize_t i = 0; i < flat_field.size(); ++i) {
            // TODO maybe collapse if else
            inverse_normalized_flat_field[i] =
                (flat_field[i] <= std::numeric_limits<double>::epsilon()
                     ? 0.0
                     : mean / flat_field[i]);

            if (inverse_normalized_flat_field[i] <
                std::numeric_limits<double>::epsilon())
                mythen_detector->get_bad_channels()[i] = true;
        }
    }

    void normalized_flatfield() {
        double mean = calculate_mean();

        normalized_flat_field = NDArray<double, 1>(flat_field.shape());

        for (ssize_t i = 0; i < flat_field.size(); ++i) {
            normalized_flat_field[i] = flat_field[i] / mean;
            if (normalized_flat_field[i] <
                std::numeric_limits<double>::epsilon())
                mythen_detector->get_bad_channels()[i] = true;
        }
    }

  private:
    double calculate_mean() const {
        auto [sum, count] = std::accumulate(
            flat_field.begin(), flat_field.end(),
            std::make_pair<double, ssize_t>(0.0, 0),
            [](std::pair<double, ssize_t> acc, const auto &element) {
                return element == 0 ? acc
                                    : std::make_pair(acc.first + element,
                                                     acc.second + 1);
            });

        return sum / count;
    }

    NDArray<uint32_t, 1> flat_field; // TODO: should be 2d
    NDArray<double, 1> normalized_flat_field;
    NDArray<double, 1> inverse_normalized_flat_field;
    std::shared_ptr<MythenDetectorSpecifications> mythen_detector;
    std::shared_ptr<SimpleFileInterface> m_custom_detector_file_ptr =
        std::make_shared<CustomMythenFile>();
};
} // namespace angcal