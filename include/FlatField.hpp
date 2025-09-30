
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
#include "MythenFileReader.hpp"
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
            m_custom_file_ptr = custom_file_ptr.value();
        }

        flat_field = NDArray<double, 2>(
            std::array<ssize_t, 2>{mythen_detector->num_strips()}, 2);
    }

    /**
     * @brief sums up the photon counts for multiple acquisitions
     * @param filelist: path to file that stores the file paths to the aquired
     * data the list should contain multiple acquisitions for different detector
     * angles
     * @warning only works if member m_custom_detectot_file_ptr supports reading
     * the format
     */
    void create_flatfield_from_filelist(const std::filesystem::path &filelist, const std::filesystem::path &incident_intensities) {
        std::ifstream file_filelist(filelist);

        MythenFileReader file_reader(incident_intensities, ""); //TODO needs to be configurable 

        double average_incident_intensity{};
        size_t count = 0; 

        try {
            std::string filename;

            while(std::getline(file_filelist, filename)) {
                MythenFrame frame = file_reader.read_frame(filename);
                
                if (frame.rows() * frame.cols() != mythen_detector->num_strips()) {
                    throw std::runtime_error(
                        fmt::format("sizes mismatch. Expect a size of "
                                    "{} - frame has a size of {}",
                                    mythen_detector->num_strips(),
                                    frame.rows() * frame.cols()));
                }
                average_incident_intensity += frame.incident_intensity + 1;
                ++count; 
                for (ssize_t row = 0; row < frame.rows(); ++row)
                    for (ssize_t col = 0; col < frame.cols(); ++col) {
                        flat_field(row * frame.cols() + col, 0) += (frame.photon_counts(row, col) + 1)/(frame.incident_intensity + 1);
                             // TODO inefficient as one has to
                                        // copy twice into frame and into
                                        // flat_field   
                        //flatfield variance 
                        flat_field(row * frame.cols() + col, 1) += std::pow(frame.photon_counts(row, col) + 1,1)/std::pow(frame.incident_intensity + 1,3)  +  (frame.photon_counts(row, col) + 1)/std::pow(frame.incident_intensity + 1, 2);        
                    }
            }
        } catch (const std::exception &e) {
            LOG(TLogLevel::logERROR) << "Error: " << e.what() << '\n';
        }

        average_incident_intensity/=count; 
        double average_incident_intensity_2 = average_incident_intensity*average_incident_intensity;
        for(ssize_t index = 0; index < flat_field.shape(0); ++index) {
            flat_field(index, 0) *= average_incident_intensity; // the photon_counts should be scaled by one
            flat_field(index, 1) *= average_incident_intensity_2; // the photon_counts should be scaled by one
        }
            
        file_filelist.close();
        
    }

    /**
     * read flatfield from file
     * @warning only works if member m_custom_detectot_file_ptr supports reading
     * the format
     */
    void read_flatfield_from_file(const std::string &filename) {
        m_custom_file_ptr->open(filename);
        m_custom_file_ptr->read_into(
            reinterpret_cast<std::byte *>(flat_field.data()), 8);
    }

    /**
     * set flatfield from NDArray
     */
    void set_flatfield(const NDArray<double, 2> &flat_field_) {
        flat_field = flat_field_;
    }

    NDArray<double, 2> get_flatfield() const { return flat_field; }

    void set_inverse_normalized_flatfield(
        const NDArray<double, 2> &inverse_normalized_flatfield_) {
        inverse_normalized_flat_field = inverse_normalized_flatfield_;
    }

    NDView<double, 2> get_inverse_normalized_flatfield() const {
        return inverse_normalized_flat_field.view();
    }

    void
    set_normalized_flatfield(const NDArray<double, 2> &normalized_flatfield_) {
        normalized_flat_field = normalized_flatfield_;
    }

    NDView<double, 2> get_normalized_flatfield() const {
        return normalized_flat_field.view();
    }

    void calculate_inverse_normalized_flatfield() {
        double mean = calculate_mean();

        inverse_normalized_flat_field = NDArray<double, 2>(flat_field.shape());

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

    void calculate_normalized_flatfield() {
        double mean = calculate_mean();

        normalized_flat_field = NDArray<double, 2>(flat_field.shape());

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

    NDArray<double, 2> flat_field; 
    NDArray<double, 2> normalized_flat_field;
    NDArray<double, 2> inverse_normalized_flat_field;
    std::shared_ptr<MythenDetectorSpecifications> mythen_detector;
    std::shared_ptr<SimpleFileInterface> m_custom_file_ptr =
        std::make_shared<CustomFlatFieldFile>();
};
} // namespace angcal