/************************************************
 * @file MythenFileReader.hpp
 * @short minimal file reader to read mythen files
 ***********************************************/

#pragma once

#include <bitset>
#include <filesystem>
#include <string>

#include "Hdf5FileReader.hpp"
#include "aare/NDArray.hpp"

using namespace aare;
namespace angcal {

struct MythenFrame {
    NDArray<uint32_t, 1> photon_counts{};
    double detector_angle{};
    double incident_intensity{}; // I_0
    std::array<uint8_t, 3> channel_mask{};

    MythenFrame() = default;
    MythenFrame(MythenFrame &&) = default; // important for move
};

/** minimal version for a mythen file reader */
class MythenFileReader : public HDF5FileReader {

  public:
    MythenFileReader() = default;

    /**
     * @brief read acquisition file from hdf5
     * @return MythenFrame storing the photon counts, channel mask,
     * detector position [degrees]
     */
    MythenFrame read_frame(const std::string &file_name) {

        MythenFrame frame;
        open_file(file_name);

        auto dataset_photon_count =
            get_dataset("/entry/instrument/detector/data");

        frame.photon_counts =
            dataset_photon_count.store_as_ndarray<uint32_t, 1>();

        auto dataset_detector_angle =
            get_dataset("/entry/instrument/NDAttributes/DetectorAngle");

        dataset_detector_angle.read_into_buffer(
            reinterpret_cast<std::byte *>(&frame.detector_angle));

        auto dataset_incident_intensity =
            get_dataset("/entry/instrument/NDAttributes/I0");

        dataset_incident_intensity.read_into_buffer(
            reinterpret_cast<std::byte *>(&frame.incident_intensity));

        auto dataset_channel_number =
            get_dataset("/entry/instrument/NDAttributes/CounterMask");

        uint32_t channel_number;

        dataset_channel_number.read_into_buffer(
            reinterpret_cast<std::byte *>(&channel_number));

        std::bitset<3> binary_channel_numbers(channel_number); // 1 0 0

        // binary_channel_numbers.flip(); // TODO not sure where most
        // significant
        //  bit is ask Anna again

        frame.channel_mask = std::array<uint8_t, 3>{binary_channel_numbers[0],
                                                    binary_channel_numbers[1],
                                                    binary_channel_numbers[2]};

        close_file();

        return frame;
    }
};

class MythenFileReader {

  public:
    MythenFileReader() = default;

    MythenFileReader(const std::string &detector_angles,
                     const std::string &incident_intensities)
        : detector_angles_filename(detector_angles),
          incident_intensities_filename(incident_intensities) {

        // maybe read directly into an NDArray
        detector_angles_file.open(detector_angles_filename, std::ios_base::in);
        incident_intensities_file.open(incident_intensities_file, std::ios_base::in);
        f.read(img.buffer(), img.size() * sizeof(T));
        f.close();
    };

    double read_detector_angle(const size_t acquisition_index) {}

  private:
    std::string detector_angles_filename{};
    std::string incident_intensities_filename{};
    // NDArray<double> detector_angles{};
    // NDArray<uint32_t> incident_intensities{}; // How many bytes do we need?
    std::ifstream detector_angles_file;
    std::ifstream incident_intensities_file;

};
   
} // namespace angcal