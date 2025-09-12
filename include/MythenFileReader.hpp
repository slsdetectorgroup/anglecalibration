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
    // double reference_intensity{}; not needed
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

        //++frame.photon_counts; // Why though?

        auto dataset_detector_angle =
            get_dataset("/entry/instrument/NDAttributes/DetectorAngle");

        dataset_detector_angle.read_into_buffer(
            reinterpret_cast<std::byte *>(&frame.detector_angle));

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

} // namespace angcal