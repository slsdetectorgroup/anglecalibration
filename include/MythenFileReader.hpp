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
#include "aare/RawMasterFile.hpp"
#include "aare/File.hpp"
#include "aare/Frame.hpp"

using namespace aare;
namespace angcal {

struct MythenFrame : public aare::Frame {
    //NDArray<uint32_t, 1> photon_counts{};
    double detector_angle{};
    double incident_intensity{}; // I_0
    std::array<uint8_t, 3> channel_mask{};

    MythenFrame() = delete;
    MythenFrame(aare::Frame &&frame) : aare::Frame(std::move(frame)) {};
    MythenFrame(MythenFrame &&) = default; // important for move

    NDView<uint32_t, 1> photon_counts() {
        std::array<ssize_t, 1> shape = {static_cast<ssize_t>(this->rows()) *
                                        static_cast<ssize_t>(this->cols())};
        uint32_t *data = reinterpret_cast<uint32_t *>(this->data());
        return NDView<uint32_t, 1>(data, shape);
    }

    uint32_t photon_counts(size_t row, size_t col = 0) const {
        return *reinterpret_cast<uint32_t *>(this->pixel_ptr(row, col));
    }
};

class MythenFileReader {

  public:
    MythenFileReader() = default;

    MythenFileReader(const std::filesystem::path &detector_angles,
                     const std::filesystem::path &incident_intensities)
        : detector_angles_filename(detector_angles),
          incident_intensities_filename(incident_intensities) {

        // maybe read directly into an NDArray
        try {
            detector_angles_file.open(detector_angles_filename,
                                      std::ios_base::in);
        } catch (const std::exception& error) {
            LOG(TLogLevel::logERROR) << "Could not open file\n";
        }
        try {
            incident_intensities_file.open(incident_intensities_filename,
                                           std::ios_base::in);
        } catch (const std::exception& error) {
            LOG(TLogLevel::logERROR) << "Could not open file\n";
        }
    };

    ~MythenFileReader() {
        detector_angles_file.close();
        incident_intensities_file.close();
    }

    void read_detector_angle(const size_t acquisition_index,
                             double &detector_angle) {
        detector_angles_file.seekg(
            acquisition_index *
            sizeof(double)); // what if its not stored as a double
        detector_angles_file.read(reinterpret_cast<char *>(&detector_angle),
                                  sizeof(double));
    }

    void read_incident_intensity(const size_t acquisition_index,
                                 uint64_t incident_intensity) {
        detector_angles_file.seekg(
            acquisition_index *
            sizeof(uint64_t)); // what if its not stored as a double
        detector_angles_file.read(reinterpret_cast<char *>(&incident_intensity),
                                  sizeof(uint64_t));
    }

    MythenFrame read_frame(const std::filesystem::path &filename) {

        aare::File file(filename); // expects raw file

        // extract acquisition_index

        MythenFrame frame = static_cast<MythenFrame>(file.read_frame());

        size_t acquisition_index = get_acquisition_index(filename);

        read_incident_intensity(acquisition_index, frame.incident_intensity);
        read_detector_angle(acquisition_index, frame.detector_angle);

        return frame;
    }

    size_t get_acquisition_index(const std::filesystem::path &filename) {
        std::string base_name = filename.stem();
        try {
            auto pos = base_name.rfind('_');
            return static_cast<size_t>(std::stoi(base_name.substr(pos + 1)));
        } catch (const std::invalid_argument &e) {
            throw std::runtime_error(LOCATION +
                                     "Could not parse acquisition index");
        }
    }

  private:
    std::string detector_angles_filename{};
    std::string incident_intensities_filename{};
    // NDArray<double> detector_angles{};
    // NDArray<uint32_t> incident_intensities{}; // How many bytes do we need?
    std::ifstream detector_angles_file;
    std::ifstream incident_intensities_file;
};

/** minimal version for a mythen file reader */
class CustomMythenFileReader : public HDF5FileReader, public MythenFileReader {

  public:
    CustomMythenFileReader() = default;

    CustomMythenFileReader(
        const std::filesystem::path &incident_intensity_filename)
        : MythenFileReader(incident_intensity_filename, "") {};

    /**
     * @brief read acquisition file from hdf5
     * @return MythenFrame storing the photon counts, channel mask,
     * detector position [degrees]
     */
    MythenFrame read_frame(const std::string &file_name) {

        open_file(file_name);

        auto dataset_photon_count =
            get_dataset("/entry/instrument/detector/data");

        /*
        frame.photon_counts =
            dataset_photon_count.store_as_ndarray<uint32_t, 1>();
        */

        MythenFrame frame =
            static_cast<MythenFrame>(dataset_photon_count.store_as_frame());

        auto dataset_detector_angle =
            get_dataset("/entry/instrument/NDAttributes/DetectorAngle");

        dataset_detector_angle.read_into_buffer(
            reinterpret_cast<std::byte *>(&frame.detector_angle));

        size_t acquisition_index = get_acquisition_index(file_name);

        MythenFileReader::read_incident_intensity(acquisition_index,
                                                  frame.incident_intensity);

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

    void read_incident_intensity(uint64_t incident_intensity) {
        auto dataset_incident_intensity =
            get_dataset("/entry/instrument/NDAttributes/I0");

        dataset_incident_intensity.read_into_buffer(
            reinterpret_cast<std::byte *>(&incident_intensity));
    }
};

/**
* @brief store as array
* @param filename full path to file, file needs to have one value per line
*/
template <typename T>
NDArray<T> store_as_array(const std::filesystem::path &filename) {

    std::ifstream file;
    file.open(filename); //do i need to open with std::ios::ate?
    
    size_t lineCount = std::count(std::istreambuf_iterator<char>(file),
                                  std::istreambuf_iterator<char>(), '\n') + 1; //TODO: buggy, think about other filesystem - raw, binary file easier

    NDArray<T, 1> array{std::array<ssize_t,1>{lineCount}};

    size_t index = 0;
    std::string line{};

    try {
        while (std::getline(file, line)) {
            array(index) = static_cast<T>(line);
            ++index; 
        }
    } catch (const std::exception &e) {
        LOG(TLogLevel::logERROR) << e.what();
    }

    file.close();

    return array;
}

} // namespace angcal