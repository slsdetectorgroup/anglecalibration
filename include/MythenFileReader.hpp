/************************************************
 * @file MythenFileReader.hpp
 * @short minimal file reader to read mythen files
 ***********************************************/

#pragma once

#include <bitset>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

#include "aare/Frame.hpp"
#include "aare/NDArray.hpp"
#include "aare/RawFile.hpp"
#include "aare/RawMasterFile.hpp"
#include "helpers/Hdf5FileReader.hpp"

using namespace aare;
namespace angcal {

/**
 * @brief converts channel mask to bit mask e.g. 4 -> {1 0 0} only channel 0 is
 * enabled
 */
inline std::array<uint8_t, 3>
counter_mask_to_channel_indices(uint8_t counter_mask) {
    std::bitset<3> binary_channel_numbers(counter_mask); // 1 0 0

    // binary_channel_numbers.flip(); // TODO not sure where most
    // significant
    //  bit is ask Anna again

    return std::array<uint8_t, 3>{binary_channel_numbers[0],
                                  binary_channel_numbers[1],
                                  binary_channel_numbers[2]};
}

struct EpicsMythenHeader {
    double detector_angle{};
    uint64_t incident_intensity{}; // I_0
    /// @brief exposure time in seconds
    double exposure_time{}; // in seconds
    /// @brief counter mask given as binary representation (counter 0 is most
    /// significant bit)
    std::array<uint8_t, 3> channel_mask{};
};

struct MythenFrame : public aare::Frame {
    /// @brief detector angle in degrees
    double detector_angle{};
    uint64_t incident_intensity{}; // I_0
    /// @brief exposure time in seconds
    double exposure_time{}; // in seconds

    /// @brief counter mask given as binary representation (counter 0 is most
    /// significant bit)
    std::array<uint8_t, 3> channel_mask{};

    MythenFrame() = delete;
    MythenFrame(aare::Frame &&frame) : aare::Frame(std::move(frame)){};
    MythenFrame(MythenFrame &&) = default; // important for move

    /**
     * @brief photon counts for each channel
     * stored as 1d array
     */
    NDView<uint32_t, 1> photon_counts() {
        std::array<ssize_t, 1> shape = {static_cast<ssize_t>(this->rows()) *
                                        static_cast<ssize_t>(this->cols())};
        uint32_t *data = reinterpret_cast<uint32_t *>(this->data());
        return NDView<uint32_t, 1>(data, shape);
    }

    /**
     * @brief photon counts at specfic channel
     * @param col channel index
     * @param row row index (default: 0)
     */
    uint32_t photon_counts(size_t col, size_t row = 0) const {
        return *reinterpret_cast<uint32_t *>(this->pixel_ptr(row, col));
    }

    /**
     * @brief size of frame (given in pixels)
     */
    ssize_t size() const { return static_cast<ssize_t>(Frame::size()); }
};

/**
 * @brief abstract base class for MythenFileReader
 */
class MythenFileReader {
  public:
    MythenFileReader() = default;
    virtual ~MythenFileReader() = default;

    virtual MythenFrame
    read_frame(const std::filesystem::path &file_name) = 0; // virtual function

    virtual double read_detector_angle(
        const std::filesystem::path &file_name) = 0; // virtual function
};

/**
 * @brief MythenFileReader to read raw files following format from
 * slsDetectorPackage
 */
class RawMythenFileReader : public MythenFileReader {

  public:
    // RawMythenFileReader() = default;

    /**
     * @brief Constructor for MythenFileReader class
     * @param detector_angles path to file storing detector angles of the
     * acquisitions (need to be stored as binary file of doubles)
     * @param incident_intensities path to file storing incident intensities of
     * the acquisitions (need to be stored as binary file of uint64_t)
     */
    RawMythenFileReader(const std::filesystem::path &detector_angles,
                        const std::filesystem::path &incident_intensities)
        : detector_angles_filename(detector_angles),
          incident_intensities_filename(incident_intensities) {

        // maybe read directly into an NDArray
        try {
            detector_angles_file.open(detector_angles_filename,
                                      std::ios_base::in);
            if (!detector_angles_file.is_open())
                throw std::runtime_error(fmt::format("Could not open file {}\n",
                                                     detector_angles_filename));
        } catch (const std::exception &error) {
            LOG(TLogLevel::logERROR) << error.what();
        }
        try {
            incident_intensities_file.open(incident_intensities_filename,
                                           std::ios_base::in);
            if (!incident_intensities_file.is_open())
                throw std::runtime_error(fmt::format(
                    "Could not open file {}\n", incident_intensities_filename));
        } catch (const std::exception &error) {
            LOG(TLogLevel::logERROR) << error.what();
        }
    };

    ~RawMythenFileReader() {
        detector_angles_file.close();
        incident_intensities_file.close();
    }

    /**
     * @brief read Mythenframe
     */
    MythenFrame read_frame(const std::filesystem::path &file_name) override {

        aare::RawFile file(file_name); // expects raw file

        MythenFrame frame = static_cast<MythenFrame>(file.read_frame());

        auto counter_mask = file.master().counter_mask();

        if (counter_mask.value())
            frame.channel_mask =
                counter_mask_to_channel_indices(counter_mask.value());

        auto exposure_time = file.master().exptime();
        if (exposure_time.has_value()) {
            frame.exposure_time =
                std::chrono::duration<double>(exposure_time.value()).count();
        } else {
            LOG(TLogLevel::logWARNING)
                << "Exposure time not found in master file, setting to 0.";
            frame.exposure_time = 0.0;
        }

        // extract acquisition_index
        size_t acquisition_index = get_acquisition_index(file_name);

        read_incident_intensity(acquisition_index, frame.incident_intensity);
        read_detector_angle(acquisition_index, frame.detector_angle);

        return frame;
    }

    double
    read_detector_angle(const std::filesystem::path &file_name) override {
        size_t acquisition_index = get_acquisition_index(file_name);
        double detector_angle{};
        read_detector_angle(acquisition_index, detector_angle);
        return detector_angle;
    }

  private:
    /**
     * @brief reads the detector angle from file
     */
    void read_detector_angle(const size_t acquisition_index,
                             double &detector_angle) {
        detector_angles_file.seekg(
            acquisition_index *
            sizeof(double)); // what if its not stored as a double
        detector_angles_file.read(reinterpret_cast<char *>(&detector_angle),
                                  sizeof(double));
    }

    /**
     * @brief reads the incident intensity from file
     */
    void read_incident_intensity(const size_t acquisition_index,
                                 uint64_t &incident_intensity) {
        // TODO: add a contract check that it is a binary file of uint64_t
        incident_intensities_file.seekg(
            acquisition_index *
            sizeof(uint64_t)); // what if its not stored as a double
        incident_intensities_file.read(
            reinterpret_cast<char *>(&incident_intensity), sizeof(uint64_t));
    }

    /**
     * @brief gets the acquisition index from the filename
     */
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

    std::string detector_angles_filename{};
    std::string incident_intensities_filename{};
    // NDArray<double> detector_angles{};
    // NDArray<uint32_t> incident_intensities{}; // How many bytes do we need?
    std::ifstream detector_angles_file{};
    std::ifstream incident_intensities_file{};
};

/**
 * @brief class to read a series of MythenFiles following Epics file format
 * (e.g. base_name_0000.h5, base_name_0001.h5, ...). Frames are indexed
 * consecutively over all files of the series, when all frames of a file are
 * read the next file is opened.
 */
class EpicsMythenFileReader : public HDF5FileReader {
  public:
    EpicsMythenFileReader(const std::filesystem::path &fname,
                          std::optional<std::filesystem::path>
                              incident_intensities_filename = std::nullopt);

    MythenFrame read_frame();                   // current frame, advances
    MythenFrame read_frame(size_t frame_index); // seek + read_frame()
    std::vector<MythenFrame> read_n(size_t n_frames);
    void read_into(std::byte *buffer); // photon counts only, advances
    void read_into(std::byte *buffer, size_t n_frames);

    double
    read_detector_angle(size_t frame_index) const; // does NOT move cursor
    uint64_t
    read_incident_intensity(size_t frame_index); // does NOT move cursor

    void seek(size_t frame_index);
    size_t tell() const;
    size_t total_frames() const;
    size_t n_files() const;
    size_t pixels_per_frame() const;
    size_t bytes_per_frame() const;
    std::filesystem::path current_file() const;
    std::string base_name() const;

  private:
    void parse_fname(const std::filesystem::path &fname);
    void scan_files();
    void open_file(size_t file_index);
    std::filesystem::path fpath(size_t file_index) const;

    /// @brief index of the file in the series storing frame_index
    size_t get_file_index(size_t frame_index) const;

    /// @brief index of frame_index within its file
    size_t frame_index_in_file(size_t frame_index) const;

    /**
     * @brief reads the value of an NDAttributes dataset for frame_index
     * @param dataset_name name of dataset in /entry/instrument/NDAttributes/
     */
    template <typename T>
    T read_attribute(const std::string &dataset_name, size_t frame_index) const;

    std::filesystem::path m_path{};
    std::string m_base_name{};
    size_t m_max_digits_file_index{};
    size_t m_offset{0};
    size_t m_current_file_index{0};
    size_t m_current_frame_index{0};
    std::vector<size_t> m_last_frame_in_file{};
    size_t m_total_frames{0};
    size_t m_cols{};                        ///< number of strips per frame
    const std::type_info *m_dtype{nullptr}; ///< data type of photon counts
    std::ifstream incident_intensities_file{};
};

/**
 * @brief store as array
 * @param filename full path to file, file needs to have one value per line
 */
template <typename T>
NDArray<T> store_as_array(const std::filesystem::path &filename) {

    std::ifstream file;
    file.open(filename); // do i need to open with std::ios::ate?

    size_t lineCount = std::count(std::istreambuf_iterator<char>(file),
                                  std::istreambuf_iterator<char>(), '\n') +
                       1; // TODO: buggy, think about other filesystem - raw,
                          // binary file easier

    NDArray<T, 1> array{std::array<ssize_t, 1>{lineCount}};

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