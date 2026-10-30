#include "MythenFileReader.hpp"

#include "aare/algorithm.hpp"

namespace angcal {

/**
 * @brief Constructor for EpicsMythenFileReader class
 * @param fname path to any hdf5 file of the series storing the acquisitions,
 * the series starts at this file
 * @param incident_intensities_filename optional path to file storing incident
 * intensities of the acquisitions (need to be stored as binary file of
 * uint64_t)
 * @note if no incident intensities file is given the incident intensities
 * written by epics in the hdf5 file are taken. For newer versions this will be
 * garbage.
 */
EpicsMythenFileReader::EpicsMythenFileReader(
    const std::filesystem::path &fname,
    std::optional<std::filesystem::path> incident_intensities_filename) {
    if (!std::filesystem::exists(fname)) {
        throw std::runtime_error(LOCATION +
                                 "File does not exist: " + fname.string());
    }
    parse_fname(fname);
    scan_files();
    open_file(m_current_file_index);
    // TODO: maybe read meta data here instead of in scan_files

    if (incident_intensities_filename) {
        // TODO maybe directly read into NDArray faster than always reading
        // from filesystem
        incident_intensities_file.open(incident_intensities_filename.value(),
                                       std::ios_base::in);
        if (!incident_intensities_file.is_open()) {
            throw std::runtime_error(
                LOCATION + "Could not open file " +
                incident_intensities_filename.value().string());
        }
    }
}

/**
 * @brief read current frame and advance to the next frame
 * @return tuple of header (channel mask, detector position [degrees],
 * incident intensity and exposure time) and MythenFrame storing the photon
 * counts
 */
std::tuple<EpicsMythenHeader, MythenFrame> EpicsMythenFileReader::read_frame() {
    MythenFrame frame(aare::Frame(
        1, static_cast<uint32_t>(m_cols),
        Dtype(*m_dtype))); // TODO: take into account multiple frames

    EpicsMythenHeader header{};
    read_into(frame.data(), &header);

    return {header, std::move(frame)};
}

/**
 * @brief read frame at frame_index and advance to the next frame
 */
std::tuple<EpicsMythenHeader, MythenFrame>
EpicsMythenFileReader::read_frame(size_t frame_index) {
    seek(frame_index);
    return read_frame();
}

/**
 * @brief read the next n_frames frames, opens the next file if all frames of
 * the current file are read
 * @return tuple of headers and frames
 */
std::tuple<std::vector<EpicsMythenHeader>, std::vector<MythenFrame>>
EpicsMythenFileReader::read_n(size_t n_frames) {
    std::vector<EpicsMythenHeader> headers;
    std::vector<MythenFrame> frames;
    headers.reserve(n_frames);
    frames.reserve(n_frames);
    for (size_t i = 0; i < n_frames; ++i) {
        auto [header, frame] = read_frame();
        headers.push_back(header);
        frames.push_back(std::move(frame));
    }
    return {std::move(headers), std::move(frames)};
}

/**
 * @brief read photon counts of current frame into buffer and advance to the
 * next frame
 * @param buffer buffer to read the frame into (caller is responsible for
 * allocating bytes_per_frame() bytes)
 * @param header optional header to fill with the meta data of the frame
 */
void EpicsMythenFileReader::read_into(std::byte *buffer,
                                      EpicsMythenHeader *header) {
    if (m_current_frame_index >= m_total_frames) {
        throw std::runtime_error(LOCATION + "No frames left to read");
    }

    if (header) {
        *header = read_header(m_current_frame_index);
    }

    auto dataset_photon_count = get_dataset("/entry/instrument/detector/data");

    if (dataset_photon_count.get_shape().size() == 1) {
        dataset_photon_count.read_into_buffer(buffer);
    } else {
        dataset_photon_count.read_into_buffer(
            buffer, Subset{{1, m_cols},
                           {frame_index_in_file(m_current_frame_index), 0}});
    }

    // prepare for next read
    // if we are at the end of the file, open the next file
    ++m_current_frame_index;
    if (m_current_frame_index >= m_last_frame_in_file[m_current_file_index] &&
        (m_current_frame_index < m_total_frames)) {
        open_file(m_current_file_index + 1);
    }
}

/**
 * @brief read photon counts of the next n_frames frames into buffer
 * @param buffer buffer to read the frames into (caller is responsible for
 * allocating n_frames * bytes_per_frame() bytes)
 * @param headers optional array of headers to fill with the meta data of the
 * frames (caller is responsible for allocating n_frames headers)
 */
void EpicsMythenFileReader::read_into(std::byte *buffer, size_t n_frames,
                                      EpicsMythenHeader *headers) {
    for (size_t i = 0; i < n_frames; ++i) {
        read_into(buffer + i * bytes_per_frame(),
                  headers ? headers + i : nullptr);
    }
}

EpicsMythenHeader EpicsMythenFileReader::read_header(size_t frame_index) {
    EpicsMythenHeader header{};
    header.detector_angle = read_detector_angle(frame_index);
    header.incident_intensity = read_incident_intensity(frame_index);
    header.exposure_time = read_attribute<double>("AcquireTime", frame_index);
    header.channel_mask = counter_mask_to_channel_indices(
        read_attribute<uint32_t>("CounterMask", frame_index));
    return header;
}

/**
 * @brief read the detector monitor position of frame at frame_index
 * @return detector position [degrees]
 */
double EpicsMythenFileReader::read_detector_angle(size_t frame_index) const {
    return read_attribute<double>("DetectorAngle", frame_index);
}

/**
 * @brief read the incident intensity of frame at frame_index either from
 * seperate file or from the epics hdf5 file
 */
uint64_t EpicsMythenFileReader::read_incident_intensity(size_t frame_index) {
    uint64_t incident_intensity{};
    if (incident_intensities_file.is_open()) {
        // acquisition index is the file number for files storing one frame
        const size_t acquisition_index = m_offset + frame_index;
        // TODO: add a contract check that it is a binary file of uint64_t
        incident_intensities_file.clear();
        incident_intensities_file.seekg(acquisition_index * sizeof(uint64_t));
        incident_intensities_file.read(
            reinterpret_cast<char *>(&incident_intensity), sizeof(uint64_t));
    } else {
        incident_intensity = static_cast<uint64_t>(
            std::lround(read_attribute<double>("Izero", frame_index)));
    }
    return incident_intensity;
}

/**
 * @brief seek to the given frame index (index over all files of the series)
 */
void EpicsMythenFileReader::seek(size_t frame_index) {
    if (frame_index >= m_total_frames) {
        throw std::runtime_error(LOCATION + "Frame index out of range: " +
                                 std::to_string(frame_index));
    }
    m_current_frame_index = frame_index;
    auto file_index = get_file_index(frame_index);

    if (file_index != m_current_file_index)
        open_file(file_index);
}

size_t EpicsMythenFileReader::tell() const { return m_current_frame_index; }

size_t EpicsMythenFileReader::total_frames() const { return m_total_frames; }

size_t EpicsMythenFileReader::n_files() const {
    return m_last_frame_in_file.size();
}

size_t EpicsMythenFileReader::pixels_per_frame() const { return m_cols; }

size_t EpicsMythenFileReader::bytes_per_frame() const {
    return m_cols * Dtype(*m_dtype).bytes();
}

std::filesystem::path EpicsMythenFileReader::current_file() const {
    return fpath(m_current_file_index + m_offset);
}

std::string EpicsMythenFileReader::base_name() const { return m_base_name; }

/**
 * @brief parse the filename to extract the base name, file index and number
 * of digits of the file index
 * @param fname path to the hdf5 file storing the acquisitions
 */
void EpicsMythenFileReader::parse_fname(const std::filesystem::path &fname) {
    m_path = fname.parent_path();
    m_base_name = fname.stem().string();

    // find file index, then remove if from the base name
    auto pos = m_base_name.find_last_of('_');
    if (pos == std::string::npos) {
        throw std::runtime_error(LOCATION + "Could not parse file index of " +
                                 fname.string());
    }
    m_max_digits_file_index = m_base_name.size() - pos - 1;
    try {
        m_offset = std::stoul(m_base_name.substr(pos + 1));
    } catch (const std::invalid_argument &e) {
        throw std::runtime_error(LOCATION + "Could not parse file index of " +
                                 fname.string());
    }
    m_base_name.erase(pos);
}

/**
 * @brief find the files of the series and the number of frames in each file
 */
void EpicsMythenFileReader::scan_files() {
    m_last_frame_in_file.clear();
    size_t file_index = m_offset;
    while (std::filesystem::exists(fpath(file_index))) {
        HDF5FileReader::open_file(fpath(file_index).string());
        auto dataset_photon_count =
            get_dataset("/entry/instrument/detector/data");
        auto shape = dataset_photon_count.get_shape();

        // one frame per file stored as 1d dataset or several frames stored as
        // 2d dataset (frames x strips)
        size_t n_frames{};
        size_t cols{};
        if (shape.size() == 1) {
            n_frames = 1;
            cols = shape[0];
        } else if (shape.size() == 2) {
            n_frames = shape[0];
            cols = shape[1]; // TODO: could be 2d for counters
        } else {
            throw std::runtime_error(LOCATION +
                                     "Expected 1d or 2d dataset in file " +
                                     fpath(file_index).string());
        }

        if (file_index == m_offset) {
            m_cols = cols;
            m_dtype = dataset_photon_count.get_cpp_type();
        }
        // TODO: maybe check counters
        if (cols != m_cols || dataset_photon_count.get_cpp_type() != m_dtype) {
            throw std::runtime_error(
                LOCATION + "Dataset shape or type mismatch in file " +
                fpath(file_index)
                    .string()); // TODO: maybe rename to meta data mismatch
        }

        m_last_frame_in_file.push_back(n_frames);
        ++file_index;
    }

    if (m_last_frame_in_file.empty()) {
        throw std::runtime_error(LOCATION + "No files found matching " +
                                 fpath(m_offset).string());
    }

    // find where we need to open the next file and total number of frames
    m_last_frame_in_file = cumsum(m_last_frame_in_file);
    m_total_frames = m_last_frame_in_file.back();
}

/**
 * @brief open file of the series
 * @param file_index index of the file in the series (starting at 0)
 */
void EpicsMythenFileReader::open_file(size_t file_index) {
    close_file(); // close current file if open
    HDF5FileReader::open_file(fpath(file_index + m_offset).string());
    m_current_file_index = file_index;
}

std::filesystem::path EpicsMythenFileReader::fpath(size_t file_index) const {
    return m_path /
           (m_base_name + "_" +
            fmt::format("{:0{}d}", file_index, m_max_digits_file_index) +
            ".h5");
}

size_t EpicsMythenFileReader::get_file_index(size_t frame_index) const {
    return first_larger(m_last_frame_in_file, frame_index);
}

size_t EpicsMythenFileReader::frame_index_in_file(size_t frame_index) const {
    const size_t file_index = get_file_index(frame_index);
    return file_index ? frame_index - m_last_frame_in_file[file_index - 1]
                      : frame_index;
}

template <typename T>
T EpicsMythenFileReader::read_attribute(const std::string &dataset_name,
                                        size_t frame_index) const {
    if (frame_index >= m_total_frames) {
        throw std::runtime_error(LOCATION + "Frame index out of range: " +
                                 std::to_string(frame_index));
    }

    const Subset subset{{1}, {frame_index_in_file(frame_index)}};

    T value{};

    const size_t file_index = get_file_index(frame_index);

    if (file_index == m_current_file_index) {
        get_dataset("/entry/instrument/NDAttributes/" + dataset_name)
            .read_into_buffer(reinterpret_cast<std::byte *>(&value), subset);
    } else {
        // read from file without moving the current frame
        HDF5FileReader temp_file_reader;

        temp_file_reader.open_file(fpath(file_index + m_offset).string());
        temp_file_reader
            .get_dataset("/entry/instrument/NDAttributes/" + dataset_name)
            .read_into_buffer(reinterpret_cast<std::byte *>(&value), subset);
    }

    return value;
}

} // namespace angcal
