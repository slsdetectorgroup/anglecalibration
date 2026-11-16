#pragma once
#include "fmt/format.h"
#include <H5Cpp.h>
#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#include "TempFileDirectory.hpp"

/// @brief meta data and photon counts of one dummy acquisition
struct DummyAcquisition {
    std::vector<uint32_t> photon_counts;
    double detector_angle{};
    double izero{};
    double exposure_time{};
    uint32_t counter_mask{};
};

namespace {
template <typename T>
void write_dataset(H5::Group &group, const std::string &name,
                   const H5::PredType &type, const std::vector<hsize_t> &dims,
                   const std::vector<T> &data) {
    H5::DataSpace dataspace(static_cast<int>(dims.size()), dims.data());
    H5::DataSet dataset = group.createDataSet(name, type, dataspace);
    dataset.write(data.data(), type);
}
} // namespace

class TempEpicsMythenFile {

  public:
    TempEpicsMythenFile(const std::string &base_filename,
                        const std::vector<DummyAcquisition> &acquisitions,
                        const size_t frames_per_file)
        : m_base_filename(base_filename) {
        if (frames_per_file == 0) {
            throw std::invalid_argument("frames_per_file must be > 0");
        }

        const size_t n_acquisitions = acquisitions.size();
        const size_t n_files =
            (n_acquisitions + frames_per_file - 1) / frames_per_file;
        for (size_t file = 0; file < n_files; ++file) {
            write_epics_file(
                m_filesystem_dir.path() /
                    fmt::format("{}_{:04}.h5", m_base_filename, file),
                acquisitions.begin() +
                    static_cast<std::ptrdiff_t>(file * frames_per_file),
                acquisitions.begin() + static_cast<std::ptrdiff_t>(std::min(
                                           (file + 1) * frames_per_file,
                                           n_acquisitions)));
        }
    }

    std::filesystem::path get_filename() const {
        return m_filesystem_dir.path() /
               fmt::format("{}_{:04}.h5", m_base_filename, 0);
    }

    std::filesystem::path get_filesystem_dir() const {
        return m_filesystem_dir.path();
    }

  private:
    using AcquisitionIterator = std::vector<DummyAcquisition>::const_iterator;

    /**
     * @brief writes an hdf5 file following the Epics layout expected by
     * EpicsMythenFileReader:
     *   /entry/instrument/detector/data      photon counts (strips) for a
     *                                        single frame or (frames x strips)
     *   /entry/instrument/NDAttributes/...   one value per frame
     * @param acq_begin, acq_end range of acquisitions stored in this file
     */
    static void write_epics_file(const std::filesystem::path &filename,
                                 AcquisitionIterator acq_begin,
                                 AcquisitionIterator acq_end) {
        const auto n_frames =
            static_cast<hsize_t>(std::distance(acq_begin, acq_end));
        const hsize_t n_strips = acq_begin->photon_counts.size();

        H5::H5File file(filename.string(), H5F_ACC_TRUNC);
        H5::Group entry = file.createGroup("/entry");
        H5::Group data = file.createGroup("/entry/data");
        H5::Group instrument = file.createGroup("/entry/instrument");
        H5::Group detector = file.createGroup("/entry/instrument/detector");
        H5::Group attributes =
            file.createGroup("/entry/instrument/NDAttributes");

        std::vector<uint32_t> photon_counts;
        photon_counts.reserve(n_frames * n_strips);
        std::vector<double> detector_angles, izeros, exposure_times;
        std::vector<uint32_t> counter_masks;
        for (auto acquisition = acq_begin; acquisition != acq_end;
             ++acquisition) {
            photon_counts.insert(photon_counts.end(),
                                 acquisition->photon_counts.begin(),
                                 acquisition->photon_counts.end());
            detector_angles.push_back(acquisition->detector_angle);
            izeros.push_back(acquisition->izero);
            exposure_times.push_back(acquisition->exposure_time);
            counter_masks.push_back(acquisition->counter_mask);
        }

        // single frame is stored as 1d dataset, several frames as 2d dataset
        const std::vector<hsize_t> data_dims =
            n_frames == 1 ? std::vector<hsize_t>{n_strips}
                          : std::vector<hsize_t>{n_frames, n_strips};

        write_dataset(detector, "data", H5::PredType::NATIVE_UINT32,
                      data_dims, photon_counts);
        write_dataset(attributes, "DetectorAngle", H5::PredType::IEEE_F64LE,
                      {n_frames}, detector_angles);
        write_dataset(attributes, "Izero", H5::PredType::IEEE_F64LE,
                      {n_frames}, izeros);
        write_dataset(attributes, "AcquireTime", H5::PredType::IEEE_F64LE,
                      {n_frames}, exposure_times);
        write_dataset(attributes, "CounterMask", H5::PredType::NATIVE_UINT32,
                      {n_frames}, counter_masks);
    }

    std::string m_base_filename;
    TemporaryDirectory m_filesystem_dir{"mythen-acquisitions"};
};
