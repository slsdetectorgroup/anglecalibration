#include "CalibrationRecorder.hpp"
#include "aare/defs.hpp"

#include <fmt/format.h>

namespace angcal {

namespace {

// creates a dataset which can grow along the first dimension
H5::DataSet create_extendable_dataset(H5::Group &group, const std::string &name,
                                      std::vector<hsize_t> step_shape) {

    std::vector<hsize_t> dims{0};
    std::vector<hsize_t> max_dims{H5S_UNLIMITED};
    std::vector<hsize_t> chunk_dims{step_shape.empty() ? hsize_t{64}
                                                       : hsize_t{1}};
    dims.insert(dims.end(), step_shape.begin(), step_shape.end());
    max_dims.insert(max_dims.end(), step_shape.begin(), step_shape.end());
    chunk_dims.insert(chunk_dims.end(), step_shape.begin(), step_shape.end());

    H5::DSetCreatPropList properties;
    properties.setChunk(static_cast<int>(chunk_dims.size()), chunk_dims.data());

    return group.createDataSet(
        name, H5::PredType::NATIVE_DOUBLE,
        H5::DataSpace(static_cast<int>(dims.size()), dims.data(),
                      max_dims.data()),
        properties);
}

template <typename T>
void write_dataset(H5::Group &group, const std::string &name,
                   const std::vector<T> &values, const H5::PredType &type) {
    const hsize_t size = values.size();
    group.createDataSet(name, type, H5::DataSpace(1, &size))
        .write(values.data(), type);
}

// appends one step to a dataset created with create_extendable_dataset
void append(H5::DataSet dataset, const double *values) {

    H5::DataSpace space = dataset.getSpace();
    std::vector<hsize_t> dims(space.getSimpleExtentNdims());
    space.getSimpleExtentDims(dims.data());

    std::vector<hsize_t> offset(dims.size(), 0);
    offset[0] = dims[0];

    dims[0] += 1;
    dataset.extend(dims.data());

    std::vector<hsize_t> count = dims;
    count[0] = 1;

    H5::DataSpace file_space = dataset.getSpace();
    file_space.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data());

    H5::DataSpace memory_space(static_cast<int>(count.size()), count.data());

    dataset.write(values, H5::PredType::NATIVE_DOUBLE, memory_space,
                  file_space);
}

// H5::Exception is not derived from std::exception
std::runtime_error to_runtime_error(const H5::Exception &error) {
    return std::runtime_error(fmt::format("hdf5 error in {}: {}",
                                          error.getFuncName(),
                                          error.getDetailMsg()));
}

H5::H5File create_file(const std::filesystem::path &filename) {
    try {
        return H5::H5File(filename.string(), H5F_ACC_TRUNC);
    } catch (const H5::Exception &error) {
        throw to_runtime_error(error);
    }
}

} // namespace

CalibrationRecorder::CalibrationRecorder(const std::filesystem::path &filename)
    : m_file(create_file(filename)) {}

std::string CalibrationRecorder::group_name(const size_t module_index,
                                            const CalibrationPhase phase) {
    return fmt::format("module_{:02}/{}", module_index,
                       phase == CalibrationPhase::CoupledParameters ? "coupled"
                                                                    : "offset");
}

void CalibrationRecorder::record(const CalibrationStep &step) {
    try {
        append_step(step);
    } catch (const H5::Exception &error) {
        throw to_runtime_error(error);
    }
}

void CalibrationRecorder::append_step(const CalibrationStep &step) {

    const std::string name = group_name(step.module_index, step.phase);

    const hsize_t num_curves = step.base_peak_curves.shape(0);
    const hsize_t num_bins = step.base_peak_curves.shape(1);

    // nameExists fails if the parent group does not exist
    const std::string module_group = name.substr(0, name.find('/'));

    if (!m_file.nameExists(module_group) || !m_file.nameExists(name)) {

        H5::LinkCreatPropList create_intermediate_groups;
        create_intermediate_groups.setCreateIntermediateGroup(true);

        H5::Group group(H5Gcreate2(m_file.getId(), name.c_str(),
                                   create_intermediate_groups.getId(),
                                   H5P_DEFAULT, H5P_DEFAULT));

        create_extendable_dataset(group, "curves", {num_curves, num_bins});
        create_extendable_dataset(group, "parameters", {3});
        create_extendable_dataset(group, "similarity", {});

        write_dataset(group, "bin_angles", step.bin_angles,
                      H5::PredType::NATIVE_DOUBLE);
        write_dataset(group, "detector_angles", step.detector_angles,
                      H5::PredType::NATIVE_DOUBLE);
        write_dataset(group, "curve_modules",
                      std::vector<uint64_t>(step.curve_modules.begin(),
                                            step.curve_modules.end()),
                      H5::PredType::NATIVE_UINT64);
    }

    H5::Group group = m_file.openGroup(name);

    H5::DataSet curves = group.openDataSet("curves");
    hsize_t dims[3];
    curves.getSpace().getSimpleExtentDims(dims);
    if (dims[1] != num_curves || dims[2] != num_bins) {
        throw std::runtime_error(
            LOCATION +
            fmt::format("step of {} has {} curves with {} bins but the file "
                        "stores {} curves with {} bins",
                        name, num_curves, num_bins, dims[1], dims[2]));
    }

    append(curves, step.base_peak_curves.data());
    append(group.openDataSet("parameters"), step.parameters.data());
    append(group.openDataSet("similarity"), &step.similarity);

    m_file.flush(H5F_SCOPE_GLOBAL);
}

void CalibrationRecorder::close() { m_file.close(); }

} // namespace angcal
