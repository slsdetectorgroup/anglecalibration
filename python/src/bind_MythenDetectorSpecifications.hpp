#include "MythenDetectorSpecifications.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;

using namespace angcal;

void define_MythenDetectorSpecifications_binding(py::module &m) {

    py::class_<MythenDetectorSpecifications,
               std::shared_ptr<MythenDetectorSpecifications>>(
        m, "MythenDetectorSpecifications", R"(Attributes
----------
min_angle : double 
    Read-only static. Minimum potential detector angle
    (measured as displacement of first strip) [degree]
    (-180.0 deg)
max_angle : double 
    Read-only static. Maximum potential detector angle
    (measured as displacement of first strip) [degree]
    (-180.0 deg) 
bad_channels : numpy.ndarray of bool, shape (n_channels,)
    Expected size: number of channels/strips in the detector.
    Each element is ``True`` if the channel is bad, otherwise ``False``.)
pitch: double 
    Read-only static. Strip/channel width of Mythen detector [mm] 
    (0.05 mm)
strips_per_module: int
    Read-only static. Strips per module for Mythen detector 
    (1280)
num_strips: int
    Read-only. Total number of strips in Mythen detector
max_modules: int 
    Read-only. Number of modules in detector. 
    Default (48)
exposure_time: double 
    Read-only. Exposure time [s]
bloffset: double 
dtt0: double 
)")

        .def(py::init<std::optional<std::shared_ptr<SimpleFileInterface>>>(),
             py::arg("file_interface") = std::nullopt,
             R"(
             Parameters
             ----------
             file_interface : Optional[SimpleFileInterface], default None
                 A file interface to read bad channels file. 
                 If none is provided bad channel file is expected to be a text file where each line stores the channel index of a bad channel. Consecutive bad channels can be stored in one line by seperating the first and last channel index of the bad channel block e.g. bad_channel_index0-bad_channel_index1.)")

        .def(py::init<const size_t, const double, const double, const double,
                      std::optional<std::shared_ptr<SimpleFileInterface>>>(),
             py::arg("max_modules"), py::arg("exposure_time"),
             py::arg("num_counters") = 1, py::arg("bloffset") = 1.532,
             py::arg("file_interface") = std::nullopt,
             R"(
             Parameters
             ----------
             max_modules: int 
                Number of modules in detector (default 48).
             exposure_time: double 
                Exposure time in seconds. 
             num_counters: int
                Number of counters active. 
             file_interface : Optional[SimpleFileInterface], default None
                 A file interface to read bad channels file. If none is provided bad channel file is expected to be a text file where each line stores the channel index of a bad channel. Consecutive bad channels can be stored in one line by seperating the first and last channel index of the bad channel block e.g. bad_channel_index0-bad_channel_index1.)")

        .def(
            "read_bad_channels_from_file",
            [](MythenDetectorSpecifications &self,
               const std::string &filename) {
                self.read_bad_channels_from_file(filename);
            },
            py::arg("filename"), R"(reads bad channels from file)")

        .def_property(
            "unconnected_modules",
            [](MythenDetectorSpecifications &self) {
                auto unconnected_modules =
                    new NDArray<ssize_t, 1>(self.get_unconnected_modules());
                return return_image_data(unconnected_modules);
            },
            [](MythenDetectorSpecifications &self,
               py::buffer unconnected_modules) {
                py::buffer_info info = unconnected_modules.request();
                if (info.ndim != 1 ||
                    info.format !=
                        py::format_descriptor<std::size_t>::format()) {
                    throw std::runtime_error(
                        "Expected 1D buffer of type unsigned int");
                }
                std::vector<ssize_t> temp_vec(info.shape[0]);
                std::memcpy(temp_vec.data(), info.ptr,
                            info.shape[0] *
                                sizeof(ssize_t)); // TODO: std::size_t
                                                  // exists in python?

                self.set_unconnected_modules(temp_vec);
            })

        .def_property(
            "bad_channels",
            [](MythenDetectorSpecifications &self) {
                auto bad_channels =
                    new NDArray<ssize_t, 1>(self.get_bad_channels());
                return return_image_data(bad_channels);
            },
            [](MythenDetectorSpecifications &self,
               py::array_t<bool, py::array::forcecast> bad_channels) {
                py::buffer_info info = bad_channels.request();
                if (info.ndim != 1 ||
                    info.format != py::format_descriptor<bool>::format()) {
                    throw std::runtime_error("Expected 1D buffer of type bool");
                }
                NDView<bool, 1> temp_array_view(
                    reinterpret_cast<bool *>(info.ptr),
                    std::array<ssize_t, 1>{info.shape[0]});
                NDArray temp_array(temp_array_view); // first copy
                self.set_bad_channels(
                    temp_array); // second copy TODO im copying twice
            })

        .def_property_readonly_static(
            "pitch", &MythenDetectorSpecifications::pitch, R"()")

        .def_property_readonly_static(
            "strips_per_module",
            &MythenDetectorSpecifications::strips_per_module)

        .def_property_readonly("max_modules",
                               &MythenDetectorSpecifications::max_modules)

        .def_property_readonly("num_counters",
                               &MythenDetectorSpecifications::num_counters)

        .def_property_readonly("exposure_time",
                               &MythenDetectorSpecifications::exposure_time)

        .def_property_readonly("bloffset",
                               &MythenDetectorSpecifications::bloffset)

        .def_property_readonly("dtt0", &MythenDetectorSpecifications::dtt0)

        .def_property_readonly_static("min_angle",
                                      &MythenDetectorSpecifications::min_angle)

        .def_property_readonly_static("max_angle",
                                      &MythenDetectorSpecifications::max_angle)

        .def_property_readonly("num_strips",
                               &MythenDetectorSpecifications::num_strips);
}