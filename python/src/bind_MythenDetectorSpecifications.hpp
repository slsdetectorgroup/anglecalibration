#include "MythenDetectorSpecifications.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;

using namespace angcal;

void define_MythenDetectorSpecifications_binding(py::module &m) {

    py::class_<MythenDetectorSpecifications,
               std::shared_ptr<MythenDetectorSpecifications>>(
        m, "MythenDetectorSpecifications")

        .def(py::init<std::optional<std::shared_ptr<SimpleFileInterface>>>(),
             py::arg("file_interface") = std::nullopt)

        .def(py::init<const size_t, const double, const double, const double,
                      std::optional<std::shared_ptr<SimpleFileInterface>>>(),
             py::arg("max_modules"), py::arg("exposure_time"),
             py::arg("num_counters") = 1, py::arg("bloffset") = 1.532,
             py::arg("file_interface") = std::nullopt)

        .def("read_bad_channels_from_file",
             [](MythenDetectorSpecifications &self,
                const std::string &filename) {
                 self.read_bad_channels_from_file(filename);
             })

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
               py::array_t<bool> bad_channels) {
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

        .def_property_readonly_static("pitch",
                                      &MythenDetectorSpecifications::pitch)

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