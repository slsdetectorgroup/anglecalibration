#include "MythenFileReader.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include "np_helpers.hpp"

namespace py = pybind11;

using namespace angcal;

void define_MythenFrame_bindings(py::module &m) {
    py::class_<MythenFrame>(m, "MythenFrame")
        .def(py::init<>())
        .def_property_readonly(
            "photon_counts",
            [](MythenFrame &self) {
                auto photon_counts = new NDArray<uint32_t, 1>(
                    self.photon_counts); // TODO: maybe dont want a copy
                                         // actually
                return return_image_data(photon_counts);
            })

        .def_property_readonly(
            "detector_angle",
            [](MythenFrame &self) { return self.detector_angle; })
        .def_property_readonly("channel_mask", [](MythenFrame &self) {
            return self.channel_mask;
        });
}

void define_MythenFileReader_bindings(py::module &m) {
    py::class_<MythenFileReader, std::shared_ptr<MythenFileReader>>(
        m, "MythenFileReader")
        .def(py::init<>())

        .def("read_frame",
             [](MythenFileReader &self, const std::string &file_name) {
                 return self.read_frame(file_name);
             });
}