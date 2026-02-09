#include "MythenDetectorSpecifications.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;

using namespace angcal;

void define_MythenDetectorSpecifications_binding(py::module &m) {

    py::class_<MythenDetectorSpecifications,
               std::shared_ptr<MythenDetectorSpecifications>>(
        m, "MythenDetectorSpecifications", 
        R"(
        Attributes
        ----------

        pitch: double
            Read-only static. Strip/channel width of Mythen detector [mm] (0.05 mm)

        transverse_width: double
            Read-only static. Transverse width of Mythen detector [mm] (8.0 mm)

        strips_per_module: int
            Read-only static. Strips per module for Mythen detector (1280)

        max_modules: int
            Number of modules in detector. Default (48)

        sample_detector_offset: double
            Offset between sample horizontal plane and detector [degrees] (default: 1.4715°)

        offset: double
            Additional offset to sample detector offset (can change in experimental setup) [degrees] (default: 0.0°)

        dead_time: double
            Measured dead-time [s] (default: 76.08e-9 s)

        average_distance_sample_pixel: double
            average euclidean distance between sample and pixel [mm] (default: 2500.0 / pi mm)

        unconnected_modules: list of int
            list of unconnected modules

        Methods
        -------

        num_strips: int
            total number of strips in Mythen detector
        )")

        .def(py::init<>())

        .def_property_readonly_static(
            "pitch",
            [](py::object) { return MythenDetectorSpecifications::pitch; })

        .def_property_readonly_static(
            "transverse_width",
            [](py::object) {
                return MythenDetectorSpecifications::transverse_width;
            })

        .def_property_readonly_static(
            "strips_per_module",
            [](py::object) {
                return MythenDetectorSpecifications::strips_per_module;
            })

        .def_readwrite("max_modules",
                       &MythenDetectorSpecifications::max_modules)

        .def_readwrite("sample_detector_offset",
                       &MythenDetectorSpecifications::sample_detector_offset)

        .def_readwrite("offset", &MythenDetectorSpecifications::offset)

        .def_readwrite("dead_time", &MythenDetectorSpecifications::dead_time)

        .def_readwrite(
            "average_distance_sample_pixel",
            &MythenDetectorSpecifications::average_distance_sample_pixel)

        .def_readwrite("unconnected_modules",
                       &MythenDetectorSpecifications::unconnected_modules)

        .def("num_strips", &MythenDetectorSpecifications::num_strips);
}