#include "MythenDetectorSpecifications.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;

using namespace angcal;

void define_MythenDetectorSpecifications_binding(py::module &m) {

    py::class_<MythenDetectorSpecifications,
               std::shared_ptr<MythenDetectorSpecifications>>(
        m, "MythenDetectorSpecifications",
        R"(Attributes
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
sample_detector_offset: double 
offset: double 
)")

        .def(py::init<>())

        .def_property_readonly_static(
            "pitch",
            [](py::object) { return MythenDetectorSpecifications::pitch; },
            R"(Strip/channel width of Mythen detector [mm] (0.05 mm))")

        .def_property_readonly_static(
            "transverse_width",
            [](py::object) {
                return MythenDetectorSpecifications::transverse_width;
            },
            R"(Transverse width of Mythen detector [mm] (8.0 mm))")

        .def_property_readonly_static(
            "strips_per_module",
            [](py::object) {
                return MythenDetectorSpecifications::strips_per_module;
            },
            R"(Strips/channels per module for Mythen detector (1280))")

        .def_readwrite("max_modules",
                       &MythenDetectorSpecifications::max_modules,
                       R"(Number of modules in detector. Default (48))")

        .def_readwrite(
            "sample_detector_offset",
            &MythenDetectorSpecifications::sample_detector_offset,
            R"(Offset between sample horizontal plane and detector [degrees] (default: 1.4715°))")

        .def_readwrite(
            "offset", &MythenDetectorSpecifications::offset,
            R"(Additional offset to sample detector offset (can change in experimental setup) [degrees] (default: 0.0°))")

        .def_readwrite("dead_time", &MythenDetectorSpecifications::dead_time,
                       R"(Measured dead-time [s] (default: 76.08e-9 s))")

        .def_readwrite(
            "average_distance_sample_pixel",
            &MythenDetectorSpecifications::average_distance_sample_pixel,
            R"(average euclidean distance between sample and pixel [mm] (default: 2500.0 / pi mm))")

        .def_readwrite("unconnected_modules",
                       &MythenDetectorSpecifications::unconnected_modules,
                       R"(list of unconnected modules)")

        .def("num_strips", &MythenDetectorSpecifications::num_strips,
             R"(Total number of strips in Mythen detector)");
}