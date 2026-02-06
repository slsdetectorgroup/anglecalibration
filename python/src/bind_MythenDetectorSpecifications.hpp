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

        .def(py::init<const double, const size_t, const size_t>(),
             py::arg("offset"), py::arg("num_counters") = 1,
             py::arg("max_modules") = 48,
             R"(
             Parameters
             ----------
             offset: double 
                Additional offset to sample detector offset. 
             num_counters: int
                Number of counters active. 
             max_modules: int 
                Number of modules in detector (default 48).
            )")

        .def_property(
            "unconnected_modules",
            [](MythenDetectorSpecifications &self) {
                return self.get_unconnected_modules();
            },
            [](MythenDetectorSpecifications &self,
               std::vector<ssize_t> &unconnected_modules) {
                self.set_unconnected_modules(unconnected_modules);
            })

        .def_property_readonly_static(
            "pitch",
            [](py::object) { return MythenDetectorSpecifications::pitch(); })

        .def_property_readonly_static(
            "strips_per_module",
            [](py::object) {
                return MythenDetectorSpecifications::strips_per_module();
            })

        .def_property_readonly("max_modules",
                               &MythenDetectorSpecifications::max_modules)

        .def_property_readonly("num_counters",
                               &MythenDetectorSpecifications::num_counters)

        .def_property_readonly_static(
            "sample_detector_offset",
            [](py::object) {
                return MythenDetectorSpecifications::sample_detector_offset();
            })

        .def_property_readonly("offset", &MythenDetectorSpecifications::offset)

        .def_property_readonly("num_strips",
                               &MythenDetectorSpecifications::num_strips);
}