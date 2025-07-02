#include "bind_AngleCalibration.hpp"
#include "bind_FlatField.hpp"
#include "bind_MythenDetectorSpecifications.hpp"
#include "bind_SimpleFileInterface.hpp"

// Pybind stuff
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(_angcal, m) {
    define_AngleCalibration_binding(m);
    define_MythenDetectorSpecifications_binding(m);
    define_FlatField_binding(m);
    define_SimpleFileInterface_binding(m);
}