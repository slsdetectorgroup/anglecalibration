#include "CalibrationRecorder.hpp"
#include "CalibrationStep.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;

using namespace angcal;

void define_CalibrationStep_binding(py::module &m) {

    py::enum_<CalibrationPhase>(m, "CalibrationPhase")
        .value("CoupledParameters", CalibrationPhase::CoupledParameters)
        .value("Offset", CalibrationPhase::Offset);

    py::class_<CalibrationStep>(m, "CalibrationStep", R"(
        state of the calibration of one module after an optimization step)")

        .def_readonly("module_index", &CalibrationStep::module_index)

        .def_readonly("phase", &CalibrationStep::phase)

        .def_readonly("iteration", &CalibrationStep::iteration, R"(
                      index of optimization step, '0' denotes the state before the first step)")

        .def_readonly("similarity", &CalibrationStep::similarity, R"(
                      similarity criterion of the base peaks)")

        .def_readonly("parameters", &CalibrationStep::parameters, R"(
                      BC parameters of the module: module center sample distance, angle between module center and module normal, angle between module center and beam)")

        .def_property_readonly(
            "bin_angles",
            [](const CalibrationStep &self) {
                return py::array_t<double>(self.bin_angles.size(),
                                           self.bin_angles.data());
            },
            R"(diffraction angle of each bin in base peak ROI [degrees])")

        .def_property_readonly(
            "detector_angles",
            [](const CalibrationStep &self) {
                return py::array_t<double>(self.detector_angles.size(),
                                           self.detector_angles.data());
            },
            R"(detector angle of the acquisition of each curve [degrees])")

        .def_property_readonly(
            "curve_modules",
            [](const CalibrationStep &self) {
                return py::array_t<size_t>(self.curve_modules.size(),
                                           self.curve_modules.data());
            },
            R"(index of the module each curve was measured with)")

        .def_property_readonly(
            "base_peak_curves",
            [](const CalibrationStep &self) {
                return py::array_t<double>(self.base_peak_curves.shape(),
                                           self.base_peak_curves.data());
            },
            R"(
            numpy.ndarray (curves, bins)
                photon counts redistributed to fixed angle width bins in base peak ROI, NaN where the base peak ROI is not covered)");

    py::class_<CalibrationRecorder>(m, "CalibrationRecorder", R"(
        records calibration steps to a hdf5 file. Can be passed to AngleCalibration.set_calibration_observer)")

        .def(py::init<const std::filesystem::path &>(), py::arg("filename"),
             R"(
             Parameters
             ----------
             filename: str
                path to hdf5 file, an existing file is overwritten
             )")

        .def("record", &CalibrationRecorder::record, py::arg("step"))

        .def("__call__", &CalibrationRecorder::record, py::arg("step"))

        .def("close", &CalibrationRecorder::close);
}
