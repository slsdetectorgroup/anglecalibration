
#include "AngleCalibration.hpp"
#include <filesystem>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include "np_helpers.hpp"

namespace py = pybind11;

using namespace angcal;

void define_AngleCalibration_binding(py::module &m) {

    py::class_<AngleCalibration>(m, "AngleCalibration")
        .def(py::init<std::shared_ptr<MythenDetectorSpecifications>,
                      std::shared_ptr<FlatField>,
                      std::optional<std::shared_ptr<MythenFileReader>>,
                      std::optional<std::shared_ptr<SimpleFileInterface>>>(),
             py::arg("MythenDetectorSpecifications"), py::arg("FlatField"),
             py::arg("MythenFileReader") = std::nullopt,
             py::arg("file_interface") = std::nullopt)

        .def_property("histogram_bin_width",
                      &AngleCalibration::get_histogram_bin_width,
                      &AngleCalibration::set_histogram_bin_width)

        .def_property_readonly("new_number_of_bins",
                               &AngleCalibration::new_number_of_bins)
        .def("read_initial_calibration_from_file",
             [](AngleCalibration &self, const std::string &filename) {
                 self.read_initial_calibration_from_file(filename);
             })

        .def_property("base_peak_angle", &AngleCalibration::get_base_peak_angle,
                      &AngleCalibration::set_base_peak_angle)

        .def("base_peak_is_in_module",
             [](AngleCalibration &self, const size_t module_index,
                const double detector_angle) {
                 return self.base_peak_is_in_module(module_index,
                                                    detector_angle);
             })
        .def("module_is_disconnected",
             [](AngleCalibration &self, const size_t module_index) {
                 return self.module_is_disconnected(module_index);
             })

        // TODO set initial calibration from file

        .def_property_readonly(
            "DGparameters",
            [](AngleCalibration &self) {
                auto DGparameters =
                    new NDArray<double, 2>(self.get_DGparameters().parameters);
                return return_image_data(
                    DGparameters); // maybe return memoryview::from_memory
            })                     // should I have a python class with method
                                   // centers, etc? - use py::buffer

        .def("calibrate",
             [](AngleCalibration &self,
                const std::vector<std::string> &file_list,
                const double base_peak_angle, const size_t module_index) {
                 self.calibrate(file_list, base_peak_angle, module_index);
             })

        .def("calibrate",
             [](AngleCalibration &self,
                const std::vector<std::string> &file_list,
                const double base_peak_angle) {
                 self.calibrate(file_list, base_peak_angle);
             })

        // TODO: is ite better to pass a string for the filename and expose
        // mythen_file_reader of AngleCalibration
        .def("redistributed_photon_counts_in_base_peak_ROI",
             [](AngleCalibration &self, const MythenFrame &frame,
                const size_t module_index) {
                 auto result = new NDArray<double, 1>(
                     self.redistributed_photon_counts_in_base_peak_ROI(
                         frame, module_index));
                 return_image_data(result);
             })

        .def("redistribute_photon_counts_to_fixed_angle_width_bins",
             [](AngleCalibration &self, const MythenFrame &frame,
                const size_t module_index) {
                 auto result = new NDArray<double, 1>(
                     self.redistribute_photon_counts_to_fixed_angle_width_bins(
                         frame, module_index));
                 return_image_data(result);
             })

        .def("redistribute_photon_counts_to_fixed_angle_width_bins",
             [](AngleCalibration &self, const MythenFrame &frame) {
                 auto result = new NDArray<double, 1>(
                     self.redistribute_photon_counts_to_fixed_angle_width_bins(
                         frame));
                 return_image_data(result);
             });
}
