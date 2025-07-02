
#include "AngleCalibration.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "np_helpers.hpp"

namespace py = pybind11;

using namespace angcal;

void define_AngleCalibration_binding(py::module &m) {

    py::class_<AngleCalibration>(m, "AngleCalibration")
        .def(py::init<std::shared_ptr<MythenDetectorSpecifications>,
                      std::shared_ptr<FlatField>,
                      std::shared_ptr<MythenFileReader>,
                      std::optional<std::shared_ptr<SimpleFileInterface>>>(),
             py::arg("MythenDetectorSpecifications"), py::arg("FlatField"),
             py::arg("MythenFileReader"),
             py::arg("file_interface") = std::nullopt)

        .def_property("histogram_bin_width",
                      &AngleCalibration::get_histogram_bin_width,
                      &AngleCalibration::set_histogram_bin_width)

        .def_property_readonly("num_bins", &AngleCalibration::get_new_num_bins)
        .def("read_initial_calibration_from_file",
             [](AngleCalibration &self, const std::string &filename) {
                 self.read_initial_calibration_from_file(filename);
             })

        .def_property_readonly("DGparameters",
                               [](AngleCalibration &self) {
                                   auto DGparameters = new NDArray<double, 2>(
                                       self.get_DGparameters().parameters);
                                   return return_image_data(DGparameters);
                               }) // should I have a python class with method
                                  // centers, etc? - use py::buffer

        .def_property_readonly("new_photon_counts",
                               [](AngleCalibration &self) {
                                   auto new_photon_counts =
                                       new NDArray<double, 1>(
                                           self.get_new_photon_counts());
                                   return return_image_data(new_photon_counts);
                               })

        .def_property_readonly(
            "new_statistical_errors",
            [](AngleCalibration &self) {
                auto new_photon_count_errors =
                    new NDArray<double, 1>(self.get_new_statistical_errors());
                return return_image_data(new_photon_count_errors);
            })

        .def("calculate_fixed_bin_angle_width_histogram",
             [](AngleCalibration &self, const size_t start_frame_index,
                const size_t end_frame_index) {
                 self.calculate_fixed_bin_angle_width_histogram(
                     start_frame_index, end_frame_index);
             });
}
