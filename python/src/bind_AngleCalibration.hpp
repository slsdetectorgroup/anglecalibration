
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

    py::class_<AngleCalibration>(m, "AngleCalibration", R"(Attributes
----------
histogram_bin_width : double 
    bin width of fixed angle width histogram [degree]
    (default: 0.0036 deg)
base_peak_angle: double
    angle of center of base peak [degree]
number_of_bins : int 
    Read-only. number of bins for new fixed angle width histogram
)")
        .def(py::init<std::shared_ptr<MythenDetectorSpecifications>,
                      std::shared_ptr<FlatField>,
                      std::shared_ptr<MythenFileReader>>(),
             py::arg("MythenDetectorSpecifications"), py::arg("FlatField"),
             py::arg("MythenFileReader"), R"(
             Parameters
             ----------
             MythenDetectorSpecifications: MythenDetectorSpecifications 
                storing all mythen specific parameters
             FlatField: FlatField 
                class storing inverse normalized flatfield
             MythenFileReader: File Reader to read Mythen acquisition files
             )")

        .def_property("histogram_bin_width",
                      &AngleCalibration::get_histogram_bin_width,
                      &AngleCalibration::set_histogram_bin_width)

        .def_property_readonly("num_fixed_angle_width_bins",
                               &AngleCalibration::num_fixed_angle_width_bins)

        .def(
            "read_initial_calibration_from_file",
            [](AngleCalibration &self, const std::string &filename) {
                self.read_initial_calibration_from_file(filename);
            },
            R"(reads the historical Detector Group (DG) parameters from file and
     transforms them to Best Computing parameters)")

        .def_property("base_peak_angle", &AngleCalibration::get_base_peak_angle,
                      &AngleCalibration::set_base_peak_angle)

        .def(
            "base_peak_is_in_module",
            [](AngleCalibration &self, const size_t module_index,
               const double detector_angle) {
                return self.base_peak_is_in_module(module_index,
                                                   detector_angle);
            },
            py::arg("module_index"), py::arg("detector_angle"),
            R"(
            check if base peak ROI is contained within module region

            Parameters
            ----------
            module_index : int
                Index of the module.
            detector_angle : double
                Detector position, measured as the offset of the first strip from
                the default detector position [degrees].

            Returns
            -------
            bool
                True if the base peak ROI lies inside the module region, False otherwise.)")

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

        .def(
            "calibrate",
            [](AngleCalibration &self,
               const std::vector<std::string> &file_list,
               const double base_peak_angle, const size_t module_index) {
                self.calibrate(file_list, base_peak_angle, module_index);
            },
            py::arg("file_list"), py::arg("base_peak_angle"),
            py::arg("module_index"),
            R"(
            calibrates BC parameters for respective module 

            Parameters
            ----------
            file_list: list 
                list of paths to acquisition files
            base_peak_angle: double 
                angle of base peak center [degree]
            module_index: int
                index of module)")

        .def(
            "calibrate",
            [](AngleCalibration &self,
               const std::vector<std::string> &file_list,
               const double base_peak_angle) {
                self.calibrate(file_list, base_peak_angle);
            },
            py::arg("file_list"), py::arg("base_peak_angle"),
            R"(
            calibrates BC parameters for all modules

            Parameters
            ----------
            file_list: list 
                list of paths to acquisition files
            base_peak_angle: double
                angle of base peak center [degree]
            )")

        // TODO: is ite better to pass a string for the filename and expose
        // mythen_file_reader of AngleCalibration
        .def(
            "redistributed_photon_counts_in_base_peak_ROI",
            [](AngleCalibration &self, const MythenFrame &frame,
               const size_t module_index) {
                auto result = new NDArray<double, 1>(
                    self.redistributed_photon_counts_in_base_peak_ROI(
                        frame, module_index));
                return return_image_data(result);
            },
            R"(
            redistribute photon counts to fixed angle width bins which are
            within base peak region

            Returns
            -------
            numpy.ndarray (,number_of_bins_in_base_peak_ROI)
                to fixed angle width redistributed, flatfield corrected and variance scaled photon counts of respective module within base peak ROI)")

        .def(
            "redistribute_photon_counts_to_fixed_angle_width_bins",
            [](AngleCalibration &self, const MythenFrame &frame,
               const size_t module_index) {
                auto result = new NDArray<double, 1>(
                    self.redistribute_photon_counts_to_fixed_angle_width_bins(
                        frame, module_index));
                return return_image_data(result);
            },
            R"(
            redistribute photon counts of respective module fixed angle width bins 

            Returns
            -------
            numpy.ndarray (,num_fixed_angle_width_bins)
                to fixed angle width redistributed, flatfield corrected and variance scaled photon counts of respective module)")

        .def(
            "convert",
            [](AngleCalibration &self,
               const std::vector<std::string> &file_list) {
                auto result = new NDArray<double, 1>(self.convert(file_list));
                return return_image_data(result);
            },
            R"(
            performs angular conversion e.g. calculates diffraction pattern from raw photon counts

            Params: 
            ----------
            file_list: list 
                list of paths to acquisition files

            Returns
            -------
            numpy.ndarray (,num_fixed_angle_width_bins)
                to fixed angle width redistributed, flatfield corrected and variance scaled photon counts)");
}
