
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
                      &AngleCalibration::set_histogram_bin_width, R"(
                      bin width of fixed angle width histogram [degrees]
                      default: '0.0036°')")

        .def_property("base_peak_ROI_width",
                      &AngleCalibration::get_base_peak_ROI_width,
                      &AngleCalibration::set_base_peak_ROI_width, R"(
                      width of base peak region of interest [degrees]
                      e.g. [base_peak - base_peak_ROI_width, base_peak + base_peak_ROI_width] given in angles
                      default: '0.05°')")

        .def_property_readonly("base_peak_ROI_num_bins",
                               &AngleCalibration::get_base_peak_ROI_num_bins,
                               R"(
                             number of bins covered by base peak region of interest
                             default: '101' bins)")

        .def_property_readonly(
            "num_fixed_angle_width_bins",
            &AngleCalibration::num_fixed_angle_width_bins,
            R"(number of bins in fixed angle width histogram)")

        .def_property("scale_factor", &AngleCalibration::get_scale_factor,
                      &AngleCalibration::set_scale_factor,
                      R"(scale factor to scale correction to reasonable values
                      )")

        .def_property(
            "angular_range",
            [](AngleCalibration &self) { return self.get_angular_range(); },
            [](AngleCalibration &self, std::pair<double, double> angle_range) {
                self.set_angular_range(angle_range.first, angle_range.second);
            },
            R"(angular range for conversion [degrees]
            Returns
            -------
            tuple of float
                (min_angle, max_angle)
            )")

        .def(
            "read_initial_calibration_from_file",
            [](AngleCalibration &self, const std::string &filename) {
                self.read_initial_calibration_from_file(filename);
            },
            R"(reads the historical Detector Group (DG) parameters from file and
     transforms them to Best Computing parameters)")

        .def(
            "read_bad_channels_from_file",
            [](AngleCalibration &self, const std::string &filename) {
                self.read_bad_channels_from_file(filename);
            },
            py::arg("filename"), R"(reads bad channels from file)")

        .def_property("base_peak_angle", &AngleCalibration::get_base_peak_angle,
                      &AngleCalibration::set_base_peak_angle,
                      R"(center of chosen base peak for calibration [degrees])")

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

        .def(
            "module_is_disconnected",
            [](AngleCalibration &self, const size_t module_index) {
                return self.module_is_disconnected(module_index);
            },
            py::arg("module_index"),
            R"(
            check if a module only has bad channels or is disconnected 

            Parameters
            ----------
            module_index : int
                Index of the module.
            Returns
            -------
            bool
                True if module is disconnected, False otherwise.
            )")

        .def(
            "rate_correction",
            [](const AngleCalibration &self, const double photon_count,
               const double photon_count_error, const double exposure_time) {
                return self.rate_correction(photon_count, photon_count_error,
                                            exposure_time);
            },
            py::arg("photon_count"), py::arg("photon_count_error"),
            py::arg("exposure_time"),
            R"(
            performs rate correction

            Parameters
            ----------
            photon_count : float
                measured photon counts
            photon_count_error : float
                propagated error of measured photon counts
            exposure_time : float
                exposure time of acquisition [s]

            Returns
            -------
            tuple of float
                pair {rate corrected photon counts, propagated_error}
            
            )")

        .def_property_readonly(
            "DGparameters",
            [](AngleCalibration &self) {
                return self.get_DGparameters(); // TODO: in c++ actually returns
                                                // a const reference, what
                                                // return value policy to use?
            },
            R"(
            historic DG parameters 
            )")

        .def_property_readonly(
            "MythenDetectorSpecifications",
            [](AngleCalibration &self) {
                return self.get_detector_specifications();
            },
            R"(
            MythenDetectorSpecifications storing all mythen specific parameters
            )")

        .def_property_readonly(
            "BCparameters",
            [](AngleCalibration &self) { return self.get_BCparameters(); },
            R"(
            current "Best Computing" BC parameters 
            )")

        .def_property(
            "bad_channels",
            [](AngleCalibration &self) {
                auto bad_channels =
                    new NDArray<ssize_t, 1>(self.get_bad_channels());
                return return_image_data(bad_channels);
            },
            [](AngleCalibration &self,
               py::array_t<bool, py::array::forcecast> bad_channels) {
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
            },
            R"(
            bad_channels : numpy.ndarray of bool, shape (n_channels,)
                Expected size: number of channels/strips in the detector.
                Each element is ``True`` if the channel is bad, otherwise ``False``.
            )")

        .def(
            "diffraction_angle_from_DG_parameters",
            [](AngleCalibration &self, const size_t module_index,
               const double detector_angle, const size_t strip_index,
               const double distance_to_strip) {
                return self.diffraction_angle_from_DG_parameters(
                    module_index, detector_angle, strip_index,
                    distance_to_strip);
            },
            py::arg("module_index"), py::arg("detector_angle"),
            py::arg("strip_index"), py::arg("distance_to_strip"),
            R"(
            calculates diffraction angle from DG parameters for given strip

            Parameters
            ----------
            module_index : int
                index of module
            detector_angle : double
                detector position, measured as the offset of the first strip from the default detector position [degrees]
            strip_index : int
                index of strip in module
            distance_to_strip : double
                distance to strip [given in strips]

            Returns
            -------
            double
                diffraction angle for given strip [degrees]
            )")

        .def(
            "calibrate",
            [](AngleCalibration &self,
               const std::vector<std::string> &file_list,
               const double base_peak_angle, const size_t module_index,
               const bool plot_calibration_process = false) {
                if (plot_calibration_process) {
                    self.calibrate<true>(file_list, base_peak_angle,
                                         module_index);
                } else {
                    self.calibrate(file_list, base_peak_angle, module_index);
                }
            },
            py::arg("file_list"), py::arg("base_peak_angle"),
            py::arg("module_index"),
            py::arg("plot_calibration_process") = false,
            R"(
            calibrates BC parameters for respective module

            file_list: list 
                List of paths to acquisition files.
            base_peak_angle: float
                Angle of base peak center [degree].
            module_index: int
                Index of module
            plot_calibration_process: bool, default: false
                Whether to plot the calibration process 
            )")

        .def(
            "calibrate",
            [](AngleCalibration &self,
               const std::vector<std::string> &file_list,
               const double base_peak_angle,
               std::optional<std::string> output_filename = std::nullopt,
               const bool plot_calibration_process = false) {
                if (plot_calibration_process) {
                    self.calibrate<true>(file_list, base_peak_angle,
                                         output_filename);
                } else {
                    self.calibrate(file_list, base_peak_angle, output_filename);
                }
            },
            py::arg("file_list"), py::arg("base_peak_angle"),
            py::arg("output_filename"),
            py::arg("plot_calibration_process") = false,
            R"(
            calibrates BC parameters for all modules
            
            file_list : list
                list of paths to acquisition files
            base_peak_angle : double
                angle of base peak center [degree]
            output_filename : str, optional 
                if given, writes calibrated DG parameters to file with given name
            plot_calibration_process: bool, default: false
                Whether to plot the calibration process
            )")

        .def(
            "convert",
            [](AngleCalibration &self,
               const std::vector<std::string> &file_list) {
                auto result = new NDArray<double, 1>(self.convert(file_list));
                return return_image_data(result);
            },
            py::arg("file_list"),
            R"(
            performs angular conversion e.g. calculates diffraction pattern from raw photon counts. If module_index is given calculates diffraction pattern for a specific module

            Parameters
            ----------

            file_list: list 
                list of paths to acquisition files
            module_index: int, optional
                index of module

            Returns
            -------
            numpy.ndarray (,num_fixed_angle_width_bins)
                photon counts redistributed to fixed angle width bins, flatfield corrected and variance scaled photon counts)")

        .def(
            "convert",
            [](AngleCalibration &self,
               const std::vector<std::string> &file_list,
               const size_t module_index) {
                auto result = new NDArray<double, 1>(
                    self.convert(file_list, module_index));
                return return_image_data(result);
            },
            py::arg("file_list"), py::arg("module_index"))

        .def(
            "write_DG_parameters_to_file",
            [](AngleCalibration &self, const std::string &filename,
               const DGParameters &parameters) {
                self.write_DG_parameters_to_file(filename, parameters);
            },
            py::arg("filename"), py::arg("parameters"),
            R"(
            writes DG parameters to file

             Parameters
             ----------
             filename: str
                path to output file
             parameters: DGParameters
                DGParameters object containing the parameters to write to file
             )");
}
