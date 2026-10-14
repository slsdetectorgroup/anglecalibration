
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
        .def(py::init<std::shared_ptr<MythenDetector>,
                      std::shared_ptr<FlatField>>(),
             py::arg("MythenDetector"), py::arg("FlatField"), R"(
             Parameters
             ----------
             MythenDetector: MythenDetector 
                storing all mythen specific parameters
             FlatField: FlatField 
                class storing inverse normalized flatfield
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

        .def_property(
            "angular_range",
            [](AngleCalibration &self) { return self.get_angular_range(); },
            [](AngleCalibration &self, std::pair<double, double> angle_range) {
                self.set_angular_range(angle_range.first, angle_range.second);
            },
            R"(
            angular range for conversion [degrees]

            Returns
            -------
            tuple of float
                (min_angle, max_angle)
            )")

        .def_property_readonly(
            "DGparameters",
            [](AngleCalibration &self) {
                auto DGparameters =
                    new NDArray<double, 2>(self.get_DGparameters().parameters);
                return return_image_data(DGparameters);
            },
            R"(
            DG parameters for each module

            Returns
            -------
            numpy.ndarray of float, shape (n_modules, 3)
                Each row corresponds to a module and contains the parameters:
                [center, conversion, offset]
            )")

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

        .def_property(
            "bad_channels",
            [](AngleCalibration &self) {
                auto bad_channels =
                    new NDArray<bool, 1>(self.get_bad_channels());
                return return_image_data(bad_channels);
            },
            [](AngleCalibration &self,
               py::array_t<bool, py::array::c_style> bad_channels) {
                auto temp_array = make_view_1d(bad_channels);
                self.set_bad_channels(
                    temp_array); // second copy TODO im copying twice
            },
            R"(
            bad_channels : numpy.ndarray of bool, shape (n_channels,)
                Expected size: number of channels/strips in the detector.
                Each element is ``True`` if the channel is bad, otherwise ``False``.
            )")

        .def(
            "calibrate",
            [](AngleCalibration &self,
               py::array_t<uint32_t, py::array::c_style> frames,
               py::array_t<double, py::array::c_style> detector_angles,
               const double base_peak_angle, const size_t module_index,
               const bool plot_calibration_process = false) {
                auto frames_view = make_view_2d(frames);
                auto detector_angles_view = make_view_1d(detector_angles);
                if (plot_calibration_process) {
                    self.calibrate<true>(frames_view, detector_angles_view,
                                         base_peak_angle, module_index);
                } else {
                    self.calibrate(frames_view, detector_angles_view,
                                   base_peak_angle, module_index);
                }
            },
            py::arg("frames"), py::arg("detector_angles"),
            py::arg("base_peak_angle"), py::arg("module_index"),
            py::arg("plot_calibration_process") = false,
            R"(
            calibrates BC parameters for respective module

            frames: numpy.ndarray of uint32, shape (n_frames, n_channels)
                Array of photon counts for each acquisition frame and channel. (acquisitions need to be ordered by motor position)
            detector_angles: numpy.ndarray of double, shape (n_frames,)
                Array of motor positions corresponding to each frame.
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
               py::array_t<uint32_t, py::array::c_style> frames,
               py::array_t<double, py::array::c_style> detector_angles,
               const double base_peak_angle,
               const bool plot_calibration_process = false) {
                auto frames_view = make_view_2d(frames);
                auto detector_angles_view = make_view_1d(detector_angles);
                if (plot_calibration_process) {
                    self.calibrate<true>(frames_view, detector_angles_view,
                                         base_peak_angle);
                } else {
                    self.calibrate(frames_view, detector_angles_view,
                                   base_peak_angle);
                }
            },
            py::arg("frames"), py::arg("detector_angles"),
            py::arg("base_peak_angle"),
            py::arg("plot_calibration_process") = false,
            R"(
            calibrates BC parameters for all modules
            
            frames : numpy.ndarray of uint32, shape (n_frames, n_channels)
                Array of photon counts for each acquisition frame and channel. (acquisitions need to be ordered by motor position)
            detector_angles : numpy.ndarray of double, shape (n_frames,)
                Array of motor positions corresponding to each frame.
            base_peak_angle : double
                angle of base peak center [degree]
            plot_calibration_process: bool, default: false
                Whether to plot the calibration process
            )");
}
