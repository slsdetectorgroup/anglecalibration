#include "AngularConversion.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "np_helpers.hpp"

namespace py = pybind11;

using namespace angcal;

void define_AngularConversion_binding(py::module &m) {

    py::class_<AngularConversion>(m, "AngularConversion")
        .def(py::init<std::shared_ptr<MythenDetector>>(),
             py::arg("MythenDetector"), R"(
             Parameters
             ----------
             MythenDetector: MythenDetector
                storing all mythen specific parameters
             )")

        .def(py::init<std::shared_ptr<MythenDetector>, double,
                      std::tuple<double, double>>(),
             py::arg("MythenDetector"), py::arg("histogram_bin_width"),
             py::arg("angular_range"), R"(
             Parameters
             ----------
             MythenDetector: MythenDetector
                storing all mythen specific parameters
             histogram_bin_width: float
                bin width of fixed angle width histogram [degrees]
             angular_range: tuple of float
                (min_angle, max_angle) angular range for conversion [degrees]
             )")

        .def_property_readonly("histogram_bin_width",
                               &AngularConversion::get_histogram_bin_width,
                               R"(
                      bin width of fixed angle width histogram [degrees]
                      default: '0.0036°')")

        .def_property_readonly(
            "angular_range", &AngularConversion::get_angular_range,
            R"(angular range for conversion [degrees] (default: (-180.0, 180.0))
            Returns
            -------
            tuple of float
                (min_angle, max_angle)
            )")

        .def_property_readonly(
            "num_fixed_angle_width_bins",
            &AngularConversion::num_fixed_angle_width_bins,
            R"(number of bins in fixed angle width histogram)")

        .def_property("scale_factor", &AngularConversion::get_scale_factor,
                      &AngularConversion::set_scale_factor,
                      R"(scale factor to scale correction to reasonable values
                      )")

        .def_property(
            "flatfield",
            [](AngularConversion &self) {
                auto result = new NDArray<double, 2>(
                    self.get_FlatField()); // makes a copy
                return return_image_data(result);
            },
            [](AngularConversion &self,
               py::array_t<double, py::array::c_style> normalized_flatfield) {
                auto shape = get_checked_shape<2>(normalized_flatfield);
                NDView<const double, 2> flatfield_view(
                    normalized_flatfield.data(), shape);
                self.set_FlatField(flatfield_view);
            },
            R"(normalized_flatfield: numpy.ndarray of float, shape (n_channels, 2)
                first column normalized flatfield value, second column its standard deviation
            )")

        .def_property(
            "bad_channels",
            [](AngularConversion &self) {
                auto result = new NDArray<bool, 1>(
                    self.get_bad_channels()); // makes a copy
                return return_image_data(result);
            },
            [](AngularConversion &self,
               py::array_t<bool, py::array::c_style> bad_channels) {
                auto bad_channels_view = make_view_1d(bad_channels);
                self.set_bad_channels(bad_channels_view);
            },
            R"(bad_channels: numpy.ndarray of bool, shape (n_channels,)
                Each element is ``True`` if the channel is bad, otherwise ``False``.
            )")

        .def(
            "convert",
            [](AngularConversion &self,
               py::array_t<uint32_t, py::array::c_style> photon_counts,
               py::array_t<double, py::array::c_style> detector_angles,
               py::array_t<uint64_t, py::array::c_style> incident_intensities,
               const double exposure_time) {
                auto result = new NDArray<double, 2>(self.convert(
                    make_view_2d(photon_counts), make_view_1d(detector_angles),
                    make_view_1d(incident_intensities), exposure_time));
                return return_image_data(result);
            },
            py::arg("photon_counts"), py::arg("detector_angles"),
            py::arg("incident_intensities"), py::arg("exposure_time"), R"(
            convert frames of photon counts to a diffraction pattern. The photon counts are corrected for flatfield, rate, incident intensity and solid angle, redistributed to fixed angle width bins and averaged

            Parameters
            ----------
            photon_counts: numpy.ndarray of uint32, shape (n_frames, n_channels)
                photon counts of each frame
            detector_angles: numpy.ndarray of float, shape (n_frames,)
                detector angle of each frame [degrees]
            incident_intensities: numpy.ndarray of uint64, shape (n_frames,)
                incident intensity of each frame
            exposure_time: float
                exposure time for each frame [s]

            Returns
            -------
            numpy.ndarray (num_fixed_angle_width_bins, 2)
                diffraction pattern in the range [min_angle, max_angle], first dimension photon counts redistributed to fixed angle width bins, second dimension respective standard deviation, values of -1.0 denote bins with insufficient coverage or which are denoted by a bad channel)")

        .def(
            "convert",
            [](AngularConversion &self,
               py::array_t<uint32_t, py::array::c_style> photon_counts,
               const double detector_angle, const uint64_t incident_intensity,
               const double exposure_time) {
                auto result = new NDArray<double, 2>(
                    self.convert(make_view_1d(photon_counts), detector_angle,
                                 incident_intensity, exposure_time));
                return return_image_data(result);
            },
            py::arg("photon_counts"), py::arg("detector_angle"),
            py::arg("incident_intensity"), py::arg("exposure_time"), R"(
            convert a frame to a diffraction pattern. The photon counts are corrected for flatfield, rate, incident intensity and solid angle, redistributed to fixed angle width bins and averaged

            Parameters
            ----------
            photon_counts: numpy.ndarray of uint32, shape (n_channels,)
                photon counts of the frame
            detector_angle: float
                detector angle of the frame [degrees]
            incident_intensity: uint64
                incident intensity of the frame
            exposure_time: float
                exposure time for the frame [s]

            Returns
            -------
            numpy.ndarray (num_fixed_angle_width_bins, 2)
                diffraction pattern in the range [min_angle, max_angle], first dimension photon counts redistributed to fixed angle width bins, second dimension respective standard deviation, values of -1.0 denote bins with insufficient coverage or which are denoted by a bad channel)");
}
