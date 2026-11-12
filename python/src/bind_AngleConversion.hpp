#include "AngularConversion.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "np_helpers.hpp"

namespace py = pybind11;

using namespace angcal;

void define_AngularConversion_binding(py::module &m) {

    py::class_<HistogramMetaData>(m, "HistogramMetaData")
        .def(py::init<>())
        .def(py::init([](std::tuple<double, double> angular_range,
                         double bin_width) {
                 return HistogramMetaData{angular_range, bin_width};
             }),
             py::arg("angular_range"), py::arg("bin_width"), R"(
             Parameters
             ----------
             angular_range: tuple of float
                (min_angle, max_angle) angular range for conversion [degrees]
             bin_width: float
                bin width of fixed angle width histogram [degrees]
             )")

        .def_readwrite(
            "angular_range", &HistogramMetaData::angular_range,
            R"(angular range for conversion [degrees] (default: (-180.0, 180.0)))")

        .def_readwrite(
            "bin_width", &HistogramMetaData::bin_width,
            R"(bin width of fixed angle width histogram [degrees] (default: 0.0036))")

        .def_property_readonly(
            "num_fixed_angle_width_bins",
            &HistogramMetaData::num_fixed_angle_width_bins,
            R"(number of bins in fixed angle width histogram)");

    m.def(
        "convert",
        [](const MythenDetector &mythen_detector,
           const HistogramMetaData &histogram_meta_data,
           py::array_t<uint32_t, py::array::c_style> photon_counts,
           const double detector_angle, const uint64_t incident_intensity,
           const double exposure_time,
           py::array_t<bool, py::array::c_style> bad_channels,
           py::array_t<double, py::array::c_style> flatfield,
           const double scale_factor) {
            std::cout << "in this overload" << std::endl;
            auto result = new NDArray<double, 2>(convert(
                mythen_detector, histogram_meta_data,
                make_view_1d(photon_counts), detector_angle, incident_intensity,
                exposure_time, make_view_1d(bad_channels),
                make_view_2d(flatfield), scale_factor));
            return return_image_data(result);
        },
        py::arg("mythen_detector"), py::arg("histogram_meta_data"),
        py::arg("photon_counts"), py::arg("detector_angle"),
        py::arg("incident_intensity"), py::arg("exposure_time"),
        py::arg("bad_channels"), py::arg("flatfield"),
        py::arg("scale_factor") = 1.0, R"(
        convert a frame to a diffraction pattern. The photon counts are corrected for flatfield, rate, incident intensity and solid angle, redistributed to fixed angle width bins and averaged

        Parameters
        ----------
        photon_counts: numpy.ndarray of uint32, shape (n_channels,)
            photon counts of the frame
        detector_angle: float
            detector angle of the frame [degrees]
        incident_intensity: int
            incident intensity of the frame
        see the multi-frame overload for the remaining parameters
        )");

    m.def(
        "convert",
        [](const MythenDetector &mythen_detector,
           const HistogramMetaData &histogram_meta_data,
           py::array_t<uint32_t, py::array::c_style> photon_counts,
           py::array_t<double, py::array::c_style> detector_angles,
           py::array_t<uint64_t, py::array::c_style> incident_intensities,
           const double exposure_time,
           py::array_t<bool, py::array::c_style> bad_channels,
           py::array_t<double, py::array::c_style> flatfield,
           const double scale_factor) {
            auto result = new NDArray<double, 2>(convert(
                mythen_detector, histogram_meta_data,
                make_view_2d(photon_counts), make_view_1d(detector_angles),
                make_view_1d(incident_intensities), exposure_time,
                make_view_1d(bad_channels), make_view_2d(flatfield),
                scale_factor));

            return return_image_data(result);
        },
        py::arg("mythen_detector"), py::arg("histogram_meta_data"),
        py::arg("photon_counts"), py::arg("detector_angles"),
        py::arg("incident_intensities"), py::arg("exposure_time"),
        py::arg("bad_channels"), py::arg("flatfield"),
        py::arg("scale_factor") = 1.0, R"(
        convert frames of photon counts to a diffraction pattern. The photon counts are corrected for flatfield, rate, incident intensity and solid angle, redistributed to fixed angle width bins and averaged

        Parameters
        ----------
        mythen_detector: MythenDetector
            storing all mythen specific parameters
        histogram_meta_data: HistogramMetaData
            angular range and bin width of fixed angle width histogram
        photon_counts: numpy.ndarray of uint32, shape (n_frames, n_channels)
            photon counts of each frame
        detector_angles: numpy.ndarray of float, shape (n_frames,)
            detector angle of each frame [degrees]
        incident_intensities: numpy.ndarray of uint64, shape (n_frames,)
            incident intensity of each frame
        exposure_time: float
            exposure time for each frame [s]
        bad_channels: numpy.ndarray of bool, shape (n_channels,)
            ``True`` if the channel is bad, otherwise ``False``
        flatfield: numpy.ndarray of float, shape (n_channels, 2)
            first column normalized flatfield value, second column its standard deviation
        scale_factor: float
            scale factor to scale correction to reasonable values (default: 1.0)

        Returns
        -------
        numpy.ndarray (num_fixed_angle_width_bins, 2)
            diffraction pattern in the range [min_angle, max_angle], first dimension photon counts redistributed to fixed angle width bins, second dimension respective standard deviation, values of -1.0 denote bins with insufficient coverage or which are denoted by a bad channel)");
}
