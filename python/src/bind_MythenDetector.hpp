#include "MythenDetector.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "np_helpers.hpp"

namespace py = pybind11;

using namespace angcal;

using DGParametersArray =
    py::array_t<double, py::array::c_style | py::array::forcecast>;

inline void set_DGparameters_from_array(MythenDetector &self,
                                        DGParametersArray dgparameters) {
    auto view = make_view_2d(dgparameters);
    if (view.shape(0) != static_cast<ssize_t>(self.num_modules()) ||
        view.shape(1) != 3) {
        throw py::value_error(
            fmt::format("Expected array of shape ({}, 3), got ({}, {})",
                        self.num_modules(), view.shape(0), view.shape(1)));
    }
    DGParameters parameters;
    parameters.parameters = NDArray<double, 2>(view);
    self.set_DGparameters(parameters);
}

void define_MythenDetector_binding(py::module &m) {

    py::class_<MythenDetector, std::shared_ptr<MythenDetector>>(
        m, "MythenDetector",
        R"(
        Attributes
        ----------

        pitch: double
            Read-only static. Strip/channel width of Mythen detector [mm] (0.05 mm)

        strip_transverse_width: double
            Read-only static. Transverse width of a Mythen strip [mm] (8.0 mm)

        strips_per_module: int
            Read-only static. Strips per module for Mythen detector (1280)

        max_modules: int
            Read-only. Number of modules in detector, set in constructor. Default (48)

        sample_detector_offset: double
            Offset between sample horizontal plane and detector [degrees] (default: 1.4715°)

        offset: double
            Additional offset to sample detector offset (can change in experimental setup) [degrees] (default: 0.0°)

        sample_x_displacement: double
            Displacement of sample in x direction of horizontal sample plane [mm] (default: 0.0 mm)
        
        sample_y_displacement: double
            Displacement of sample in y direction of horizontal sample plane [mm] (default: 0.0 mm)

        elastic_correction_factor: double
            Elastic correction factor for torsional compliance correction model (default: 0.0)
        
        detector_vertical_axis_offset: double
            Offset of detector to vertical axis (used for elastic correction) [degrees] (default: 0.0°)

        dead_time: double
            Measured dead-time [s] (default: 76.08e-9 s)

        average_distance_sample_pixel: double
            average euclidean distance between sample and pixel [mm] (default: 2500.0 / pi mm)

        unconnected_modules: list of int
            list of unconnected modules
        )")

        .def(py::init<size_t>(), py::arg("num_modules") = 48)

        .def_property_readonly_static(
            "pitch", [](py::object) { return MythenDetector::pitch(); })

        .def_property_readonly_static(
            "strip_transverse_width",
            [](py::object) { return MythenDetector::strip_transverse_width(); })

        .def_property_readonly_static(
            "strips_per_module",
            [](py::object) { return MythenDetector::strips_per_module(); })

        .def_property_readonly("num_modules", &MythenDetector::num_modules)

        .def_property(
            "sample_detector_offset",
            py::overload_cast<>(&MythenDetector::sample_detector_offset,
                                py::const_),
            py::overload_cast<double>(&MythenDetector::sample_detector_offset),
            R"(Offset between sample horizontal plane and detector [degrees])")

        .def_property(
            "offset", py::overload_cast<>(&MythenDetector::offset, py::const_),
            py::overload_cast<double>(&MythenDetector::offset),
            R"(Additional offset to sample detector offset (can change in experimental setup) [degrees])")

        .def_property(
            "elastic_correction_factor",
            py::overload_cast<>(&MythenDetector::elastic_correction_factor,
                                py::const_),
            py::overload_cast<double>(
                &MythenDetector::elastic_correction_factor),
            R"(Elastic correction factor for torsional compliance correction model)")

        .def_property(
            "detector_vertical_axis_offset",
            py::overload_cast<>(&MythenDetector::detector_vertical_axis_offset,
                                py::const_),
            py::overload_cast<double>(
                &MythenDetector::detector_vertical_axis_offset),
            R"(Offset of detector to vertical axis (used for elastic correction) [degrees])")

        .def_property(
            "sample_x_displacement",
            py::overload_cast<>(&MythenDetector::sample_x_displacement,
                                py::const_),
            py::overload_cast<double>(&MythenDetector::sample_x_displacement),
            R"(Displacement of sample in x direction of horizontal sample plane [mm])")

        .def_property(
            "sample_y_displacement",
            py::overload_cast<>(&MythenDetector::sample_y_displacement,
                                py::const_),
            py::overload_cast<double>(&MythenDetector::sample_y_displacement),
            R"(Displacement of sample in y direction of horizontal sample plane [mm])")

        .def_property(
            "dead_time",
            py::overload_cast<>(&MythenDetector::dead_time, py::const_),
            py::overload_cast<double>(&MythenDetector::dead_time),
            R"(Measured dead-time [s])")

        .def_property(
            "average_distance_sample_pixel",
            py::overload_cast<>(&MythenDetector::average_distance_sample_pixel,
                                py::const_),
            py::overload_cast<double>(
                &MythenDetector::average_distance_sample_pixel),
            R"(Average euclidean distance between sample and pixel [mm] (default: 2500.0 / pi mm))")

        .def_property_readonly("num_strips", &MythenDetector::num_strips, R"(
            Total number of strips in Mythen detector
        )")

        .def_property(
            "DGparameters",
            py::overload_cast<>(&MythenDetector::get_DGparameters, py::const_),
            [](MythenDetector &self, py::object dgparameters) {
                if (py::isinstance<DGParameters>(dgparameters)) {
                    self.set_DGparameters(dgparameters.cast<DGParameters &>());
                } else {
                    set_DGparameters_from_array(
                        self, dgparameters.cast<DGParametersArray>());
                }
            },
            R"(
            Get or set the historical DG parameters.
            Accepts a DGParameters object or a numpy.ndarray of float with shape (num_modules, 3)
        )")

        .def("set_DGparameters", &MythenDetector::set_DGparameters,
             py::arg("dgparameters"), R"(
             Set the historical DG parameters

             Parameters
             ----------
                dgparameters : DGParameters
                    DGParameters object
             )")

        .def("set_DGparameters", &set_DGparameters_from_array,
             py::arg("dgparameters"), R"(
             Set the historical DG parameters from a numpy array

             Parameters
             ----------
                dgparameters : numpy.ndarray of float, shape (num_modules, 3)
                    Each row corresponds to a module, columns are (center, conversion, offset)
             )")

        .def("diffraction_angle_from_DG_parameters",
             &MythenDetector::diffraction_angle_from_DG_parameters,
             py::arg("module_index"), py::arg("detector_angle"),
             py::arg("strip_index"), py::arg("distance_to_strip") = 0.0, R"(
             Calculates the diffraction angle of a strip from DG parameters
             Parameters
             ----------
                module_index : int
                    Index of the module.
                detector_angle : float
                    Detector position, measured as the offset of the first strip from the default detector position [degrees]
                strip_index : int
                    Local strip index of module e.g. 0-1279
                distance_to_strip : float, optional
                    Distance to strip [given in strips]. Default is 0.0.
             Returns
             -------
                float
                    Diffraction angle of strip calculated from DG parameters [degrees]
             )")

        .def("diffraction_angle_from_BC_parameters",
             &MythenDetector::diffraction_angle_from_BC_parameters,
             py::arg("module_index"), py::arg("detector_angle"),
             py::arg("strip_index"), py::arg("distance_to_strip") = 0.0, R"(
             Calculates the diffraction angle of a strip from BC parameters
             Parameters
             ----------
                module_index : int
                    Index of the module.
                detector_angle : float 
                    Detector position, measured as the offset of the first strip from the default detector position [degrees]
                strip_index : int
                    Local strip index of module e.g. 0-1279
                distance_to_strip : float, optional
                    Distance to strip [given in strips]. Default is 0.0.
             Returns
             -------
                float
                    Diffraction angle of strip calculated from BC parameters [degrees]
             )")

        .def("diffraction_angle_from_EE_parameters",
             &MythenDetector::diffraction_angle_from_EE_parameters,
             py::arg("module_center_distance"), py::arg("normal_distance"),
             py::arg("angle"), py::arg("detector_angle"),
             py::arg("strip_index"), py::arg("distance_to_strip") = 0.0, R"(
             Calculates the diffraction angle of a strip from EE parameters
             Parameters
             ----------
                module_center_distance : float
                    Center distance of the module.
                normal_distance : float
                    Normal distance of the module.
                angle : float
                    Angle of the module.
                detector_angle : float
                    Detector position, measured as the offset of the first strip from the default detector position [degrees]
                strip_index : int
                    Local strip index of module e.g. 0-1279
                distance_to_strip : float, optional
                    Distance to strip [given in strips]. Default is 0.0.
             Returns
             -------
                float
                    Diffraction angle of strip calculated from EE parameters [degrees]
             )")

        .def("angular_strip_width_from_DG_parameters",
             &MythenDetector::angular_strip_width_from_DG_parameters,
             py::arg("module_index"), py::arg("local_strip_index"), R"(
             Calculated the strip width expressed as angle from DG module parameters

             Parameters
             -----------
                module_index : int
                    Index of the module.
                local_strip_index : int
                    Local strip index of module e.g. 0-1279

            Returns
            -------
                float
                    Strip width in angles [degrees]
            )")

        .def("angular_strip_width_from_BC_parameters",
             &MythenDetector::angular_strip_width_from_BC_parameters,
             py::arg("module_index"), py::arg("local_strip_index"), R"(
             Calculated the strip width expressed as angle from BC module parameters

             Parameters
             -----------
                module_index : int
                    Index of the module.
                local_strip_index : int
                    Local strip index of module e.g. 0-1279

            Returns
            -------
            float
                Strip width in angles [degrees]
            )")

        .def("angular_strip_width_from_EE_parameters",
             &MythenDetector::angular_strip_width_from_EE_parameters,
             py::arg("module_center_distance"), py::arg("normal_distance"),
             py::arg("angle"), py::arg("local_strip_index"), R"(
             Calculated the strip width expressed as angle from EE module parameters

             Parameters
             -----------
                module_center_distance : float
                    Center distance of the module.
                normal_distance : float
                    Normal distance of the module.
                angle : float
                    Angle of the module.
                local_strip_index : int
                    Local strip index of module e.g. 0-1279
                    
            Returns
            -------
            float
                Strip width in angles [degrees]
            )")

        .def("elastic_correction", &MythenDetector::elastic_correction,
             py::arg("detector_angle"), R"(
             Calculates the elastic correction factor for given detector angle
             based on torsional compliance correction model

             Parameters
             ----------
                detector_angle : float
                    Detector axis motor position [degrees]
            Returns
            -------
                float
                    Elastic correction factor for given detector angle [degrees]
            )")

        .def("transverse_width", &MythenDetector::transverse_width,
             py::arg("module_index"), py::arg("strip_index"),
             py::arg("from_BC_parameters"), R"(
             Calculates the transverse width of the strip projected onto a plane

             Parameters
             ----------
                module_index : int
                    Index of the module.
                strip_index : int
                    Local strip index of module e.g. 0-1279
                from_BC_parameters : bool
                    Whether to use BC parameters (true) or DG parameters (false)

             Returns
             -------
                float
                    Transverse width of strip projected onto a plane [mm]
            )")

        .def("solid_angle", &MythenDetector::solid_angle,
             py::arg("module_index"), py::arg("strip_index"), R"(
             Calculates the solid angle of a strip

             Parameters
             ----------
                module_index : int
                    Index of the module.
                strip_index : int
                    Local strip index of module e.g. 0-1279

             Returns
             -------
                float
                    Solid angle of strip [sr]
            )");
}