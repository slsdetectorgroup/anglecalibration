#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Parameters.hpp"

namespace py = pybind11;

void define_Parameters_binding(py::module &m) {

    py::class_<Parameters>(m, "Parameters")
        .def(py::init<>())

        .def("num_modules", &Parameters::num_modules,
             R"(

            get number of modules
            
            )");

    py::class_<DGParameters, Parameters>(m, "DGParameters")
        .def(py::init<>())
        .def(py::init<const ssize_t>(), py::arg("num_modules"), R"(
             Constructor for DGParameters

             Parameters
            ----------
            num_modules: int
                number of modules to initialize parameters for
             )")

        .def(
            "convert_to_BCParameters",
            [](DGParameters &self, BCParameters &bcparameters) {
                self.convert_to_BCParameters(bcparameters);
            },
            py::arg("bcparameters"),
            R"(

            converts DG parameters to BC parameters and stores them in bcparameters

            Parameters
            ----------

            bcparameters: BCParameters
                BCParameters object to store converted BC parameters
            )")

        .def(
            "convert_to_EEParameters",
            [](DGParameters &self, EEParameters &eeparameters) {
                self.convert_to_EEParameters(eeparameters);
            },
            py::arg("eeparameters"),
            R"(

            converts DG parameters to EE parameters and stores them in eeparameters

            Parameters
            ----------

            eeparameters: EEParameters
                EEParameters object to store converted EE parameters
            )")

        .def(
            "parameters",
            [](DGParameters &self) {
                return return_view_data(
                    self.parameters.view(),
                    py::cast(self)); // maybe return memory view
            },
            R"(
                parameters as numpy array
    
                Returns
                -------

                numpy.ndarray (,3)
                    parameters stored as numpy array with shape (num_modules, 3) where the columns correspond to the respective parameters (e.g. center, conversion, offset)
                )");

    py::class_<EEParameters, Parameters>(m, "EEParameters")
        .def(py::init<>())

        .def(py::init<const ssize_t>(), py::arg("num_modules"), R"(
             Constructor for EEParameters

             Parameters
            ----------
            num_modules: int
                number of modules to initialize parameters for
             )")

        .def(
            "parameters",
            [](EEParameters &self) {
                return return_view_data(
                    self.parameters.view(),
                    py::cast(self)); // maybe return memory view
            },
            R"(
                parameters as numpy array
    
                Returns
                -------

                numpy.ndarray (,3)
                    parameters stored as numpy array with shape (num_modules, 3) where the columns correspond to the respective parameters (e.g. normal_distances, module_center_distances, angles)
                )");

    py::class_<BCParameters, Parameters>(m, "BCParameters")
        .def(py::init<>())

        .def(py::init<const ssize_t>(), py::arg("num_modules"), R"(
             Constructor for BCParameters

             Parameters
            ----------
            num_modules: int
                number of modules to initialize parameters for
             )")

        .def(
            "convert_to_DGParameters",
            [](BCParameters &self, DGParameters &dgparameters) {
                self.convert_to_DGParameters(dgparameters);
            },
            py::arg("dgparameters"),
            R"(
            converts BC parameters to DG parameters and stores them in dgparameters

            Parameters
            ----------

            dgparameters: DGParameters
                DGParameters object to store converted DG parameters
            )")

        .def(
            "convert_to_EEParameters",
            [](BCParameters &self, EEParameters &eeparameters) {
                self.convert_to_EEParameters(eeparameters);
            },
            py::arg("eeparameters"),
            R"(
            converts BC parameters to EE parameters and stores them in eeparameters

            Parameters
            ----------

            eeparameters: EEParameters
                EEParameters object to store converted EE parameters
            )")

        .def(
            "parameters",
            [](BCParameters &self) {
                return return_view_data(
                    self.parameters.view(),
                    py::cast(self)); // maybe return memory view
            },
            R"(
                parameters as numpy array
    
                Returns
                -------

                numpy.ndarray (,3)
                    parameters stored as numpy array with shape (num_modules, 3) where the columns correspond to the respective parameters (e.g. angle_center_module_normal, module_center_sample_distances, angle_center_beam)
                )");
}
