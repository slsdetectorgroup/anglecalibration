#include "FlatField.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include "np_helpers.hpp"

namespace py = pybind11;

using namespace angcal;

void define_FlatField_binding(py::module &m) {

    py::class_<FlatField, std::shared_ptr<FlatField>>(m, "FlatField")
        .def(py::init<std::shared_ptr<MythenDetectorSpecifications>>(),
             py::arg("MythenDetectorSpecifications"), R"(

             Parameters
             ----------

             MythenDetectorSpecifications: MythenDetectorSpecifications 
                storing all mythen specific parameters
             )")

        .def(
            "create_normalized_flatfield_from_filelist",
            [](FlatField &self,
               const std::vector<std::filesystem::path> &file_list,
               std::shared_ptr<MythenFileReader> mythen_file_reader) {
                self.create_normalized_flatfield_from_filelist(
                    file_list, mythen_file_reader);
            },
            py::arg("file_list"), py::arg("mythen_file_reader"), R"(

                Create normalized flatfield from list of files containing flatfield acquisitions

                Parameters
                ----------
                
                file_list: list of str
                    list of paths to acquisition files containing flatfield acquisitions
                mythen_file_reader: MythenFileReader
                    file reader to read mythen acquisition files
                )")

        /*
        .def("read_flatfield_from_file",
             [](FlatField &self, const std::string &filename) {
                 self.read_flatfield_from_file(filename);
             })
        */

        .def_property(
            "normalized_flatfield",
            [](FlatField &self) {
                auto result =
                    new NDArray<double, 2>(self.get_normalized_flatfield());
                return return_image_data(result); // maybe return memory view
            },
            [](FlatField &self,
               py::array_t<double, py::array::c_style | py::array::forcecast>
                   &normalized_flatfield) {
                py::buffer_info info = normalized_flatfield.request();
                NDView<double, 2> temp_array_view(
                    reinterpret_cast<double *>(info.ptr),
                    std::array<ssize_t, 2>{info.shape[0], info.shape[1]});
                NDArray<double, 2> temp_array(
                    temp_array_view); // first copy //maybe modify buffer //or
                                      // iterate
                self.set_normalized_flatfield(
                    temp_array); // second copy TODO im copying twice
            },
            R"(
            numpy.ndarray of float, shape (n_channels, 2) : Each row corresponds to a strip/channel in the detector. The first column contains the normalized flatfield value for that strip, and the second column contains the standard deviation of the flatfield value for that strip.
             )")

        .def_property("scale_factor", &FlatField::get_scale_factor,
                      &FlatField::set_scale_factor, R"(
                      float : scale factor to scale incident intensity to reasonable values (default 1.0)
                      )")

        .def_property("soft_window", &FlatField::get_soft_window,
                      &FlatField::set_soft_window, R"(
                      tuple of floats : bounds [degrees] to exclude strips with lower exposure (default (3.0, 34.0)).
                      )")

        .def_property(
            "bad_channels",
            [](FlatField &self) {
                auto bad_channels =
                    new NDArray<ssize_t, 1>(self.get_bad_channels());
                return return_image_data(bad_channels);
            },
            [](FlatField &self,
               py::array_t<bool, py::array::forcecast> bad_channels) {
                py::buffer_info info = bad_channels.request();
                if (info.ndim != 1 ||
                    info.format != py::format_descriptor<bool>::format()) {
                    throw std::runtime_error("Expected 1D buffer of type bool");
                }
                NDView<bool, 1> temp_array_view(
                    reinterpret_cast<bool *>(info.ptr),
                    std::array<ssize_t, 1>{info.shape[0]});
                NDArray<bool, 1> temp_array(temp_array_view); // first copy
                self.set_bad_channels(
                    temp_array); // second copy TODO im copying twice
            },
            R"(
            numpy.ndarray of bool, shape (n_channels,) : Each element is ``True`` if the channel is bad, otherwise ``False``.
            )")

        .def(
            "read_module_parameters_from_file",
            [](FlatField &self, const std::filesystem::path &filename) {
                self.read_module_parameters_from_file(filename);
            },
            R"(
             read module parameters from file
             (expects following format module [module_index] center [center] +- [error] conversion [conversion] +- [error] offset [offset] +- [error])

             Parameters
             ----------
             
             filename: str
                path to file containing module parameters
                )")

        .def(
            "read_bad_channels_from_file",
            [](FlatField &self, const std::filesystem::path &filename) {
                self.read_bad_channels_from_file(filename);
            },
            R"(
             read bad channels from file 

             Parameters
             ----------
             
             filename: str
                path to file containing bad channels
                )");
}