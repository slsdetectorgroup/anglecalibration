#include "FlatField.hpp"
#include <pybind11/pybind11.h>

#include "np_helpers.hpp"

namespace py = pybind11;

using namespace angcal;

void define_FlatField_binding(py::module &m) {

    py::class_<FlatField, std::shared_ptr<FlatField>>(m, "FlatField")
        .def(py::init<std::shared_ptr<MythenDetectorSpecifications>,
                      std::optional<std::shared_ptr<SimpleFileInterface>>>(),
             py::arg("MythenDetectorSpecifications"),
             py::arg("file_interface") = std::nullopt)

        .def("create_flatfield_from_filelist",
             [](FlatField &self,
                std::shared_ptr<MythenFileReader> mythen_file_reader,
                const std::filesystem::path &file_list) {
                 self.create_flatfield_from_filelist(mythen_file_reader,
                                                     file_list);
             })

        .def("read_flatfield_from_file",
             [](FlatField &self, const std::string &filename) {
                 self.read_flatfield_from_file(filename);
             })

        .def_property(
            "flatfield",
            [](FlatField &self) {
                auto flatfield = new NDArray<double, 2>(self.get_flatfield());
                return return_image_data(flatfield);
            },
            [](FlatField &self, py::array flatfield) {
                py::buffer_info info = flatfield.request();
                if (info.ndim != 2 || info.strides[0] != sizeof(double) ||
                    info.itemsize != sizeof(double)) {
                    throw std::runtime_error(
                        "Expected 2D buffer of type double and stride 1");
                }
                NDView<double, 2> temp_array_view(
                    reinterpret_cast<double *>(info.ptr),
                    std::array<ssize_t, 2>{info.shape[0], info.shape[1]});
                NDArray temp_array(temp_array_view); // first copy
                self.set_flatfield(
                    temp_array); // second copy TODO im copying twice
            })

        .def_property(
            "inverse_normalized_flatfield",
            [](FlatField &self) {
                auto result = new NDArray<double, 2>(
                    self.get_inverse_normalized_flatfield());
                return return_image_data(result); // maybe return memory view
            },
            [](FlatField &self,
               py::array_t<double, py::array::c_style | py::array::forcecast>
                   &inverse_normalized_flatfield) {
                py::buffer_info info = inverse_normalized_flatfield.request();
                NDView<double, 2> temp_array_view(
                    reinterpret_cast<double *>(info.ptr),
                    std::array<ssize_t, 2>{info.shape[0], info.shape[1]});
                NDArray<double, 2> temp_array(
                    temp_array_view); // first copy //maybe modify buffer //or
                                      // iterate
                self.set_inverse_normalized_flatfield(
                    temp_array); // second copy TODO im copying twice
            })

        .def("calculate_normalized_flatfield",
             &FlatField::calculate_normalized_flatfield)

        .def("calculate_inverse_normalized_flatfield",
             &FlatField::calculate_inverse_normalized_flatfield)

        .def_property_readonly("get_normalized_flatfield", [](FlatField &self) {
            auto result =
                new NDArray<double, 2>(self.get_normalized_flatfield());
            return return_image_data(result);
        });
}