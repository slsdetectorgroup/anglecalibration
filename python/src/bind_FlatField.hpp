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

        .def("create_flatfield_from_rawfilesystem",
             [](FlatField &self, const std::filesystem::path &file_path) {
                 self.create_flatfield_from_rawfilesystem(file_path);
             })

        .def("create_flatfield_from_filelist",
             [](FlatField &self, const std::filesystem::path &file_list) {
                 self.create_flatfield_from_filelist(file_list);
             }) // TODO: Can i

        .def("read_flatfield_from_file",
             [](FlatField &self, const std::string &filename) {
                 self.read_flatfield_from_file(filename);
             })

        .def_property(
            "flatfield",
            [](FlatField &self) {
                auto flatfield = new NDArray<uint32_t, 1>(self.get_flatfield());
                return return_image_data(flatfield);
            },
            [](FlatField &self, py::array_t<uint32_t> flatfield) {
                py::buffer_info info = flatfield.request();
                if (info.ndim != 1) {
                    throw std::runtime_error(
                        "Expected 1D buffer of type uint32_t");
                }
                NDView<uint32_t, 1> temp_array_view(
                    reinterpret_cast<uint32_t *>(info.ptr),
                    std::array<ssize_t, 1>{info.shape[0]});
                NDArray temp_array(temp_array_view); // first copy
                self.set_flatfield(
                    temp_array); // second copy TODO im copying twice
            })

        .def("inverse_normalized_flatfield",
             [](FlatField &self) {
                 auto result = new NDArray<double, 1>(
                     self.inverse_normalized_flatfield());
                 return return_image_data(result);
             })

        .def("normalized_flatfield", [](FlatField &self) {
            auto result = new NDArray<double, 1>(self.normalized_flatfield());
            return return_image_data(result);
        });
}