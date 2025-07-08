#include "helpers/FileInterface.hpp"
#include <pybind11/pybind11.h>

namespace py = pybind11;

using namespace angcal;


struct PySimpleFileInterface : SimpleFileInterface {
    // using SimpleFileInterface::SimpleFileInterface; // inherit constructors
    // if
    //  any

    // maybe inherit open - but how to implement read_into then?

    // Override virtual methods, forwarding to Python overrides
    void open(const std::string &filename) override {
        PYBIND11_OVERRIDE(void, SimpleFileInterface, open, filename);
    } // this is actually not a pure virtual function

    void read_into(std::byte *image_buf,
                   const ssize_t data_type_bytes = 1) override {
        PYBIND11_OVERRIDE_PURE(void, SimpleFileInterface, read_into, image_buf,
                               data_type_bytes);
    }
};

void define_SimpleFileInterface_binding(py::module &m) {
    py::class_<SimpleFileInterface, PySimpleFileInterface,
               std::shared_ptr<SimpleFileInterface>>(m, "SimpleFileInterface")
        .def(py::init<>())
        .def("open", [](SimpleFileInterface &self,
                        const std::string &filename) { self.open(filename); })
        .def("read_into", [](SimpleFileInterface &self, py::buffer buf,
                             const ssize_t data_type_bytes = 1) {
            py::buffer_info info = buf.request();
            if (info.itemsize != data_type_bytes) {
                throw std::runtime_error(
                    "Buffer itemsize doesn't match data_type_bytes");
            }

            self.read_into(reinterpret_cast<std::byte *>(info.ptr),
                           data_type_bytes);
        });
}