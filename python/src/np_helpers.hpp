// TODO: maybe include in aare CMAKE

#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "aare/NDArray.hpp"
#include "aare/NDView.hpp"

namespace py = pybind11;
using namespace aare;

// Pass image data back to python as a numpy array
template <typename T, ssize_t Ndim>
py::array return_image_data(aare::NDArray<T, Ndim> *image) {

    py::capsule free_when_done(image, [](void *f) {
        aare::NDArray<T, Ndim> *foo =
            reinterpret_cast<aare::NDArray<T, Ndim> *>(f);
        delete foo;
    });

    return py::array_t<T>(
        image->shape(),        // shape
        image->byte_strides(), // C-style contiguous strides for double
        image->data(),         // the data pointer
        free_when_done);       // numpy array references this parent
}

template <class T, int Flags> auto make_view_1d(py::array_t<T, Flags> &arr) {
    return aare::NDView<T, 1>(arr.mutable_data(), aare::Shape<1>(arr.shape(0)));
}

template <typename T, ssize_t Ndim>
py::array return_view_data(aare::NDView<T, Ndim> view) {

    auto byte_strides = view.strides();
    std::for_each(byte_strides.begin(), byte_strides.end(),
                  [](auto &stride) { stride *= sizeof(T); });

    return py::array_t<T>(view.shape(), // shape
                          byte_strides, // C-style contiguous strides for double
                          view.data()); // the data pointer
}