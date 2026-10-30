#include "MythenFileReader.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include "np_helpers.hpp"

namespace py = pybind11;

using namespace angcal;

void define_EpicsMythenHeader_bindings(py::module &m) {
    py::class_<EpicsMythenHeader>(m, "EpicsMythenHeader")
        .def_readonly("detector_angle", &EpicsMythenHeader::detector_angle,
                      R"(float: Detector angle in degrees)")
        .def_readonly("incident_intensity",
                      &EpicsMythenHeader::incident_intensity,
                      R"(int: Incident intensity)")
        .def_readonly("exposure_time", &EpicsMythenHeader::exposure_time,
                      R"(float: Exposure time in seconds)")
        .def_readonly("channel_mask", &EpicsMythenHeader::channel_mask,
                      R"(list: Channel mask)");

    PYBIND11_NUMPY_DTYPE(EpicsMythenHeader, detector_angle, incident_intensity,
                         exposure_time, channel_mask);
}

void define_MythenFrame_bindings(py::module &m) {
    py::class_<MythenFrame>(m, "MythenFrame")

        .def_property_readonly(
            "photon_counts",
            [](MythenFrame &self) {
                return return_view_data(self.photon_counts(), py::cast(self));
            },
            R"(ndarray: Photon counts stored as 2D array with shape (counters, num_strips))")

        .def(
            "photon_count",
            [](MythenFrame &self, size_t row, size_t col = 0) {
                return self.photon_counts(row, col);
            },
            py::arg("row"), py::arg("col") = 0,
            R"(Get photon count for a specific strip index (row))");
}

void define_MythenFileReader_bindings(py::module &m) {

    py::class_<EpicsMythenFileReader>(m, "EpicsMythenFileReader")
        .def(py::init<const std::filesystem::path &,
                      std::optional<const std::filesystem::path>>(),
             py::arg("fname"),
             py::arg("incident_intensities_filename") = std::nullopt, R"(
             Parameters:
                fname : path 
                    path to any hdf5 file of the series storing the acquisitions,
                    the series starts at this file
                incident_intensities_filename: str (optional)
                    path to file containing incident intensities (I0) values need to be stored as uint64_t in binary format. 
        )")
        .def(
            "read_frame",
            [](EpicsMythenFileReader &self, const size_t frame_index) {
                return self.read_frame(frame_index);
            },
            py::arg("frame_index"), R"(
            read a MythenFrame from file at the given frame index (current frame index advances to the next frame).

            Parameters:
            -----------

            frame_index : int
                index of the frame to read

            Returns:
            --------
            tuple : (EpicsMythenHeader, MythenFrame)
                tuple of EpicsMythenHeader and MythenFrame
            )")

        .def(
            "read_frame",
            [](EpicsMythenFileReader &self) { return self.read_frame(); },
            R"(read the frame at the current index from file.
            The current frame index advances to the next frame.

            Returns:
            --------

            tuple : (EpicsMythenHeader, MythenFrame)
                tuple of EpicsMythenHeader and MythenFrame
            )")

        .def(
            "read_n",
            [](EpicsMythenFileReader &self, const size_t n_frames) {
                auto [headers, frames] = self.read_n(n_frames);
                auto header_ptr =
                    new NDArray<EpicsMythenHeader, 1>(std::move(headers));
                auto frames_ptr = new NDArray<uint32_t, 2>(std::move(frames));
                return std::make_tuple(return_image_data(header_ptr),
                                       return_image_data(frames_ptr));
            },
            py::arg("n_frames"), R"(
            read n MythenFrames from file starting at the current frame index (current frame index advances to the next frame).

            Parameters:
            -----------
            n_frames : int
                number of frames to read

            Returns:
            --------
            tuple : (NDArray<EpicsMythenHeader, 1>, NDArray<uint32_t, 2>)
                tuple of NDArray of EpicsMythenHeader and NDArray of photon counts
            )")

        .def(
            "read_header",
            [](EpicsMythenFileReader &self, size_t frame_index) {
                return self.read_header(frame_index);
            },
            py::arg("frame_index"), R"(
            read header from file at the given frame index. (does not advance the current frame index)
            Parameters:
            -----------
            frame_index : int
                index of the frame to read
                
            Returns:
            --------
            EpicsMythenHeader : header of the frame at the given index
            )")

        .def(
            "read_detector_angle",
            [](EpicsMythenFileReader &self, size_t frame_index) {
                return self.read_detector_angle(frame_index);
            },
            py::arg("frame_index"), R"(
            read detector angle from file at the given frame index. (does not advance the current frame index)
            Parameters:
            -----------
            frame_index : int
                index of the frame to read

            Returns:
            --------
            float : detector angle at the given frame index [degrees]
            )")

        .def(
            "read_incident_intensity",
            [](EpicsMythenFileReader &self, size_t frame_index) {
                return self.read_incident_intensity(frame_index);
            },
            py::arg("frame_index"), R"(
            read incident intensity from file at the given frame index. (does not advance the current frame index)
            Parameters:
            -----------
            frame_index : int
                index of the frame to read

            Returns:
            --------
            float : incident intensity at the given frame index
            )")

        .def(
            "seek",
            [](EpicsMythenFileReader &self, size_t frame_index) {
                self.seek(frame_index);
            },
            py::arg("frame_index"),
            R"(seek to the given frame index (index over all files of the series))")

        .def(
            "tell", [](EpicsMythenFileReader &self) { return self.tell(); },
            R"(tell the current frame index (index over all files of the series))")

        .def(
            "total_frames",
            [](EpicsMythenFileReader &self) { return self.total_frames(); },
            R"(return the total number of frames in the series of files)")

        .def(
            "n_files",
            [](EpicsMythenFileReader &self) { return self.n_files(); },
            R"(return the number of files in the series)")

        .def(
            "pixels_per_frame",
            [](EpicsMythenFileReader &self) { return self.pixels_per_frame(); },
            R"(return the number of pixels per frame)")

        .def(
            "bytes_per_frame",
            [](EpicsMythenFileReader &self) { return self.bytes_per_frame(); },
            R"(return the number of bytes per frame)")

        .def(
            "current_file",
            [](EpicsMythenFileReader &self) { return self.current_file(); },
            R"(return the path of the current file in the series)")

        .def(
            "base_name",
            [](EpicsMythenFileReader &self) { return self.base_name(); },
            R"(return the base name of the series of files)");

    /*
    .def(
        "read_into",
        [](EpicsMythenFileReader &self,
           py::array_t<uint32_t, py::array::c_style> buffer,
           const std::string &file_name) {
            self.read_into(
                reinterpret_cast<std::byte *>(buffer.mutable_data()),
                file_name);
        },
        py::arg("buffer"), py::arg("file_name"), R"(
         read the entire content of the file into a byte buffer.)");
    */
}