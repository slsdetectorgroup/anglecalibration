#include "MythenFileReader.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include "np_helpers.hpp"

namespace py = pybind11;

using namespace angcal;

void define_MythenFrame_bindings(py::module &m) {
    py::class_<MythenFrame>(m, "MythenFrame")

        .def_property_readonly(
            "photon_counts",
            [](MythenFrame &self) {
                return return_view_data(self.photon_counts(), py::cast(self));
            },
            R"(ndarray: Photon counts stored as 2D array with shape (counters, num_strips))")
        // TODO: handle view

        .def(
            "photon_count",
            [](MythenFrame &self, size_t row, size_t col = 0) {
                return self.photon_counts(row, col);
            },
            py::arg("row"), py::arg("col") = 0,
            R"(Get photon count for a specific strip index (row))")

        .def_property_readonly(
            "detector_angle",
            [](MythenFrame &self) { return self.detector_angle; },
            R"(float: Detector angle in degrees)")
        .def_property_readonly(
            "channel_mask", [](MythenFrame &self) { return self.channel_mask; },
            R"(list: Channel mask)")

        .def_property_readonly(
            "incident_intensity",
            [](MythenFrame &self) { return self.incident_intensity; },
            R"(int: Incident intensity)")

        .def_property_readonly(
            "exposure_time",
            [](MythenFrame &self) { return self.exposure_time; },
            R"(float: Exposure time in seconds)");
}

void define_MythenFileReader_bindings(py::module &m) {

    py::class_<MythenFileReader, std::shared_ptr<MythenFileReader>>(
        m, "MythenFileReader");

    py::class_<RawMythenFileReader, MythenFileReader,
               std::shared_ptr<RawMythenFileReader>>(m, "RawMythenFileReader")
        .def(py::init<const std::filesystem::path &,
                      const std::filesystem::path &>(),
             py::arg("detector_positions_filename"),
             py::arg("incident_intensities_filename"), R"(
             Parameters: 
                detector_positions_filename: str
                    path to file containing detector positions (in degrees) need to be stored as double in binary format. 
                incident_intensities_filename: str
                    path to file containing incident intensities (I0) values need to be stored as uint64_t in binary format.
        )")
        .def(
            "read_frame",
            [](RawMythenFileReader &self, const std::string &file_name) {
                return self.read_frame(file_name);
            },
            py::arg("file_name"), R"(
             read a MythenFrame from file.)")

        .def(
            "read_detector_angle",
            [](RawMythenFileReader &self, const std::string &file_name) {
                return self.read_detector_angle(file_name);
            },
            py::arg("file_name"), R"(read detector angle from file.)");

    py::class_<EpicsMythenFileReader, MythenFileReader,
               std::shared_ptr<EpicsMythenFileReader>>(m,
                                                       "EpicsMythenFileReader")
        .def(py::init<>())
        .def(py::init<const std::filesystem::path &>(),
             py::arg("incident_intensities_filename"), R"(
             Parameters:
                incident_intensities_filename: str
                    path to file containing incident intensities (I0) values need to be stored as uint64_t in binary format. 
        )")
        .def(
            "read_frame",
            [](EpicsMythenFileReader &self, const std::string &file_name) {
                return self.read_frame(file_name);
            },
            py::arg("file_name"), R"(
             read a MythenFrame from file.)")

        .def(
            "read_detector_angle",
            [](EpicsMythenFileReader &self, const std::string &file_name) {
                return self.read_detector_angle(file_name);
            },
            py::arg("file_name"), R"(read detector angle from file.)");
}