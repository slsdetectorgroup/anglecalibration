/************************************************
 * @file MythenFileReader.test.cpp
 * @short test cases for the EpicsMythenFileReader
 ***********************************************/

#include "MythenFileReader.hpp"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <string>
#include <vector>

#include <H5Cpp.h>

#include "TempEpicsMythenFile.hpp"
#include "TempFileDirectory.hpp"
#include "test_config.hpp"

#include <catch2/catch_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace aare;
using namespace angcal;

namespace {

void check_header(const EpicsMythenHeader &header,
                  const DummyAcquisition &expected) {
    CHECK(header.detector_angle == expected.detector_angle);
    CHECK(header.incident_intensity ==
          static_cast<uint64_t>(std::lround(expected.izero)));
    CHECK(header.exposure_time == expected.exposure_time);
    CHECK(header.channel_mask ==
          counter_mask_to_channel_indices(
              static_cast<uint8_t>(expected.counter_mask)));
}

void check_frame(MythenFrame &frame, const DummyAcquisition &expected) {
    const auto n_strips = static_cast<ssize_t>(expected.photon_counts.size());
    REQUIRE(frame.size() == n_strips);
    auto photon_counts = frame.photon_counts();
    REQUIRE(photon_counts.size() == n_strips);
    for (ssize_t i = 0; i < n_strips; ++i) {
        CHECK(photon_counts(i) == expected.photon_counts[i]);
    }
}

DummyAcquisition make_acquisition(size_t file_index, size_t n_strips) {
    DummyAcquisition acquisition;
    acquisition.detector_angle = static_cast<double>(file_index * 10);
    acquisition.izero = static_cast<double>(file_index * 100);
    acquisition.exposure_time = 1.0 + static_cast<double>(file_index);
    acquisition.counter_mask = static_cast<uint32_t>(file_index % 8);
    acquisition.photon_counts.resize(n_strips);
    std::iota(acquisition.photon_counts.begin(),
              acquisition.photon_counts.end(), file_index * 1000);
    return acquisition;
}

} // namespace

TEST_CASE("counter mask to channel indices", "[mythenfilereader]") {
    CHECK(counter_mask_to_channel_indices(0) ==
          std::array<uint8_t, 3>{0, 0, 0});
    CHECK(counter_mask_to_channel_indices(1) ==
          std::array<uint8_t, 3>{1, 0, 0});
    CHECK(counter_mask_to_channel_indices(4) ==
          std::array<uint8_t, 3>{0, 0, 1});
    CHECK(counter_mask_to_channel_indices(7) ==
          std::array<uint8_t, 3>{1, 1, 1});
}

TEST_CASE("read single frame file with 1d dataset", "[mythenfilereader]") {
    constexpr size_t n_strips = 64;

    const auto acquisition = make_acquisition(3, n_strips);
    TempEpicsMythenFile temp_file("scan", {acquisition}, 1);

    const auto filename = temp_file.get_filename();

    EpicsMythenFileReader file_reader(filename);

    CHECK(file_reader.base_name() == "scan");
    CHECK(file_reader.n_files() == 1);
    CHECK(file_reader.total_frames() == 1);
    CHECK(file_reader.pixels_per_frame() == n_strips);
    CHECK(file_reader.bytes_per_frame() == n_strips * sizeof(uint32_t));
    CHECK(file_reader.current_file() == filename);
    CHECK(file_reader.tell() == 0);

    auto [header, frame] = file_reader.read_frame();

    check_header(header, acquisition);
    check_frame(frame, acquisition);
    CHECK(frame.photon_counts(5) == acquisition.photon_counts[5]);

    CHECK(file_reader.tell() == 1);
    CHECK_THROWS(file_reader.read_frame());
}

TEST_CASE("read series of single frame files", "[mythenfilereader]") {
    constexpr size_t n_strips = 32;

    size_t num_files = 5;
    size_t num_frames_per_file = 1;
    std::vector<DummyAcquisition> acquisitions;
    acquisitions.reserve(num_files);
    for (size_t file_index = 0; file_index < num_files; ++file_index) {
        acquisitions.push_back(make_acquisition(file_index, n_strips));
    }

    TempEpicsMythenFile temp_file("scan", acquisitions, num_frames_per_file);

    EpicsMythenFileReader file_reader(temp_file.get_filename());

    CHECK(file_reader.n_files() == num_files);
    CHECK(file_reader.total_frames() == num_files * num_frames_per_file);

    SECTION("read consecutively opens next file") {
        for (size_t i = 0; i < acquisitions.size(); ++i) {
            size_t file_number = i / num_frames_per_file;
            CHECK(file_reader.current_file() ==
                  temp_file.get_filesystem_dir() /
                      fmt::format("scan_{:04d}.h5", file_number));
            auto [header, frame] = file_reader.read_frame();
            check_header(header, acquisitions[i]);
            check_frame(frame, acquisitions[i]);
        }
        CHECK_THROWS(file_reader.read_frame());
    }

    SECTION("seek to frame in other file") {
        auto [header, frame] = file_reader.read_frame(2);
        check_header(header, acquisitions[2]);
        check_frame(frame, acquisitions[2]);

        auto [first_header, first_frame] = file_reader.read_frame(0);
        check_header(first_header, acquisitions[0]);
        check_frame(first_frame, acquisitions[0]);
        CHECK(file_reader.tell() == 1);
    }
}

TEST_CASE("read series of multi frame files", "[mythenfilereader]") {
    constexpr size_t n_strips = 16;

    // file 0: frames 0-2, file 1: frames 3-4, file 2: frames 5-8
    const size_t frames_per_file = 4;

    const size_t num_files = 3;

    std::vector<DummyAcquisition> acquisitions;
    acquisitions.reserve(num_files * frames_per_file);
    for (size_t file_index = 0; file_index < num_files * frames_per_file;
         ++file_index) {
        acquisitions.push_back(make_acquisition(file_index, n_strips));
    }

    TempEpicsMythenFile temp_file("multi", acquisitions, frames_per_file);

    EpicsMythenFileReader file_reader(temp_file.get_filename());

    CHECK(file_reader.n_files() == num_files);
    CHECK(file_reader.total_frames() == num_files * frames_per_file);
    CHECK(file_reader.pixels_per_frame() == n_strips);

    SECTION("read all frames consecutively") {
        for (size_t i = 0; i < acquisitions.size(); ++i) {
            CHECK(file_reader.tell() == i);
            auto [header, frame] = file_reader.read_frame();
            check_header(header, acquisitions[i]);
            check_frame(frame, acquisitions[i]);
        }
        CHECK_THROWS(file_reader.read_frame());
    }

    SECTION("read_n across file boundaries") {
        const size_t files_to_read = 6;
        const size_t start_frame = 1;
        file_reader.seek(start_frame);

        auto [headers, photon_counts] =
            file_reader.read_n(files_to_read); // frames 1-6

        REQUIRE(headers.size() == files_to_read);
        REQUIRE(photon_counts.shape(0) == static_cast<ssize_t>(files_to_read));
        REQUIRE(photon_counts.shape(1) == static_cast<ssize_t>(n_strips));

        for (ssize_t i = 0; i < static_cast<ssize_t>(files_to_read); ++i) {
            const auto &expected = acquisitions[static_cast<size_t>(i) + 1];
            check_header(headers(i), expected);
            for (ssize_t strip = 0; strip < static_cast<ssize_t>(n_strips);
                 ++strip) {
                CHECK(photon_counts(i, strip) ==
                      expected.photon_counts[static_cast<size_t>(strip)]);
            }
        }

        CHECK(file_reader.tell() == start_frame + files_to_read);
        size_t file_index = start_frame + files_to_read / frames_per_file - 1;
        CHECK(file_reader.current_file() ==
              temp_file.get_filesystem_dir() /
                  fmt::format("multi_{:04d}.h5", file_index));
    }

    SECTION("read_into raw buffer") {
        std::vector<uint32_t> buffer(2 * n_strips);
        std::array<EpicsMythenHeader, 2> headers{};
        file_reader.seek(2);
        file_reader.read_into(reinterpret_cast<std::byte *>(buffer.data()), 2,
                              headers.data());

        check_header(headers[0], acquisitions[2]);
        check_header(headers[1], acquisitions[3]);
        CHECK(buffer[0] == acquisitions[2].photon_counts[0]);
        CHECK(buffer[n_strips - 1] == acquisitions[2].photon_counts.back());
        CHECK(buffer[n_strips] == acquisitions[3].photon_counts[0]);
        CHECK(buffer.back() == acquisitions[3].photon_counts.back());
    }

    SECTION("random access with read_frame(index)") {
        for (size_t index : {7, 0, 4, 3, 8, 2}) {
            auto [header, frame] = file_reader.read_frame(index);
            check_header(header, acquisitions[index]);
            check_frame(frame, acquisitions[index]);
            CHECK(file_reader.tell() == index + 1);
        }
    }

    SECTION("read meta data without moving cursor") {
        file_reader.seek(1);

        // frame in a different file than the currently opened one
        check_header(file_reader.read_header(6), acquisitions[6]);
        CHECK(file_reader.read_detector_angle(4) ==
              acquisitions[4].detector_angle);
        CHECK(file_reader.read_incident_intensity(8) ==
              static_cast<uint64_t>(acquisitions[8].izero));

        CHECK(file_reader.tell() == 1);
        CHECK(file_reader.current_file() ==
              temp_file.get_filesystem_dir() / "multi_0000.h5");

        auto [header, frame] = file_reader.read_frame();
        check_header(header, acquisitions[1]);
        check_frame(frame, acquisitions[1]);
    }

    SECTION("out of range access throws") {
        CHECK_THROWS(file_reader.seek(20));
        CHECK_THROWS(file_reader.read_frame(20));
        CHECK_THROWS(file_reader.read_header(20));
        CHECK_THROWS(file_reader.read_detector_angle(100));
    }
}

TEST_CASE("test Epics mythenfilereader with measured data",
          "[mythenfilereader][.files]") {

    auto fpath = test_data_path();

    std::cout << "Test data path: " << fpath << std::endl;

    REQUIRE(std::filesystem::exists(fpath));

    EpicsMythenFileReader file_reader(fpath / "Fructose_0p2_60_0060.h5");

    auto [header, frame] = file_reader.read_frame();

    CHECK(frame.size() == 1280 * 48);
    CHECK(frame.photon_counts().size() == 61440);
    CHECK(header.detector_angle == Catch::Approx(-0.71505));
    CHECK(header.channel_mask == std::array<uint8_t, 3>{0, 0, 1});
    CHECK(header.incident_intensity == 75766547);
    CHECK(header.exposure_time == 15.0);
}
