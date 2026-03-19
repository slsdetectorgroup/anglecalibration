
/**
 * stores flatfield for angle calibration
 */

#pragma once
#include <cmath>
#include <cstdint>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>

#include "CustomFiles.hpp"
#include "MythenDetectorSpecifications.hpp"
#include "MythenFileReader.hpp"
#include "Parameters.hpp"
#include "aare/File.hpp"
#include "aare/NDArray.hpp"
#include "helpers/FileInterface.hpp"
#include "helpers/LogFiles.hpp"

using namespace aare;
namespace angcal {

class FlatField {

  public:
    /**
     * @tparam custom_file_ptr: Fileclass that inherits from
     * DetectorFileInterface class needs to overload read_frame and constructor
     * needs to have the signature Fileclass(std::string filename, ssize_t rows,
     * ssize_t cols, std::string reading_mode) rows, and cols define the
     * intended image size to read
     */
    FlatField(std::shared_ptr<MythenDetectorSpecifications> mythen_detector_)
        : mythen_detector(mythen_detector_),
          DGparameters(mythen_detector->max_modules) {

        bad_channels = NDArray<bool, 1>(
            std::array<ssize_t, 1>{mythen_detector->num_strips()},
            false); // default initialized to false
    }

    void set_scale_factor(const double scale_factor_) {
        scale_factor = scale_factor_;
    }

    double get_scale_factor() const { return scale_factor; }

    double diffraction_angle_from_DG_parameters(
        const size_t module_index, const double detector_angle,
        size_t strip_index, const double distance_to_strip) const {

        double offset = DGparameters.offsets(module_index);
        double center = DGparameters.centers(module_index);
        double conversion = DGparameters.conversions(module_index);

        strip_index =
            std::signbit(conversion)
                ? MythenDetectorSpecifications::strips_per_module -
                      strip_index - 1
                : strip_index; // TODO: are the values sored in reserve?

        return offset +
               180.0 / M_PI *
                   (center * std::abs(conversion) -
                    std::atan((center - (strip_index + distance_to_strip)) *
                              std::abs(conversion)));
    }

    double transverse_width_correction_factor(const size_t module_index,
                                              size_t strip_index) const {

        // convert to EE parameters
        const auto [normal_distance, module_center_distance, angle] =
            DGparameters.convert_to_EEParameters(module_index);

        /*
        strip_index = std::signbit(normal_distance)
                          ? MythenDetectorSpecifications::strips_per_module -
                                strip_index - 1
                          : strip_index;
        */

        const double distance_sample_pixel =
            std::sqrt(std::pow(normal_distance, 2) +
                      std::pow(module_center_distance -
                                   mythen_detector->pitch * strip_index,
                               2));

        const double transverse_width_correction_factor =
            (2 * std::atan(mythen_detector->transverse_width /
                           (2 * distance_sample_pixel)));

        return transverse_width_correction_factor;
    }

    std::pair<double, double> solid_angle_correction(
        const double photon_counts, const double photon_counts_error,
        const size_t module_index, const size_t strip_index) const {

        // convert to EE parameters
        const auto [normal_distance, module_center_distance, angle] =
            DGparameters.convert_to_EEParameters(module_index);

        const double distance_sample_pixel =
            std::sqrt(std::pow(normal_distance, 2) +
                      std::pow(module_center_distance -
                                   mythen_detector->pitch * strip_index,
                               2));

        // orthogonal projection of strip area onto plane normal to beam
        double projection_strip_area =
            mythen_detector->pitch * mythen_detector->transverse_width *
            normal_distance /
            distance_sample_pixel; // normal_distance/distance_sample_pixel is
                                   // the cosine of the angle between the strip
                                   // normal and the diffracted beam

        double solid_angle =
            projection_strip_area / std::pow(distance_sample_pixel, 2);

        // scale solid angle to unit
        solid_angle /= mythen_detector->average_solid_angle;

        double solid_angle_corrected_photon_counts =
            photon_counts / solid_angle;
        double solid_angle_corrected_photon_counts_error =
            photon_counts_error / std::pow(solid_angle, 2); // error propagation

        return std::pair(solid_angle_corrected_photon_counts,
                         solid_angle_corrected_photon_counts_error);
    }

    std::pair<double, double> Antonio_solid_angle_correction(
        const double photon_counts, const double photon_counts_error,
        const size_t module_index, const size_t strip_index) const {

        /*
        double solid_angle =
            projection_strip_area / std::pow(distance_sample_pixel, 2);
        */

        double transwidth_correction_factor =
            transverse_width_correction_factor(module_index, strip_index);

        double angular_strip_width =
            std::abs(diffraction_angle_from_DG_parameters(module_index, 0.0,
                                                          strip_index, 0.5) -
                     diffraction_angle_from_DG_parameters(module_index, 0.0,
                                                          strip_index, -0.5));

        double solid_angle =
            M_PI / 180.0 * angular_strip_width * transwidth_correction_factor;

        LOG(TLogLevel::logDEBUG) << fmt::format("module_index {}, strip {}",
                                                module_index, strip_index);

        size_t index =
            module_index * mythen_detector->strips_per_module + strip_index;

        SolidAngleCorrectionLogFile2.append(
            fmt::format("{},{},{},{}\n", index, solid_angle,
                        transwidth_correction_factor, angular_strip_width));

        // scale solid angle to unit
        solid_angle /= mythen_detector->average_solid_angle;

        double solid_angle_corrected_photon_counts =
            photon_counts / solid_angle;
        double solid_angle_corrected_photon_counts_error =
            photon_counts_error / std::pow(solid_angle, 2); // error propagation

        return std::pair(solid_angle_corrected_photon_counts,
                         solid_angle_corrected_photon_counts_error);
    }

    void create_flatfield_from_filelist(
        const std::vector<std::filesystem::path> &filelist,
        std::shared_ptr<MythenFileReader> file_reader) {

        flat_field = NDArray<double, 2>(
            std::array<ssize_t, 2>{
                static_cast<ssize_t>(mythen_detector->num_strips()), 2},
            0.0);

        NDArray<double, 1> sum_statistical_weights(
            std::array<ssize_t, 1>{mythen_detector->num_strips()}, 0.0);

        uint32_t mighells_correction_constant = 1;

        SolidAngleCorrectionLogFile2.append(
            "strip_index,solid_angle,transverse_width,angular_strip_width\n");

        for (const auto &file : filelist) {
            auto frame = file_reader->read_frame(file);

            uint64_t incident_intensity = frame.incident_intensity;

            LOG(TLogLevel::logDEBUG1)
                << fmt::format("incident intensity {}", incident_intensity);

            if (frame.incident_intensity == 0) {
                throw std::runtime_error(
                    fmt::format("incident intensity is zero for file {}. "
                                "Cannot apply I0 correction",
                                file.string()));
            }

            const double I0_correction_factor =
                scale_factor / static_cast<double>(incident_intensity);

            for (ssize_t strip_index = 0;
                 strip_index < mythen_detector->num_strips(); ++strip_index) {

                if (bad_channels(strip_index)) {
                    continue; // skip bad channels
                }

                // mighells correction
                double photon_counts =
                    frame.photon_counts(strip_index) +
                    std::min(mighells_correction_constant,
                             frame.photon_counts(strip_index));
                double variance_photon_counts =
                    frame.photon_counts(strip_index) +
                    std::min(mighells_correction_constant,
                             frame.photon_counts(strip_index));

                if (photon_counts == 0) {
                    bad_channels(strip_index) = true; // mark as bad channel
                    continue;                         // skip bad channels
                }

                // incident intensity correction
                photon_counts *= I0_correction_factor;
                variance_photon_counts *=
                    I0_correction_factor * I0_correction_factor;

                size_t module_index =
                    strip_index /
                    MythenDetectorSpecifications::strips_per_module;
                size_t local_strip_index =
                    strip_index %
                    MythenDetectorSpecifications::strips_per_module;

                // solid angle correction
                auto [efficiency_photon_counts,
                      variance_efficiency_photon_counts] =
                    Antonio_solid_angle_correction(
                        photon_counts, variance_photon_counts, module_index,
                        local_strip_index);

                // weighted average for flatfield
                flat_field(strip_index, 0) += efficiency_photon_counts /
                                              variance_efficiency_photon_counts;
                sum_statistical_weights(strip_index) +=
                    1.0 / variance_efficiency_photon_counts;
            }
        }

        size_t num_files = filelist.size();

        double sum_good_strips = 0.0; // for normalization
        size_t num_good_strips = 0;   // for normalization

        for (ssize_t strip_index = 0;
             strip_index < mythen_detector->num_strips(); ++strip_index) {
            if (sum_statistical_weights(strip_index) == 0 ||
                bad_channels(strip_index)) {
                flat_field(strip_index, 0) = 0; // bad channel
                flat_field(strip_index, 1) = 0; // bad channel
            } else {
                flat_field(strip_index, 0) /=
                    sum_statistical_weights(strip_index);
                flat_field(strip_index, 1) =
                    std::pow(num_files / sum_statistical_weights(strip_index),
                             2); // TODO: not what Antonio suggested but
                                 // how did he handle error propagation?
                sum_good_strips += flat_field(strip_index, 0);
                LOG(TLogLevel::logDEBUG1)
                    << fmt::format("strip {}: sum_good_strips  {}", strip_index,
                                   sum_good_strips);
                ++num_good_strips;
            }
        }

        // normalize flatfield to one
        double normalization_factor = sum_good_strips / num_good_strips;
        LOG(TLogLevel::logDEBUG) << fmt::format(
            "normalization factor for flatfield: {}", normalization_factor);
        for (ssize_t strip_index = 0;
             strip_index < mythen_detector->num_strips(); ++strip_index) {
            flat_field(strip_index, 0) /= normalization_factor;
        }
    }

    // TODO: maybe dont support reading at all
    // caveat user has to pass by value or move otherwise copied twice from
    // filesystem and then into member variable !!!!
    /**
     * @brief read normalized flatfield from file
     * @param filename name of file
     * @param file_reader custom file_reader to read flatfield file default:
     * ascii file with //TODO finish documentation
     */
    void read_normalized_flatfield_from_file(
        const std::string &filename,
        const std::shared_ptr<SimpleFileInterface> file_reader =
            std::make_shared<CustomFlatFieldFile>()) {
        flat_field = NDArray<double, 2>(std::array<ssize_t, 2>{
            static_cast<ssize_t>(mythen_detector->num_strips()), 2});
        file_reader->open(filename);
        file_reader->read_into(reinterpret_cast<std::byte *>(flat_field.data()),
                               8); // TODO: should be int though !!!

        file_reader->close();
    }

    void read_module_parameters_from_file(
        const std::filesystem::path &filename,
        const std::shared_ptr<SimpleFileInterface> file_reader =
            std::make_shared<InitialAngCalParametersFile>()) {

        file_reader->open(filename);
        file_reader->read_into(DGparameters.parameters.buffer(), 8);

        file_reader->close();
    }

    void read_bad_channels_from_file(
        const std::filesystem::path &filename,
        const std::shared_ptr<SimpleFileInterface> file_reader =
            std::make_shared<CustomBadChannelsFile>()) {

        file_reader->open(filename);
        file_reader->read_into(
            reinterpret_cast<std::byte *>(bad_channels.data()));
    }

    void set_bad_channels(const NDArray<bool, 1> &bad_channels_) {
        bad_channels = bad_channels_;
    }

    NDArray<bool, 1> get_bad_channels() const { return bad_channels; }

    /**
     * set flatfield from NDArray
     */
    void set_normalized_flatfield(const NDArray<double, 2> &flat_field_) {
        flat_field = flat_field_;
    }

    void set_normalized_flatfield(NDArray<double, 2> &&flat_field_) {
        flat_field = std::move(flat_field_);
    }

    // TODO: should it be view? or copy?
    NDArray<double, 2> get_normalized_flatfield() const { return flat_field; }

    NDView<double, 2> get_normalized_flatfield_view() const {
        return flat_field.view();
    }

  private:
    /// @brief scale_factor is used to scale incident intensity to ~ one
    double scale_factor{1.0};

    std::shared_ptr<MythenDetectorSpecifications> mythen_detector;

    DGParameters DGparameters{};

    NDArray<double, 2> flat_field;

    NDArray<bool, 1> bad_channels{};
};

} // namespace angcal