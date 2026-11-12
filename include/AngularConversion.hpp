#pragma once

#include <cmath>
#include <tuple>
#include <utility>

#include "MythenDetector.hpp"
#include "aare/NDArray.hpp"

namespace angcal {

using namespace aare;

/**
 * @brief describes the fixed angle width histogram the photon counts are
 * redistributed to
 */
struct HistogramMetaData {
    /// @brief (min_angle, max_angle) angular range for conversion [degrees]
    std::tuple<double, double> angular_range{-180.0, 180.0};

    /// @brief bin width of fixed angle width histogram [degrees]
    double bin_width{0.0036};

    /// @brief number of bins in fixed angle width histogram
    ssize_t num_fixed_angle_width_bins() const {
        return static_cast<ssize_t>(
            std::ceil(std::get<1>(angular_range) / bin_width) -
            std::floor(std::get<0>(angular_range) / bin_width));
    }
};

/** @brief calculate the corrected photon counts (flatfield, rate,
 * incident intensity, transverse width)
 * @param mythen_detector detector geometry and specifications
 * @param flatfield normalized flatfield (num_strips, 2) first column
 * flatfield value, second column its standard deviation
 * @param scale_factor scale factor to scale everything to a reasonable scale
 * e.g. incident intensity of first acquisition
 * @param photon_counts photon counts
 * @param exposure_time exposure time of acquisition [s]
 * @param I0 incident intensity
 * @param global_strip_index strip index of photon counts
 * @return pair {corrected photon counts, propagated variance}
 */
std::pair<double, double>
photon_count_correction(const MythenDetector &mythen_detector,
                        const NDView<const double, 2> flatfield,
                        const double scale_factor, double photon_counts,
                        const double exposure_time, const uint64_t I0,
                        const size_t global_strip_index);

/**
 * @brief correct and redistributes photon counts of one module to fixed angle
 * width bins
 * @param mythen_detector detector geometry and specifications
 * @param histogram_meta_data angular range and bin width of fixed angle width
 * histogram
 * @param bad_channels bad channels of detector (num_strips,)
 * @param flatfield normalized flatfield (num_strips, 2)
 * @param scale_factor scale factor for incident intensity correction
 * @param module_index index of module
 * @param frame photon counts of acquisition (num_strips,)
 * @param detector_angle detector axis motor position [degrees]
 * @param exposure_time exposure time [seconds]
 * @param incident_intensity incident intensity [photons/second]
 * @param weighted_moments accumulated weighted photon counts of order 0-2
 * per fixed angle width bin
 * @param num_contributing_photons accumulated number of photons
 * contributing to each fixed angle width bin
 */
void process_photon_counts(const MythenDetector &mythen_detector,
                           const HistogramMetaData &histogram_meta_data,
                           const NDView<const bool, 1> bad_channels,
                           const NDView<const double, 2> flatfield,
                           const double scale_factor, const size_t module_index,
                           const NDView<uint32_t, 1> frame,
                           const double detector_angle,
                           const double exposure_time,
                           const uint64_t incident_intensity,
                           NDView<double, 2> weighted_moments,
                           NDView<size_t, 1> num_contributing_photons);

/**
 * @brief compute weighted mean and its standard deviation per fixed angle
 * width bin
 * @param weighted_moments accumulated weighted photon counts of order 0-2
 * (num_bins, 3)
 * @param num_contributing_photons number of photons contributing to each bin
 * @return (num_bins, 2) first column weighted mean, second column standard
 * deviation (-1.0 for bins without contributions)
 */
NDArray<double, 2>
computeErrorandMean(const NDView<double, 2> weighted_moments,
                    const NDView<size_t, 1> num_contributing_photons);

/**
 * @brief convert frames of photon counts to a diffraction pattern
 * @param mythen_detector detector geometry and specifications
 * @param histogram_meta_data angular range and bin width of fixed angle width
 * histogram
 * @param photon_counts photon counts of each frame (num_frames, num_strips)
 * @param detector_angles detector angle of each frame [degrees]
 * @param incident_intensities incident intensity of each frame
 * @param exposure_time exposure time of the frames [s]
 * @param bad_channels bad channels of detector (num_strips,)
 * @param flatfield normalized flatfield (num_strips, 2)
 * @param scale_factor scale factor to scale everything to a reasonable scale
 * e.g. incident intensity of first acquisition (default: 1.0)
 * @note the photon counts are corrected for flatfield, rate and incident
 * intensity and solid angle, redistributed to fixed angle width bins and
 * averaged
 * @return returns resulting diffraction pattern given in the range
 * [min_angle, max_angle] (values of -1.0 denote bins with insufficient
 * coverage or which are denoted by a bad channel) (first dimension
 * corresponds to redistributed photon counts, second dimension corresponds
 * to std of redistributed photon counts)
 */
NDArray<double, 2> convert(const MythenDetector &mythen_detector,
                           const HistogramMetaData &histogram_meta_data,
                           const NDView<uint32_t, 2> photon_counts,
                           const NDView<double, 1> detector_angles,
                           const NDView<uint64_t, 1> incident_intensities,
                           const double exposure_time,
                           const NDView<const bool, 1> bad_channels,
                           const NDView<const double, 2> flatfield,
                           const double scale_factor = 1.0);

/**
 * @brief convert a single frame of photon counts to a diffraction pattern
 * @see convert
 */
NDArray<double, 2> convert(const MythenDetector &mythen_detector,
                           const HistogramMetaData &histogram_meta_data,
                           const NDView<uint32_t, 1> photon_counts,
                           const double detector_angle,
                           const uint64_t incident_intensity,
                           const double exposure_time,
                           const NDView<const bool, 1> bad_channels,
                           const NDView<const double, 2> flatfield,
                           const double scale_factor = 1.0);

} // namespace angcal
