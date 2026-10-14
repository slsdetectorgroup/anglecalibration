#pragma once
#include "aare/NDView.hpp"
#include "helpers/ErrorPropagation.hpp"

#include <cmath>

namespace angcal {

using namespace aare;

/**
 * @brief redistributes photon counts around designated region to
 * fixed angle width bins
 * @param left_strip_boundary_angle left boundary of strip [degrees]
 * @param right_strip_boundary_angle right boundary of strip [degrees]
 * @param histogram_bin_width width of fixed angle width bins [degrees]
 * @param roi region of interest to redistribute photon counts [degrees]
 * @param corrected_photon_count
 * @param inverse_corrected_photon_count_variance
 * @param weighted_moments weighted photon counts of order 0-2 eg. order k:
 * sum of photon_count_i^k*weight_i counts for fixed angle width bin and the
 * corresponding variance
 * @param num_contributing_photons number of photons contributing to the
 * fixed angle width bin
 */
inline void redistribute_photon_count_to_fixed_angle_width_bins(
    const double left_strip_boundary_angle,
    const double right_strip_boundary_angle, const double histogram_bin_width,
    std::tuple<double, double> roi, const double corrected_photon_count,
    const double inverse_corrected_photon_count_variance,
    NDView<double, 2> weighted_moments,
    NDView<size_t, 1> num_contributing_photons) {

    double strip_width_angle =
        std::abs(right_strip_boundary_angle -
                 left_strip_boundary_angle); // strip width in angles

    // redistribute to fixed angle width bins
    double photon_counts_per_bin =
        corrected_photon_count * histogram_bin_width / strip_width_angle;

    double inverse_photon_counts_variance_per_bin =
        inverse_corrected_photon_count_variance *
        std::pow(strip_width_angle / histogram_bin_width,
                 2); // Var(aX) = a^2 Var(X)

    if (photon_counts_per_bin < std::numeric_limits<double>::epsilon() ||
        inverse_photon_counts_variance_per_bin <
            std::numeric_limits<double>::epsilon()) {

        // bad_channels(global_strip_index) =
        // true; // mark channel as bad if no
        return;
    }

    ssize_t left_bin_index_covered_by_strip = static_cast<ssize_t>(
        (std::max(std::get<0>(roi), left_strip_boundary_angle) /
         histogram_bin_width));
    ssize_t right_bin_index_covered_by_strip = static_cast<ssize_t>(
        (std::min(std::get<1>(roi), right_strip_boundary_angle) /
         histogram_bin_width));

    size_t proper_bin_index =
        0; // the computed bin indices dont start at zero but are relative
           // to detector angle and strip index

    for (ssize_t bin_index = left_bin_index_covered_by_strip;
         bin_index <= right_bin_index_covered_by_strip; ++bin_index) {

        // some strips dont cover an entire bin or extend over multiple bins

        double bin_coverage =
            std::abs(std::min(right_strip_boundary_angle,
                              (bin_index + 1) * histogram_bin_width) -
                     std::max(left_strip_boundary_angle,
                              (bin_index)*histogram_bin_width));

        double bin_coverage_factor =
            bin_coverage / histogram_bin_width; // how much of the strip is
                                                // covered by the bin
        // convert to bin index
        proper_bin_index =
            bin_index - static_cast<ssize_t>(
                            std::get<0>(roi) /
                            histogram_bin_width); // bin index starts at zero

        double statistical_weight =
            bin_coverage_factor * inverse_photon_counts_variance_per_bin;

        // set to zero if actually covered by beam
        if (weighted_moments(proper_bin_index, 0) == -1.0) {
            weighted_moments(proper_bin_index, 0) = 0.0;
            weighted_moments(proper_bin_index, 1) = 0.0;
            weighted_moments(proper_bin_index, 2) = 0.0;
        }

        weighted_moments(proper_bin_index, 0) += statistical_weight;
        weighted_moments(proper_bin_index, 1) +=
            statistical_weight * photon_counts_per_bin;
        weighted_moments(proper_bin_index, 2) +=
            statistical_weight * std::pow(photon_counts_per_bin, 2);

        num_contributing_photons(proper_bin_index) +=
            statistical_weight > 0.0 ? 1.0 : 0.0; // count number of photons
                                                  // contributing to bin
    }
}

/**
 * @brief calculate the rate corrected photon counts taking into account the
 * dead time
 * @param photon_count photon counts
 * @param photon_count_error error of photon counts (generally variance
 * error is used)
 * @param exposure_time exposure time of acquisition [s]
 * @param dead_time dead time of the detector [s]
 * @return pair {corrected photon counts, propagated_error}
 */
inline std::pair<double, double>
rate_correction(const double photon_count, const double photon_count_error,
                const double exposure_time, const double dead_time) {

    if (exposure_time == 0) {
        throw std::runtime_error(LOCATION +
                                 "Exposure time is zero - cannot perform rate "
                                 "correction!");
    }

    const double maximum_count_rate =
        std::exp(-1); // theoretical maximum count rate for
                      // dead_time*measured_photon_counts_per_second

    double photon_counts_per_second = photon_count / exposure_time;

    photon_counts_per_second = std::min(
        maximum_count_rate,
        photon_counts_per_second *
            dead_time); // multiply with dead time - for numerical algorithm

    double error_photon_counts_per_second =
        photon_count_error * std::pow(dead_time / exposure_time, 2);

    // -actual_count_rate*dead_time = W -> -dead_time*measured_count_rate =
    // We^{W} -> W: Lambert W function

    // numerical calculation of inverse of Lambart W - function
    // actually calculate dead_time*count_rate =
    // dead_time*actual_count_rate*exp(-dead_time*actual_count_rate)
    double W_prev_iter = photon_counts_per_second; // initial guess
    double W_next_iter{};
    double propagated_error = error_photon_counts_per_second;
    bool method_converged = false;
    while (!method_converged) {
        W_next_iter = photon_counts_per_second * std::exp(W_prev_iter);
        method_converged = std::abs(1 - W_next_iter / W_prev_iter) <
                           10 * std::numeric_limits<double>::epsilon();
        W_prev_iter = W_next_iter;
        propagated_error =
            error_propagation(std::exp(W_prev_iter),
                              photon_counts_per_second * std::exp(W_prev_iter),
                              error_photon_counts_per_second, propagated_error);
    }

    double actual_count_rate = W_prev_iter; // actual count rate * dead_time

    double rate_correction_factor =
        actual_count_rate / photon_counts_per_second;

    double rate_corrected_photon_counts = photon_count * rate_correction_factor;

    // double rate_corrected_photon_counts_error =
    // error_propagation(rate_correction_factor, photon_count,
    // photon_count_error, error_rate_correction_factor);

    double rate_corrected_photon_counts_error =
        photon_count_error *
        std::pow(rate_correction_factor, 2); // TODO: In Antonios code like this

    return std::pair<double, double>{rate_corrected_photon_counts,
                                     rate_corrected_photon_counts_error};
}

/** @brief calculate the incident intensity corrected photon counts
 * @param photon_counts photon counts
 * @param photon_count_error error of photon counts (generally variance
 * error is used)
 * @param incident_intensity incident intensity
 * @param scale_factor scale factor to scale everything to a reasonable scale
 * @return pair {corrected photon counts, propagated_error}
 */
inline std::pair<double, double> incident_intensity_correction(
    const double photon_counts, const double photon_count_error,
    const uint64_t incident_intensity, const double scale_factor) {

    if (incident_intensity == 0) {
        throw std::runtime_error(LOCATION +
                                 "Incident intensity is zero - cannot perform "
                                 "incident intensity correction");
    }

    double I0_correction_factor = scale_factor / incident_intensity;

    double I0_corrected_photon_counts = photon_counts * I0_correction_factor;
    double I0_corrected_photon_counts_error =
        photon_count_error * std::pow(I0_correction_factor, 2);

    return std::pair<double, double>{I0_corrected_photon_counts,
                                     I0_corrected_photon_counts_error};
}

inline std::pair<double, double>
transverse_width_correction(const double photon_counts,
                            const double photon_counts_error,
                            const double transverse_width) {

    const double transverse_width_correction_factor = 1.0 / transverse_width;

    double transverse_width_normalized_photon_counts =
        photon_counts * transverse_width_correction_factor;
    double transverse_width_normalized_photon_counts_error =
        photon_counts_error * std::pow(transverse_width_correction_factor, 2);

    return std::pair(transverse_width_normalized_photon_counts,
                     transverse_width_normalized_photon_counts_error);
}

/**
 * @brief calculated the flatfield corrected photon counts and the error
 * @param photon_counts photon counts
 * @param photon_counts_error error of photon counts (generally variance
 * error is used)
 * @param flatfield_value flatfield value
 * @param flatfield_error error of flatfield value
 * @return pair {corrected photon counts, propagated_error}
 */
inline std::pair<double, double> flatfield_correction(
    const double photon_counts, const double photon_counts_error,
    const double flatfield_value, const double flatfield_error) {

    if (flatfield_value <=
        0.1) { // std::numeric_limits<double>::epsilon() in Antonios code 0.1

        return std::pair(0.0, 0.0);
    }

    // flatfield normalization
    const double flatfield_normalized_photon_counts =
        photon_counts / flatfield_value;

    const double normalized_flatfield_variance = std::pow(flatfield_error, 2);

    const double flatfield_normalized_photon_counts_error = error_propagation(
        1. / flatfield_value,
        -1.0 * photon_counts / std::pow(flatfield_value, 2),
        photon_counts_error,
        normalized_flatfield_variance); // Variance // TODO: change should be
                                        // squared per default always talk about
                                        // variance error

    return std::pair(flatfield_normalized_photon_counts,
                     flatfield_normalized_photon_counts_error);
}

} // namespace angcal
