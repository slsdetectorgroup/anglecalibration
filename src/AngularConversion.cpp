#include "AngularConversion.hpp"
#include "PhotonCountProcessing.hpp"

namespace angcal {

ssize_t AngularConversion::num_fixed_angle_width_bins() const {
    ssize_t num_fixed_angle_width_bins =
        std::ceil(m_max_angle / histogram_bin_width) -
        std::floor(m_min_angle / histogram_bin_width);
    return num_fixed_angle_width_bins;
}

inline std::pair<double, double> AngularConversion::photon_count_correction(
    double photon_counts, const double exposure_time, const uint64_t I0,
    const size_t global_strip_index) const {

    if (photon_counts < std::numeric_limits<double>::epsilon()) {
        return std::pair(0.0, 0.0); // sanity check
    }

    // mighells statistics
    photon_counts += 1;

    auto [rate_corrected_photon_counts, rate_corrected_photon_counts_error] =
        rate_correction(
            photon_counts, photon_counts, exposure_time,
            m_mythen_detector
                ->dead_time()); // same variance as poisson distributed
                                // - maybe sqaured?

    auto [incident_intensity_corrected_photon_counts,
          incident_intensity_corrected_photon_counts_variance] =
        incident_intensity_correction(rate_corrected_photon_counts,
                                      rate_corrected_photon_counts_error, I0,
                                      m_scale_factor);

    const size_t module_index =
        global_strip_index / MythenDetector::strips_per_module();
    const size_t strip_index =
        global_strip_index % MythenDetector::strips_per_module();

    const double transverse_width =
        m_mythen_detector->transverse_width(module_index, strip_index, false);

    auto [transverse_width_corrected_photon_counts,
          transverse_width_corrected_photon_counts_variance] =
        transverse_width_correction(
            incident_intensity_corrected_photon_counts,
            incident_intensity_corrected_photon_counts_variance,
            transverse_width);

    auto [flatfield_normalized_photon_counts,
          flatfield_normalized_photon_counts_variance] =
        flatfield_correction(
            transverse_width_corrected_photon_counts,
            transverse_width_corrected_photon_counts_variance,
            flatfield(global_strip_index, 0),
            flatfield(global_strip_index,
                      1)); // make sure to include bad chanenls of
                           // flatfield into bad channel list

    return std::pair<double, double>{
        flatfield_normalized_photon_counts,
        flatfield_normalized_photon_counts_variance};
}

void AngularConversion::process_photon_counts(
    const size_t module_index, const NDView<uint32_t, 1> frame,
    const double detector_angle, const double exposure_time,
    const uint64_t incident_intensity, NDView<double, 2> weighted_moments,
    NDView<size_t, 1> num_contributing_photons) {

    double left_boundary_roi_base_peak = m_min_angle;  // dummy values
    double right_boundary_roi_base_peak = m_max_angle; // dummy values

    // iterate over all strips of module
    for (ssize_t strip_index = 0;
         strip_index < MythenDetector::strips_per_module(); ++strip_index) {

        size_t global_strip_index =
            module_index * MythenDetector::strips_per_module() + strip_index;

        if (bad_channels(global_strip_index)) {
            continue; // skip bad channels
        }

        double left_strip_boundary_angle =
            m_mythen_detector->diffraction_angle_from_DG_parameters(
                module_index, detector_angle, strip_index,
                0.0); // left strip boundary in angles [degrees]

        double right_strip_boundary_angle =
            m_mythen_detector->diffraction_angle_from_DG_parameters(
                module_index, detector_angle, strip_index,
                +1.0); // right strip boundary in angles [degrees]

        // skip strip if not in ROI
        if (left_strip_boundary_angle > m_max_angle ||
            right_strip_boundary_angle < m_min_angle) {
            LOG(TLogLevel::logDEBUG1) << fmt::format(
                "strip {} of module {} covers an angle of [{}, {}] and is "
                "outside "
                "of the given angular range of [{}, "
                "{}] and will be skipped",
                strip_index, module_index, left_strip_boundary_angle,
                right_strip_boundary_angle, m_min_angle, m_max_angle);
            continue;
        }

        // get flatfield corrected photon counts and propagate error
        auto [corrected_photon_counts, corrected_photon_counts_variance] =
            photon_count_correction(frame(global_strip_index), exposure_time,
                                    incident_intensity, global_strip_index);

        double inverse_corrected_photon_counts_variance =
            corrected_photon_counts_variance <
                    std::numeric_limits<double>::epsilon()
                ? 0.0
                : 1.0 / corrected_photon_counts_variance;

        redistribute_photon_count_to_fixed_angle_width_bins(
            left_strip_boundary_angle, right_strip_boundary_angle,
            histogram_bin_width,
            std::make_tuple(left_boundary_roi_base_peak,
                            right_boundary_roi_base_peak),
            corrected_photon_counts, inverse_corrected_photon_counts_variance,
            weighted_moments, num_contributing_photons);
    }
}

NDArray<double, 2>
AngularConversion::convert(const NDView<uint32_t, 2> photon_counts,
                           const NDView<double, 1> detector_angles,
                           const NDView<uint64_t, 1> incident_intensities,
                           const double exposure_time) {

    if (photon_counts.shape(0) != detector_angles.shape(0) ||
        photon_counts.shape(0) != incident_intensities.shape(0)) {
        throw std::runtime_error(
            "Input arrays must have the same number of frames");
    }

    const ssize_t new_num_bins = num_fixed_angle_width_bins();

    // weighted photon counts of order 0-2 eg. order k: sum of
    // photon_count_i^k*weight_i counts for fixed angle width bin
    NDArray<double, 2> weighted_moments(std::array<ssize_t, 2>{new_num_bins, 3},
                                        -1.0);

    // number of photons contributing to the fixed angle width bin
    NDArray<size_t, 1> num_contributing_photons(
        std::array<ssize_t, 1>{new_num_bins}, 0);

    for (ssize_t frame_index = 0; frame_index < photon_counts.shape(0);
         ++frame_index) {

        auto frame_photon_counts = photon_counts.sub_view(frame_index);
        for (size_t module_index = 0;
             module_index < m_mythen_detector->num_modules(); ++module_index) {

            if (m_mythen_detector->disconnected_modules()(module_index)) {
                continue;
            }

            process_photon_counts(
                module_index, frame_photon_counts, detector_angles(frame_index),
                exposure_time, incident_intensities(frame_index),
                weighted_moments.view(), num_contributing_photons.view());
        }
    }

    return computeErrorandMean(weighted_moments.view(),
                               num_contributing_photons.view());
}

NDArray<double, 2> AngularConversion::convert(
    const NDView<uint32_t, 1> photon_counts, const double detector_angle,
    const uint64_t incident_intensity, const double exposure_time) {

    const ssize_t new_num_bins = num_fixed_angle_width_bins();

    // weighted photon counts of order 0-2 eg. order k: sum of
    // photon_count_i^k*weight_i counts for fixed angle width bin
    NDArray<double, 2> weighted_moments(std::array<ssize_t, 2>{new_num_bins, 3},
                                        -1.0);

    // number of photons contributing to the fixed angle width bin
    NDArray<size_t, 1> num_contributing_photons(
        std::array<ssize_t, 1>{new_num_bins}, 0);

    for (size_t module_index = 0;
         module_index < m_mythen_detector->num_modules(); ++module_index) {

        if (m_mythen_detector->disconnected_modules()(module_index)) {
            continue;
        }

        process_photon_counts(module_index, photon_counts, detector_angle,
                              exposure_time, incident_intensity,
                              weighted_moments.view(),
                              num_contributing_photons.view());
    }

    return computeErrorandMean(weighted_moments.view(),
                               num_contributing_photons.view());
}

NDArray<double, 2> AngularConversion::computeErrorandMean(
    const NDView<double, 2> weighted_moments,
    const NDView<size_t, 1> num_contributing_photons) const {

    const ssize_t new_num_bins = num_fixed_angle_width_bins();

    NDArray<double, 2> fixed_angle_width_bins_photon_counts(
        std::array<ssize_t, 2>{new_num_bins, 2},
        -1.0); // first dimension photon counts, second dimension error

    // divide by statistial weight
    for (ssize_t i = 0; i < weighted_moments.shape(0); ++i) {

        if (weighted_moments(i, 1) > -1.0 && num_contributing_photons(i) > 0) {
            fixed_angle_width_bins_photon_counts(i, 0) =
                weighted_moments(i, 0) < std::numeric_limits<double>::epsilon()
                    ? 0.0
                    : weighted_moments(i, 1) /
                          weighted_moments(i, 0); // weighted mean

            double inverse_weights = 1. / weighted_moments(i, 0);

            double goodness_of_fit =
                num_contributing_photons(i) == 1
                    ? 0.0
                    : (weighted_moments(i, 2) - weighted_moments(i, 1) *
                                                    weighted_moments(i, 1) *
                                                    inverse_weights) /
                          (num_contributing_photons(i) -
                           1); // calculates chi value for optimal

            fixed_angle_width_bins_photon_counts(i, 1) = std::sqrt(
                goodness_of_fit * inverse_weights); // scale by variance
        }
    }

    return fixed_angle_width_bins_photon_counts;
}

} // namespace angcal