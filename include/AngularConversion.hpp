#include "FlatField.hpp"
#include "MythenDetector.hpp"
#include "aare/NDArray.hpp"

namespace angcal {

using namespace aare;

class AngularConversion {

  public:
    AngularConversion(std::shared_ptr<MythenDetector> mythen_detector_)
        : m_mythen_detector(mythen_detector_){};

    AngularConversion(std::shared_ptr<MythenDetector> mythen_detector_,
                      double histogram_bin_width_,
                      std::tuple<double, double> angular_range_)
        : m_mythen_detector(mythen_detector_),
          histogram_bin_width(histogram_bin_width_),
          m_min_angle(std::get<0>(angular_range_)),
          m_max_angle(std::get<1>(angular_range_)){};

    /** @brief calculate the corrected photon counts (flatfield, rate,
     * incident intensity)
     * @param photon_counts photon counts
     * @param global_strip_index strip index of photon counts
     * @param I0 incident intensity
     * @param exposure_time exposure time of acquisition [s]
     * @return pair {corrected photon counts, propagated_error}
     */
    std::pair<double, double>
    photon_count_correction(double photon_counts, const double exposure_time,
                            const uint64_t I0,
                            const size_t global_strip_index) const;

    /**
     * @brief correct and redistributes photon counts around designated region
     * to fixed angle width bins
     * @param module_index index of module
     * @param frame data from acquisition (storing detector position and photon
     * counts)
     * @param detector_angle detector axis motor position [degrees]
     * @param exposure_time exposure time [seconds]
     * @param incident_intensity incident intensity [photons/second]
     * @param weighted_moments accumulated weighted photon counts of order 0-2
     * per fixed angle width bin
     * @param num_contributing_photons accumulated number of photons
     * contributing to each fixed angle width bin
     */
    void process_photon_counts(const size_t module_index,
                               const NDView<uint32_t, 1> frame,
                               const double detector_angle,
                               const double exposure_time,
                               const uint64_t incident_intensity,
                               NDView<double, 2> weighted_moments,
                               NDView<size_t, 1> num_contributing_photons);

    ~AngularConversion() = default;

    /**
     * @brief convert frames of photon counts to a diffraction pattern
     * @param photon_counts photon counts of each frame (num_frames, num_strips)
     * @param detector_angles detector angle of each frame [degrees]
     * @param incident_intensities incident intensity of each frame
     * @param exposure_times exposure time of each frame [s]
     * @note the photon counts are corrected for flatfield, rate and incident
     * intensity and solid angle, redistributed to fixed angle width bins and
     * averaged
     * @return returns resulting diffraction pattern given in the range
     * [min_angle, max_angle] (values of -1.0 denote bins with insufficient
     * coverage or which are denoted by a bad channel) (first dimension
     * corresponds to redistributed photon counts, second dimension corresponds
     * to std of redistributed photon counts)
     */
    NDArray<double, 2> convert(const NDView<uint32_t, 2> photon_counts,
                               const NDView<double, 1> detector_angles,
                               const NDView<uint64_t, 1> incident_intensities,
                               const double exposure_time);

    NDArray<double, 2> convert(const NDView<uint32_t, 1> photon_counts,
                               const double detector_angle,
                               const uint64_t incident_intensity,
                               const double exposure_time);

    ssize_t num_fixed_angle_width_bins() const;

    double get_histogram_bin_width() const { return histogram_bin_width; }

    double get_scale_factor() const { return m_scale_factor; }

    void set_scale_factor(const double scale_factor_) {
        m_scale_factor = scale_factor_;
    }

    std::tuple<double, double> get_angular_range() const {
        return std::make_tuple(m_min_angle, m_max_angle);
    }

    void set_FlatField(NDView<const double, 2> flatfield_) {
        flatfield = NDArray<double, 2>(flatfield_);
    }

    NDView<double, 2> get_FlatField() { return flatfield.view(); }

    void set_bad_channels(const NDView<bool, 1> bad_channels_) {
        bad_channels = NDArray<bool, 1>(bad_channels_);

        // update disconnected modules
        for (size_t module_index = 0;
             module_index < m_mythen_detector->num_modules(); ++module_index) {
            if (std::all_of(bad_channels.begin() +
                                module_index *
                                    m_mythen_detector->strips_per_module(),
                            bad_channels.begin() +
                                (module_index + 1) *
                                    m_mythen_detector->strips_per_module(),
                            [](const auto &elem) { return elem; })) {
                m_mythen_detector->disconnected_modules()(module_index) = true;
            }
        }
    }

    NDView<bool, 1> get_bad_channels() { return bad_channels.view(); }

  private:
    NDArray<double, 2>
    computeErrorandMean(const NDView<double, 2> weighted_moments,
                        const NDView<size_t, 1> num_contributing_photons) const;

    std::shared_ptr<MythenDetector> m_mythen_detector;

    /**
     * bin width of fixed angle bin histogram [degrees]
     */
    double histogram_bin_width = 0.0036;

    /// @brief  @brief minimum angle for angular conversion [degrees]
    double m_min_angle = -180.0;

    /// @brief  maximum angle for angular conversion [degrees]
    double m_max_angle = 180.0;

    /**
     * @brief scale factor to scale everything to a reasonable scale
     * e.g. can be incident_intensity of first acquisition
     * @default 1.0
     */
    double m_scale_factor{1.0};

    /// @brief flatfield correction for each strip
    NDArray<double, 2> flatfield{};

    NDArray<bool, 1> bad_channels{}; // bad channels of detector
};

} // namespace angcal