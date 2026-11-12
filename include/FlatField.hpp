
/**
 * stores flatfield for angle calibration
 */

#pragma once

#include "MythenDetector.hpp"
#include "aare/NDArray.hpp"
#include <memory>
// #include "helpers/LogFiles.hpp"

using namespace aare;
namespace angcal {

class FlatField {

  public:
    /**
     * @brief Constructor for FlatField class
     * @param mythen_detector_ ptr to MythenDetector storing all
     mythen specific parameters
     */
    FlatField(std::shared_ptr<MythenDetector> mythen_detector_);

    /**
     * @brief set scale factor to scale incident intensity to reasonable values
     * e.g. can be set to incident intensity of first acquisition
     */
    void set_scale_factor(const double scale_factor_);

    /**
     * @brief get scale factor
     */
    double get_scale_factor() const;

    /**
     * @brief set soft window bounds [degrees] to exclude strips with lower
     * exposure.
     */
    void set_soft_window(const std::pair<double, double> soft_window_);

    /**
     * @brief get soft window bounds [degrees]
     */
    std::pair<double, double> get_soft_window() const;

    /**
     * @brief adds flatfield from single acquisition to the accumulated
     * flatfield
     * @param photon_counts photon counts for each strip
     * @param incident_intensity incident intensity for normalization
     * @param detector_angle position of detector [degrees]
     */
    void add_flatfield(const NDView<uint32_t, 1> photon_counts,
                       const uint32_t incident_intensity,
                       const double detector_angle);

    /**
     * @brief normalizes the accumulated flatfield
     */
    void normalize_flatfield();

    void set_bad_channels(const NDArray<bool, 1> &bad_channels_);

    NDArray<bool, 1> get_bad_channels() const;

    /**
     * @brief set flatfield from NDArray
     */
    void set_normalized_flatfield(const NDArray<double, 2> &flat_field_);

    void set_normalized_flatfield(NDArray<double, 2> &&flat_field_);

    // TODO: should it be view? or copy?
    NDArray<double, 2> get_normalized_flatfield() const;

    NDView<const double, 2> get_normalized_flatfield_view() const;

  private:
    /// @brief scale_factor is used to scale incident intensity to ~ one
    double scale_factor{1.0};

    /// @brief detector specifications of Mythen detector
    std::shared_ptr<MythenDetector> mythen_detector;

    /// @brief normalized flatfield stored as 2D array of shape (num_strips, 2)
    /// where the two values are normalized flatfield value and standard
    /// deviation
    /// value sof -1.0 denote strips which are not covered by any flatfield
    /// acquisition or are marked as bad.
    NDArray<double, 2> flat_field;

    /// @brief boolean array of size of total number of strips/channels in
    /// detector, stores 'TRUE' if strip is a bad channel otherwise 'FALSE'
    NDArray<bool, 1> bad_channels{};

    /// @brief  soft window bounds [degrees] to exclude strips with lower
    /// exposure.
    std::pair<double, double> soft_window{3.0, 34.0};
};

} // namespace angcal