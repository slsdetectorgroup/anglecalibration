AngularConversion
=================

Utilities to convert raw photon counts from the Mythen detector into a diffraction pattern. As conversion parameters it expects the "historic detector group" DG parameters stored in the ``MythenDetector``.
See :ref:`parametersets` for a more detailed description of the different detector parameters.

The conversion is declared in ``AngularConversion.hpp`` and consists of two building blocks:

* :cpp:struct:`angcal::HistogramMetaData` – describes the target histogram of the diffraction pattern (angular range and bin width).
* :cpp:func:`angcal::convert` – corrects the photon counts of one or several frames and redistributes them to the histogram.

Converting frames: ``convert``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:cpp:func:`angcal::convert` has two overloads:

* **single frame** – ``photon_counts`` as ``NDView<uint32_t, 1>`` of shape ``(n_channels)`` together with a scalar ``double detector_angle`` and ``uint64_t incident_intensity``.
* **multiple frames** – ``photon_counts`` as ``NDView<uint32_t, 2>`` of shape ``(n_frames, n_channels)`` together with ``NDView<double, 1> detector_angles`` and
  ``NDView<uint64_t, 1> incident_intensities`` of shape ``(n_frames)``.

Because of the gaps between the modules, a single frame does not cover the full angular range. Usually several frames are acquired at slightly
shifted detector angles and converted together with the multi-frame overload, which averages all frames into a single diffraction pattern.

For every frame, every connected module and every strip the following steps are applied (see :ref:`theory <pearsonchisquare>` for the formulas):

1. **Bad channels** – strips marked ``true`` in ``bad_channels`` are skipped. Make sure the bad channels of the flatfield (value ``-1.0``) are included in this list.
2. **Disconnected modules** – modules marked as disconnected in the ``MythenDetector`` are skipped.
3. **Strip boundaries** – the left and right boundary angle of the strip are computed from the DG parameters, the ``detector_angle`` of the frame and the
   offsets stored in the ``MythenDetector``. Strips lying completely outside of ``angular_range`` are skipped.
4. **Photon count correction** – strips with zero photon counts contribute nothing. All other strips are corrected as follows, while the variance is propagated alongside:

   * **Mighell correction** – a constant of 1 is added to the photon counts.
   * **Rate correction** – the counts are corrected for the dead time of the detector (:cpp:func:`angcal::MythenDetector::dead_time`) using the ``exposure_time``.
   * **Incident intensity correction** – the counts are scaled by ``scale_factor / incident_intensity``.
   * **Transverse width (solid angle) correction** – the counts are corrected for the angular height of the strip, which depends on the strip's distance to the sample.
   * **Flatfield correction** – the counts are divided by the normalized flatfield value of the strip, and the flatfield's standard deviation is propagated into the variance.

5. **Redistribution** – the corrected counts of a strip are distributed over all histogram bins the strip covers, scaled by ``bin_width / strip_width``.
   Each contribution is weighted with its inverse variance.
6. **Averaging** – after all frames are processed, every bin holds the inverse-variance weighted mean of all contributions. The error is computed from the
   reduced Pearson :math:`\chi^2` of the contributions, scaled by the inverse of the summed weights.

The result is an ``NDArray<double, 2>`` of shape ``(num_fixed_angle_width_bins(), 2)``:

* column ``0`` – the corrected photon counts per bin,
* column ``1`` – the respective standard deviation.

Bins which were not covered by any good strip keep the value ``-1.0`` and should be masked before plotting or further analysis.

.. note::
    ``exposure_time`` and ``scale_factor`` are the same for all frames. ``scale_factor`` defaults to ``1.0``.

.. note::
    ``bad_channels`` must have shape ``(n_channels)`` and ``flatfield`` shape ``(n_channels, 2)``, where ``n_channels`` is the total number of strips
    of the ``MythenDetector`` (:cpp:func:`angcal::MythenDetector::num_strips`); otherwise a ``std::runtime_error`` is thrown. For the multi-frame overload,
    ``photon_counts``, ``detector_angles`` and ``incident_intensities`` must have the same number of frames.

**Example Usage:**

.. code-block:: cpp

    #include "AngularConversion.hpp"
    #include "MythenDetector.hpp"
    #include "MythenFileReader.hpp"
    #include "utils.hpp"

    using namespace angcal;

    int main() {
        const size_t num_modules = 48;

        // (calibrated) DG parameters of the detector modules
        DGParameters dgparameters(num_modules);
        read_csv("dg_parameters.csv", dgparameters.parameters.view(), true);

        MythenDetector mythen_detector(dgparameters);
        mythen_detector.elastic_correction_factor(0.0045);
        mythen_detector.detector_vertical_axis_offset(24.8);
        mythen_detector.offset(0.0); // additional offset to the sample detector offset

        // histogram the photon counts are redistributed to
        HistogramMetaData histogram_meta_data{std::make_tuple(0.0, 90.5), 0.0036};

        // normalized flatfield, shape (n_channels, 2), e.g. created with FlatField
        NDArray<double, 2> flatfield =
            read_flatfield("normalized_flatfield.raw", mythen_detector.num_modules());

        // shape (n_channels), true marks a bad channel
        NDArray<bool, 1> bad_channels =
            read_bad_channels("bad_channels.chans", mythen_detector.num_strips());
        for (ssize_t strip = 0; strip < bad_channels.size(); ++strip) {
            // include bad channels of the flatfield
            bad_channels[strip] = bad_channels[strip] || flatfield(strip, 0) == -1.0;
        }

        EpicsMythenFileReader file_reader("sample_scan_0000.h5");

        // exposure time and scale factor are shared by all frames
        EpicsMythenHeader first_header = file_reader.read_header(0);
        const double exposure_time = first_header.exposure_time;
        const double scale_factor = first_header.incident_intensity / exposure_time;

        // --- single frame ---
        auto [header, frame] = file_reader.read_frame(0);
        NDArray<double, 2> single_frame_pattern = convert(
            mythen_detector, histogram_meta_data, frame.photon_counts(),
            header.detector_angle, header.incident_intensity, exposure_time,
            bad_channels.view(), flatfield.view(), scale_factor);

        // --- multiple frames at different detector angles (fills the module gaps) ---
        file_reader.seek(0);
        auto [headers, frames] = file_reader.read_n(4); // frames: shape (4, n_channels)

        NDArray<double, 1> detector_angles(std::array<ssize_t, 1>{headers.size()});
        NDArray<uint64_t, 1> incident_intensities(std::array<ssize_t, 1>{headers.size()});
        for (ssize_t i = 0; i < headers.size(); ++i) {
            detector_angles[i] = headers[i].detector_angle;
            incident_intensities[i] = headers[i].incident_intensity;
        }

        NDArray<double, 2> diffraction_pattern = convert(
            mythen_detector, histogram_meta_data, frames.view(),
            detector_angles.view(), incident_intensities.view(), exposure_time,
            bad_channels.view(), flatfield.view(), scale_factor);
    }


See ``examples/example_conversion.cpp`` for a complete program.

..
    TODO: maybe move detailed description to another page.

.. doxygenstruct:: angcal::HistogramMetaData
   :members:
   :undoc-members:
   :private-members:


.. doxygenfunction:: angcal::convert(const MythenDetector&, const HistogramMetaData&, const NDView<uint32_t, 2>, const NDView<double, 1>, const NDView<uint64_t, 1>, const double, const NDView<const bool, 1>, const NDView<const double, 2>, const double)


.. doxygenfunction:: angcal::convert(const MythenDetector&, const HistogramMetaData&, const NDView<uint32_t, 1>, const double, const uint64_t, const double, const NDView<const bool, 1>, const NDView<const double, 2>, const double)
