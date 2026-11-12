AngularConversion
=================

Utilities to convert raw photon counts from the Mythen detector into a diffraction pattern. As conversion parameters it expects the "historic detector group" DG parameters stored in the ``MythenDetector``.
See :ref:`parametersets` for a more detailed description of the different detector parameters.

The conversion consists of two building blocks:

* :py:class:`~angcal._angcal.HistogramMetaData` – describes the target histogram of the diffraction pattern (angular range and bin width).
* :py:func:`~angcal._angcal.convert` – corrects the photon counts of one or several frames and redistributes them to the histogram.

Converting frames: ``convert``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:py:func:`~angcal._angcal.convert` has two overloads:

* **single frame** – ``photon_counts`` of shape ``(n_channels,)`` together with a scalar ``detector_angle`` and ``incident_intensity``.
* **multiple frames** – ``photon_counts`` of shape ``(n_frames, n_channels)`` together with arrays ``detector_angles`` and ``incident_intensities`` of shape ``(n_frames,)``.

Because of the gaps between the modules, a single frame does not cover the full angular range. Usually several frames are acquired at slightly
shifted detector angles and converted together with the multi-frame overload, which averages all frames into a single diffraction pattern.

For every frame, every connected module and every strip the following steps are applied (see :ref:`theory <pearsonchisquare>` for the formulas):

1. **Bad channels** – strips marked ``True`` in ``bad_channels`` are skipped. Make sure the bad channels of the flatfield (value ``-1.0``) are included in this list.
2. **Disconnected modules** – modules marked as disconnected in the ``MythenDetector`` are skipped.
3. **Strip boundaries** – the left and right boundary angle of the strip are computed from the DG parameters, the ``detector_angle`` of the frame and the
   offsets stored in the ``MythenDetector``. Strips lying completely outside of ``angular_range`` are skipped.
4. **Photon count correction** – strips with zero photon counts contribute nothing. All other strips are corrected as follows, while the variance is propagated alongside:

   * **Mighell correction** – a constant of 1 is added to the photon counts.
   * **Rate correction** – the counts are corrected for the dead time of the detector (``MythenDetector.dead_time``) using the ``exposure_time``.
   * **Incident intensity correction** – the counts are scaled by ``scale_factor / incident_intensity``.
   * **Transverse width (solid angle) correction** – the counts are corrected for the angular height of the strip, which depends on the strip's distance to the sample.
   * **Flatfield correction** – the counts are divided by the normalized flatfield value of the strip, and the flatfield's standard deviation is propagated into the variance.

5. **Redistribution** – the corrected counts of a strip are distributed over all histogram bins the strip covers, scaled by ``bin_width / strip_width``.
   Each contribution is weighted with its inverse variance.
6. **Averaging** – after all frames are processed, every bin holds the inverse-variance weighted mean of all contributions. The error is computed from the
   reduced Pearson :math:`\chi^2` of the contributions, scaled by the inverse of the summed weights.

The result is a ``numpy.ndarray`` of shape ``(num_fixed_angle_width_bins, 2)``:

* column ``0`` – the corrected photon counts per bin,
* column ``1`` – the respective standard deviation.

Bins which were not covered by any good strip keep the value ``-1.0`` and should be masked before plotting or further analysis.

.. note::
    ``exposure_time`` and ``scale_factor`` are the same for all frames. 
    
.. note::
    ``bad_channels`` must have shape ``(n_channels,)`` and ``flatfield`` shape ``(n_channels, 2)``, where ``n_channels`` is the total number of strips
    of the ``MythenDetector``; otherwise a ``RuntimeError`` is raised. For the multi-frame overload, ``photon_counts``, ``detector_angles`` and
    ``incident_intensities`` must have the same number of frames.

**Example Usage:**

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from angcal import MythenDetector, HistogramMetaData, convert, EpicsMythenFileReader

    # (calibrated) DG parameters of the detector modules
    dg_parameters = np.loadtxt("dg_parameters.csv", delimiter=",", skiprows=1)

    mythen_detector = MythenDetector(dg_parameters)
    mythen_detector.elastic_correction_factor = 0.0045
    mythen_detector.detector_vertical_axis_offset = 24.8
    mythen_detector.offset = 0.0  # additional offset to the sample detector offset

    # histogram the photon counts are redistributed to
    histogram_meta_data = HistogramMetaData(angular_range=(0.0, 90.5), bin_width=0.0036)

    # normalized flatfield, shape (n_channels, 2), e.g. created with FlatField
    flatfield = np.load("normalized_flatfield.npy")

    # boolean array of shape (n_channels,), True marks a bad channel
    bad_channels = np.load("bad_channels.npy")
    bad_channels |= flatfield[:, 0] == -1.0  # include bad channels of the flatfield

    file_reader = EpicsMythenFileReader("sample_scan_0000.h5")

    # exposure time and scale factor are shared by all frames
    first_header = file_reader.read_header(0)
    exposure_time = first_header.exposure_time
    scale_factor = first_header.incident_intensity / exposure_time

    # --- single frame ---
    header, frame = file_reader.read_frame(0)
    photon_counts = frame.photon_counts[0]  # photon_counts has shape (counters, n_channels), use first counter
    single_frame_pattern = convert(mythen_detector, histogram_meta_data,
                                   photon_counts, header.detector_angle, header.incident_intensity,
                                   exposure_time, bad_channels, flatfield, scale_factor)

    # --- multiple frames at different detector angles (fills the module gaps) ---
    file_reader.seek(0)
    headers, frames = file_reader.read_n(4)  # frames: shape (4, n_channels)
    diffraction_pattern = convert(mythen_detector, histogram_meta_data,
                                  frames, headers["detector_angle"], headers["incident_intensity"],
                                  exposure_time, bad_channels, flatfield, scale_factor)


See ``python/examples/conversion.py`` for a complete script.

.. py:currentmodule:: angcal

.. autoclass:: angcal._angcal.HistogramMetaData
    :special-members: __init__
    :show-inheritance:
    :inherited-members:
    :member-order: groupwise


.. autofunction:: angcal._angcal.convert
