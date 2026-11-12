Normalized FlatField
=====================

The ``FlatField`` class is responsible for storing the flatfield correction factors. It provides methods to create a normalized flatfield from a series of flatfield acquisitions,
e.g. a scan of the detector over a homogeneously illuminated region (fluorescence or amorphous scatterer).

A flatfield is created from multiple acquisitions from a flatfield scan. 
Each photon count of each scan is corrected by the incident intensity and the solid angle of the strip. The multiple scans are then accumulated and normalized such that the mean over all good strips is one. The resulting flatfield can be stored and reused later.

1. **Accumulate** – every acquisition of the flatfield scan is passed to :py:meth:`~angcal._angcal.FlatField.add_flatfield`. The photon counts of each strip are corrected and summed up.
2. **Normalize** – once all acquisitions are added, :py:meth:`~angcal._angcal.FlatField.normalize_flatfield` applies the solid angle correction (same for each scan) and normalizes the accumulated flatfield such that the mean over all good strips is one.

The result is accessible via :py:attr:`~angcal._angcal.FlatField.normalized_flatfield`.

Accumulating acquisitions: ``add_flatfield``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For each strip of the detector the following steps are applied to the photon counts of a single acquisition:

* **Bad channels** – strips marked in :py:attr:`~angcal._angcal.FlatField.bad_channels` are skipped.
* **Soft window** – the angular boundaries of the strip are computed from the detector angle and the DG parameters stored in the ``MythenDetector``. Strips lying completely outside of
  :py:attr:`~angcal._angcal.FlatField.soft_window` are skipped. Strips which are only partially inside the soft window are weighted by the fraction of the strip covered by the soft window.
  This excludes the low exposure regions at the edges of the scan.
* **Mighell correction** – a constant of 1 is added to non-zero photon counts.
* **Incident intensity correction** – the photon counts are scaled by ``scale_factor / incident_intensity``. The :py:attr:`~angcal._angcal.FlatField.scale_factor` should be chosen
  such that the scaled values are of a reasonable magnitude, e.g. the incident intensity of the first acquisition.
* **Error propagation** – assuming Poisson statistics, the variance of each strip is accumulated alongside the photon counts.

Normalizing: ``normalize_flatfield``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* **Solid angle correction** – each good strip is multiplied by ``average_solid_angle / solid_angle`` of the strip, to correct for the different solid angles covered by the strips.
* The standard deviation is computed from the accumulated variance.
* **Normalization** – flatfield values and standard deviations are divided by the mean flatfield value over all good strips, so the normalized flatfield has a mean of one.

Strips that are bad channels or were never inside the soft window keep the value ``-1.0``.

.. note::
    ``normalize_flatfield`` should be called exactly once, after all acquisitions are added.

If a flatfield was already computed and stored (e.g. with ``np.save``), it can be loaded and assigned directly to :py:attr:`~angcal._angcal.FlatField.normalized_flatfield` without recomputing it.

**Example Usage:**

.. code-block:: python

    import numpy as np
    from angcal import FlatField, MythenDetector, EpicsMythenFileReader

    # initial DG parameters of the detector modules
    dg_parameters = np.loadtxt("dg_parameters.csv", delimiter=",", skiprows=1)

    mythen_detector = MythenDetector(dg_parameters)

    flatfield = FlatField(mythen_detector)

    # boolean array of shape (n_channels,), True marks a bad channel
    flatfield.bad_channels = np.load("bad_channels.npy")

    # only consider strips within 3 and 30 degrees
    flatfield.soft_window = (3.0, 30.0)

    file_reader = EpicsMythenFileReader("flatfield_scan_0000.h5")

    # scale incident intensity to reasonable values
    flatfield.scale_factor = file_reader.read_header(0).incident_intensity

    # accumulate all acquisitions of the flatfield scan
    num_frames = 1601
    for frame_index in range(num_frames):
        header, frame = file_reader.read_frame(frame_index)
        flatfield.add_flatfield(frame.photon_counts, header.incident_intensity, header.detector_angle)

    # apply solid angle correction and normalize
    flatfield.normalize_flatfield()

    # shape (n_channels, 2): column 0 flatfield value, column 1 standard deviation
    normalized_flatfield = flatfield.normalized_flatfield

    # store and reuse the flatfield later
    np.save("normalized_flatfield.npy", normalized_flatfield)

    flatfield_reloaded = FlatField(mythen_detector)
    flatfield_reloaded.normalized_flatfield = np.load("normalized_flatfield.npy")


.. py:currentmodule:: angcal

.. autoclass:: angcal._angcal.FlatField
    :special-members: __init__
    :show-inheritance:
    :inherited-members:
    :member-order: groupwise
