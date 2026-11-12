AngleCalibration
=================

``AngleCalibration`` calibrates the module parameters of the Mythen detector, i.e. the parameters used to compute the diffraction angle of each strip.
Internally the calibration is done on the "Best computing" BC parameters. Initial parameters are passed in as "historic detector group" DG parameters when the
``MythenDetector`` is constructed, and the calibrated result is available as DG parameters, the set :py:func:`~angcal._angcal.convert` expects.
See :ref:`parametersets` for a more detailed description of the different detector parameters and :ref:`parametercalibration` for the theory behind the calibration.

How the calibration works
_________________________

The calibration uses a single peak of a known sample, the **base peak**, measured by many acquisitions taken at slightly shifted detector positions
(e.g. a scan where the detector is rotated around the sample). Ideally the base peak shows up at the same diffraction angle (corrected by the detector position) in every acquisition. Deviations reveal errors in the module parameters.

For each module the calibration runs in two steps:

1. **Coupled parameters** – the module center distance :math:`L_m` and the angle :math:`\delta_m` are optimized so that the base peaks of all
   acquisitions *within one module* overlap as well as possible (minimizing a variance weighted :math:`\chi^2` similarity).
2. **Offset** – the angle :math:`\Psi_m` between the module center and the beam is optimized so that the base peak of module :math:`m` aligns with the
   base peak of the previous module :math:`m-1`. Module 0 serves as the reference and its offset stays fixed.

Only acquisitions whose base peak region of interest ``[base_peak_angle - base_peak_ROI_width, base_peak_angle + base_peak_ROI_width]`` lies fully inside the module are used.

.. note::
    The calibration updates the "best computing" parameters of the ``MythenDetector`` passed to the constructor **in place**. However, the "historic detector group" parameters are not changed. The calibrated DG parameters can be accessed via :py:attr:`~angcal._angcal.AngleCalibration.DGparameters` and used to construct a new ``MythenDetector`` for conversion.

What you need
_____________

* **Initial DG parameters** – a reasonable starting guess for every module, e.g. a previous calibration. They are passed to the ``MythenDetector``.
* **A normalized flatfield** – see :doc:`pyFlatField`.
* **Bad channels** – a boolean ``numpy.ndarray`` of shape ``(n_channels,)`` where ``True`` marks a bad strip. Disconnected modules (all channels bad) cannot be calibrated.
* **A calibration scan** – frames of photon counts of shape ``(n_frames, n_channels)`` and the corresponding detector angles of shape ``(n_frames,)``.
  The frames **must be ordered by detector angle**. The scan has to be wide enough that every module sees the base peak in several acquisitions.
* **A base peak angle** – a well isolated, sharp peak of the sample. Pick one that every module can reach within the detector rotation range,
  for the current Mythen setup at SLS a peak between :math:`6^{\circ}` and :math:`33^{\circ}`.

Tunable settings
________________

* :py:attr:`~angcal._angcal.AngleCalibration.base_peak_ROI_width` – half width of the region around the base peak used for the comparison (default ``0.05°``).
  It should cover the whole peak including some background on both sides, but no neighbouring peaks.
* :py:attr:`~angcal._angcal.AngleCalibration.histogram_bin_width` – bin width of the fixed angle width histogram the photon counts are redistributed to (default ``0.0036°``).
* :py:attr:`~angcal._angcal.AngleCalibration.angular_range` – angular range of the histogram (degrees).

Example Usage for Calibration
_____________________________

.. code-block:: python

    import numpy as np
    from angcal import (MythenDetector, FlatField, AngleCalibration, HistogramMetaData,
                        convert, EpicsMythenFileReader, PlotHelper)

    # initial DG parameters, one row per module: [center, conversion, offset]
    dg_parameters = np.loadtxt("initial_dg_parameters.csv", delimiter=",", skiprows=1)

    # detector setup - holds the module parameters that get calibrated
    mythen_detector = MythenDetector(dg_parameters)

    # normalized flatfield, shape (n_channels, 2)
    flatfield = FlatField(mythen_detector)
    flatfield.normalized_flatfield = np.load("normalized_flatfield.npy")

    anglecalibration = AngleCalibration(mythen_detector, flatfield)

    # boolean array of shape (n_channels,), True marks a bad channel
    anglecalibration.bad_channels = np.load("bad_channels.npy")

    # calibration settings
    anglecalibration.histogram_bin_width = 0.0036  # degrees (default)
    anglecalibration.base_peak_ROI_width = 0.18    # degrees, half width of the base peak ROI
    anglecalibration.base_peak_angle = 19.0678     # degrees, center of the selected base peak

    # read the calibration scan (frames ordered by detector angle)
    file_reader = EpicsMythenFileReader("calibration_scan_0000.h5")
    num_frames = 1501
    headers, frames = file_reader.read_n(num_frames)
    detector_angles = headers["detector_angle"]

    # calibrate a single module (useful to check the settings, with plots of the fit)
    anglecalibration.calibrate(frames, detector_angles, anglecalibration.base_peak_angle,
                               module_index=1, plot_calibration_process=True)

    # calibrate all modules
    anglecalibration.calibrate(frames, detector_angles, anglecalibration.base_peak_angle)

    # calibrated DG parameters, shape (n_modules, 3): [center, conversion, offset]
    calibrated_parameters = anglecalibration.DGparameters
    np.savetxt("calibrated_dg_parameters.csv", calibrated_parameters, delimiter=",",
               header="center, conversion, offset", comments="", fmt="%.6f")

The calibrated parameters can be loaded into a new ``MythenDetector`` later and used for conversion, see :doc:`pyAngularConversion`.

API
___

.. py:currentmodule:: angcal

.. autoclass:: angcal._angcal.AngleCalibration
    :special-members: __init__
    :show-inheritance:
    :inherited-members:
    :member-order: groupwise
