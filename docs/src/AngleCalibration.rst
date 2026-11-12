
AngleCalibration
==================

``AngleCalibration`` calibrates the module parameters of the Mythen detector, i.e. the parameters used to compute the diffraction angle of each strip.
Internally the calibration is done on the "Best computing" BC parameters. Initial parameters are passed in as "historic detector group" DG parameters when the
``MythenDetector`` is constructed, and the calibrated result is available as DG parameters, the set :cpp:func:`angcal::convert` expects.
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
    The calibration updates the "best computing" parameters of the ``MythenDetector`` passed to the constructor **in place**. However, the "historic detector group" parameters are not changed. The calibrated DG parameters can be accessed via :cpp:func:`angcal::AngleCalibration::get_DGparameters` and used to construct a new ``MythenDetector`` for conversion.

What you need
_____________

* **Initial DG parameters** – a reasonable starting guess for every module, e.g. a previous calibration. They are passed to the ``MythenDetector``.
* **A normalized flatfield** – see :doc:`FlatField`.
* **Bad channels** – a boolean ``NDView<bool, 1>`` of size ``num_strips`` where ``true`` marks a bad strip. Disconnected modules (all channels bad) cannot be calibrated.
* **A calibration scan** – frames of photon counts ``NDView<uint32_t, 2>`` of shape ``(n_frames, n_channels)`` and the corresponding detector angles ``NDView<double, 1>`` of shape ``(n_frames,)``.
  The frames **must be ordered by detector angle**. The scan has to be wide enough that every module sees the base peak in several acquisitions.
* **A base peak angle** – a well isolated, sharp peak of the sample. Pick one that every module can reach within the detector rotation range,
  for the current Mythen setup at SLS a peak between :math:`6^{\circ}` and :math:`33^{\circ}`.

Tunable settings
________________

* :cpp:func:`angcal::AngleCalibration::set_base_peak_ROI_width` – half width of the region around the base peak used for the comparison (default ``0.05°``).
  It should cover the whole peak including some background on both sides, but no neighbouring peaks.
* :cpp:func:`angcal::AngleCalibration::set_histogram_bin_width` – bin width of the fixed angle width histogram the photon counts are redistributed to (default ``0.0036°``).
* :cpp:func:`angcal::AngleCalibration::set_angular_range` – angular range of the histogram (degrees).

Example Usage for Calibration
_____________________________

The helper functions ``read_csv``, ``read_flatfield`` and ``read_bad_channels`` are defined in ``utils.hpp``. A complete example can be found in ``examples/example_calibration.cpp``.

.. code-block:: cpp

    #include "AngleCalibration.hpp"
    #include "utils.hpp"

    #include <fstream>

    using namespace angcal;

    const size_t num_modules = 48;

    // initial DG parameters, one row per module: [center, conversion, offset]
    DGParameters dgparameters(num_modules);
    read_csv("initial_dg_parameters.csv", dgparameters.parameters.view(), true);

    // detector setup - holds the module parameters that get calibrated
    auto mythen_detector_ptr = std::make_shared<MythenDetector>(dgparameters);

    // normalized flatfield, shape (n_channels, 2)
    auto flat_field_ptr = std::make_shared<FlatField>(mythen_detector_ptr);
    flat_field_ptr->set_normalized_flatfield(
        read_flatfield("normalized_flatfield.raw", num_modules));

    AngleCalibration anglecalibration(mythen_detector_ptr, flat_field_ptr);

    // boolean array of size num_strips, true marks a bad channel
    auto bad_channels = read_bad_channels("bad_channels.chans",
                                          mythen_detector_ptr->num_strips());
    anglecalibration.set_bad_channels(bad_channels.view());

    // calibration settings
    anglecalibration.set_histogram_bin_width(0.0036); // degrees (default)
    anglecalibration.set_base_peak_ROI_width(0.18);   // degrees, half width of the base peak ROI
    const double base_peak_angle = 19.0678;           // degrees, center of the selected base peak
    anglecalibration.set_base_peak_angle(base_peak_angle);

    // read the calibration scan (frames ordered by detector angle)
    EpicsMythenFileReader file_reader("calibration_scan_0000.h5");
    const size_t num_frames = 1501;
    auto [headers, frames] = file_reader.read_n(num_frames);

    NDArray<double, 1> detector_angles(
        std::array<ssize_t, 1>{static_cast<ssize_t>(num_frames)});
    std::transform(headers.begin(), headers.end(), detector_angles.begin(),
                   [](const EpicsMythenHeader &header) { return header.detector_angle; });

    // calibrate a single module (useful to check the settings, with plots of the fit)
    anglecalibration.calibrate<true>(frames.view(), detector_angles.view(),
                                     base_peak_angle, 1);

    // calibrate all modules
    anglecalibration.calibrate(frames.view(), detector_angles.view(), base_peak_angle);

    // calibrated DG parameters, one row per module: [center, conversion, offset]
    const DGParameters &calibrated_parameters = anglecalibration.get_DGparameters();


The calibrated parameters can be loaded into a new ``MythenDetector`` later and used for conversion, see :doc:`AngularConversion`.

API:
____

.. doxygenclass:: angcal::AngleCalibration
   :members:
   :undoc-members:
   :private-members:
