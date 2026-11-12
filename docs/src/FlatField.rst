Normalized FlatField
=====================

The ``FlatField`` class is responsible for storing the flatfield correction factors. It provides methods to create a normalized flatfield from a series of flatfield acquisitions,
e.g. a scan of the detector over a homogeneously illuminated region (fluorescence or amorphous scatterer).

A flatfield is created from multiple acquisitions from a flatfield scan. 
Each photon count of each scan is corrected by the incident intensity and the solid angle of the strip. The multiple scans are then accumulated and normalized such that the mean over all good strips is one. The resulting flatfield can be stored and reused later.

1. **Accumulate** – every acquisition of the flatfield scan is passed to :cpp:func:`FlatField::add_flatfield`. The photon counts of each strip are corrected and summed up.
2. **Normalize** – once all acquisitions are added, :cpp:func:`FlatField::normalize_flatfield` applies the solid angle correction (same for each scan) and normalizes the accumulated flatfield such that the mean over all good strips is one.

.. doxygenclass:: angcal::FlatField
   :members:
   :undoc-members:
   :private-members: