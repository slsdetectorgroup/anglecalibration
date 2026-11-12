AngularConversion
=================

Utilities to perform angular conversion. As conversion parameters it expects the "historic detector group" DG parameters. 
See :ref:`parametersets` for a more detailed description of the different detector parameters. 

..
    TODO: maybe move detailed description to another page. 

.. doxygenstruct:: angcal::HistogramMetaData
   :members:
   :undoc-members:
   :private-members:


.. doxygenfunction:: angcal::convert(const MythenDetector&, const HistogramMetaData&, const NDView<uint32_t, 2>, const NDView<double, 1>, const NDView<uint64_t, 1>, const double, const NDView<const bool, 1>, const NDView<const double, 2>, const double)


.. doxygenfunction:: angcal::convert(const MythenDetector&, const HistogramMetaData&, const NDView<uint32_t, 1>, const double, const uint64_t, const double, const NDView<const bool, 1>, const NDView<const double, 2>, const double)
