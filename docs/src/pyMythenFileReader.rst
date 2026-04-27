MythenFileReader
================

File reader to read mythen files. Supported are two file formats. ``RawMythenFileReader`` reads raw files as written by the `slsDetectorPackage <https://github.com/slsdetectorgroup/slsDetectorPackage>`_. The class ``EpicsMythenFileReader`` reads hdf5 files written by ``Epics``. 
The ``EpicsMythenFileReader`` inherits from a wrapper class based on the C++ API of the HDF5 C Library. The provided dataset need to have fields `data` storing the photon counts, 
`DetectorAngle` storing the detector position in degrees and `CounterMask` storing which counters are active. The field `CounterMask` stores an integer representation of a 3-digit bit string, where each bit is set to 1 if the counter is active. 

Example
^^^^^^^

.. code-block:: text 

    CounterMask: 4 // binary representation 1 0 0 
    only channel 0 is active
    CounterMask: 3 // binary representation 0 1 1
    channel 1 and 2 are active

.. py::currentmodule:: angcal  
    
.. autoclass:: angcal._angcal.MythenFileReader
   :members:
   :undoc-members:
   :private-members:

.. autoclass:: angcal._angcal.RawMythenFileReader
   :members:
   :undoc-members:
   :private-members:

.. autoclass:: angcal._angcal.EpicsMythenFileReader
   :members:
   :undoc-members:
   :private-members:


MythenFrame
===========

.. autoclass:: angcal._angcal.MythenFrame
   :members:
   :undoc-members:
   :private-members: