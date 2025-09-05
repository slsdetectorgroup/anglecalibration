
.. _mythenfilereader: 

MythenFileReader
================

File reader to read acquisition files stored as hdf5. The filereader inherits from a wrapper class based on the C++ API of the HDF5 C Library. The provided dataset need to have fields `data` storing the photon counts, 
`DetectorAngle` storing the detector position in degrees and `CounterMask` storing which counters are active. The field `CounterMask` stores an integer representation of a 3-digit bit string, where each bit is set to 1 if the counter is active. 

Example
-------

.. code-block:: text 

    CounterMask: 4 // binary representation 1 0 0 
    only channel 0 is active
    CounterMask: 3 // binary representation 0 1 1
    channel 1 and 2 are active

.. doxygenclass:: angcal::MythenFileReader
   :members:
   :undoc-members:
   :private-members: