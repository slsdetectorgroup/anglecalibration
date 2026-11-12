MythenFileReader
================

File reader to read mythen hdf5 files written by ``Epics``. 
The ``EpicsMythenFileReader`` inherits from a wrapper class based on the C++ API of the HDF5 C Library. The provided dataset need to have fields ``data`` storing the photon counts, 
``DetectorAngle`` storing the detector position in degrees and ``CounterMask`` storing which counters are active. The field ``CounterMask`` stores an integer representation of a 3-digit bit string, where each bit is set to 1 if the counter is active. 

Example
^^^^^^^

.. code-block:: text 

    CounterMask: 4 // binary representation 1 0 0 
    only channel 0 is active
    CounterMask: 3 // binary representation 0 1 1
    channel 1 and 2 are active


To create a new instance of the file reader pass the name of the first file of the dataset to read. Frames are indexed consecutively over all files of the series, when all frames of a file are read the next file in the series is opened automatically.

**Example Usage:**

.. code-block:: python 
   
    from angcal import EpicsMythenFileReader

    # create a new file reader instance
    mythen_file_reader = EpicsMythenFileReader("file_00.hdf5")

    # read the first frame of the dataset e.g. frame 0 
    header, frame = mythen_file_reader.read_frame(0)

    header, frame = mythen_file_reader.read_frame()  # reads the next frame (e.g. frame 1)


.. py:currentmodule:: angcal  
    

.. autoclass:: angcal._angcal.EpicsMythenFileReader
   :special-members: __init__
   :show-inheritance:
   :inherited-members:
   :member-order: groupwise
   :private-members:


MythenFrame
===========

.. autoclass:: angcal._angcal.MythenFrame
   :members:
   :member-order: groupwise


MythenHeader
============

.. autoclass:: angcal._angcal.EpicsMythenHeader
   :members:
   :member-order: bysource
