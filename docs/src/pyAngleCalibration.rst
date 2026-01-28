AngleCalibration
=================

Class to perform angular conversionand to calibrate the mythen detector parameters used for conversion. For the calibration the "Best computing" parameters are used. However, as initial parameters or for conversion it expects the "historic detector group" parameters. 
See :ref:`parametersets` for a more detail description of the different detector parameters. 

..
    TODO: maybe move detailed description to another page. 

Example Usage for Conversion: 
______________________________

.. code-block:: python

    from angcal import MythenDetectorSpecifications, FlatField, AngleCalibration, EpicsMythenFileReader, MythenFrame
    from pathlib import Path
    import numpy as np

    data_path = Path("/path/to/your/data/files")  # adjust this path

    # setup mythen detector specifications - stores all relevant parameters of the detector setup
    mythendetectorspecifications = MythenDetectorSpecifications(offset=0, num_counters=1) 

    # setup flatfield 
    flatfield = FlatField(mythendetectorspecifications)

    flatfield.normalized_flatfield = np.loadtxt(data_path / "Flatfield_E17p5keV_T12500eV_up_AUGCAL2_Sep2023_open_WS_C_X_X.raw", dtype=np.double, usecols=[1,2])

    # setup mythen data file reader to read EPICS mythen hdf5 files
    mythenfilereader = EpicsMythenFileReader()

    # setup angle calibration - has everything to do conversion and calibration 
    anglecalibration = AngleCalibration(mythendetectorspecifications, flatfield, mythenfilereader)

    anglecalibration.read_initial_calibration_from_file(str(data_path / "Angcal_2E_Feb2023_P29.off"))

    anglecalibration.read_bad_channels_from_file(str(data_path / "bc2023_003_RING.chans"))

    #set scale factor to get reasonable scales 
    frame = mythenfilereader.read_frame(str(data_path / "Fructose_0p2_60_0060.h5"))
    anglecalibration.scale_factor = frame.incident_intensity

    file_list = [str(data_path / f"Fructose_0p2_60_006{i}.h5") for i in range(0,4)] 

    redistributed_photon_counts = anglecalibration.convert(file_list)

API:
____

.. note:: 
    The method ``read_initial_calibration_from_file`` expects a text file in the format below. Each line denotes the DG parameters for each module. 

    .. code-block:: text

        module 0 center  643.320033188550 +- 0.0000 conversion 0.657661167459868E-04 +- 0.0000 offset  0.00000000000000     +- 0.0000 
        module 1 center  633.044070300816 +- 0.0000 conversion 0.657618957538105E-04 +- 0.0000 offset  5.00486981153634     +- 0.0000


    If youre file format deviates from the above format read your file into a ``np.ndarray(,3)`` where the first dimension stores the centers the second the conversions and the third the offsets and set ``DGparameters`` directly. 


.. note:: 
    The method ``read_bad_channels_from_file`` expects a text file using the format below. Each line denotes the index of a bad channel. Consecutive bad channels can be stored on one line e.g. channels with 20-42 are bad channels.

    .. code-block:: text

        1
        2
        20-42
        100
        1002-1004
    
    If youre file format deviates from the above format read your file into a boolean ``np.ndarray`` and set ``bad_channels`` directly. 

.. py::currentmodule:: angcal  
    
.. autoclass:: angcal._angcal.AngleCalibration
    :special-members: __init__
    :show-inheritance:
    :inherited-members:
    :member-order: groupwise




