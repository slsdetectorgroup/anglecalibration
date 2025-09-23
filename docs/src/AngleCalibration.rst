
AngleCalibration
==================

Class to calibrate the mythen detector parameters. For the calibration the "Best computing" BC parameters are used. However, as initial parameters it expects the "historic detector group" DG parameters. 
See :ref:`parametersets` for a more detailed description of the different detector parameters. 

.. note:: 
    The calibration process will automatically be visualized if you compiled with option ``ANGCAL_PLOT`` (:ref:`installation`). 

..
    TODO: maybe move detailed description to another page. 

Example Usage: 
_______________

.. code-block:: cpp

    std::shared_ptr<MythenDetectorSpecifications> mythen_detector_ptr =
        std::make_shared<MythenDetectorSpecifications>();

    mythen_detector_ptr->read_bad_channels_from_file("bad_channels.txt");

    std::shared_ptr<FlatField> flat_field_ptr =
        std::make_shared<FlatField>(mythen_detector_ptr);

    flat_field_ptr->read_flatfield_from_file("flatfield.raw"); 

    flat_field_ptr->calculate_inverse_normalized_flatfield(); 

    AngleCalibration anglecalibration(mythen_detector_ptr, flat_field_ptr);

    anglecalibration.read_initial_calibration_from_file(
        "initial_DG_parameters.off");

    //dummy acquisitions
    std::vector<std::string> filelist{"file_01.hdf5", "file_02.hdf5", "file_03.hdf5"}; 

    double base_peak_angle = -49.75; 

    //calibrate parameters for module 10
    anglecalibration.calibrate(filelist, base_peak_angle, 10); 


.. note:: 
    The method ``read_initial_calibration_from_file`` expects a text file in the format below. Each line denotes the DG parameters for each module. 

    .. code-block:: text

        module 0 center  643.320033188550 +- 0.0000 conversion 0.657661167459868E-04 +- 0.0000 offset  0.00000000000000     +- 0.0000 
        module 1 center  633.044070300816 +- 0.0000 conversion 0.657618957538105E-04 +- 0.0000 offset  5.00486981153634     +- 0.0000


    If youre file format deviates from the above format implement your custom file reader inheriting from the class ``SimpleFileInterface`` (:ref:`simplefileinterface`) and pass it to the constructor. 
    Similarly for the acquisition files, ``AngleCalibration`` uses per default the class ``MythenFileReader`` (:ref:`mythenfilereader`) expecting hdf5 files with entries 'data', 'DetectorAngle' and 'CounterMask'. 


.. doxygenclass:: angcal::AngleCalibration
   :members:
   :undoc-members:
   :private-members: