Normalized FlatField
=====================

The ``FlatField`` class is responsible for storing the flatfield correction factors. It provides methods to read and create normalized flatfield from a list of files containing flatfield acquisitions.

.. note:: 
    The method ``read_module_parameters_from_file`` expects a text file in the format below. Each line denotes the DG parameters for each module. 

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

.. autoclass:: angcal._angcal.FlatField
    :special-members: __init__
    :show-inheritance:
    :inherited-members:
    :member-order: groupwise
