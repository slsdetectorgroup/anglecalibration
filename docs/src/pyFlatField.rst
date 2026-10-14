Normalized FlatField
=====================

The ``FlatField`` class is responsible for storing the flatfield correction factors. It provides methods to read and create normalized flatfield from a list of files containing flatfield acquisitions.

.. note:: 
    The method ``read_bad_channels_from_file`` expects a text file using the format below. Each line denotes the index of a bad channel. Consecutive bad channels can be stored on one line e.g. channels with 20-42 are bad channels.

    .. code-block:: text

        1
        2
        20-42
        100
        1002-1004
    
    If youre file format deviates from the above format read your file into a boolean ``np.ndarray`` and set ``bad_channels`` directly. 

.. py:currentmodule:: angcal  

.. autoclass:: angcal._angcal.FlatField
    :special-members: __init__
    :show-inheritance:
    :inherited-members:
    :member-order: groupwise
