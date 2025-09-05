
MythenDetectorSpecifications
============================

MythenDetectorSpecifications contains all detector specific parameters, such as the number of modules, internal detector parameters like the strip width. 

.. note:: 
    The method ``read_bad_channels_from_file`` expects a text file using the format below. Each line denotes the index of a bad channel. Consecutive bad channels can be stored on one line e.g. channels with 20-42 are bad channels.

    .. code-block:: text

        1
        2
        20-42
        100
        1002-1004
    
    If youre file format deviates from the above format read youre file into a boolean ``aare::NDArray`` and pass to ``set_bad_channels`` directly. Alternatively implement your custom Filereader inheriting from ``SimpleFileInterface`` (:ref:`simplefileinterface`) and pass it to the constructor. 

.. doxygenclass:: angcal::MythenDetectorSpecifications
   :members:
   :undoc-members:
   :private-members:
