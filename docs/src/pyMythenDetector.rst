MythenDetector 
============================ 

MythenDetector contains all detector specific parameters, such as internal detector parameters like the strip width, the measured dead-time of the detector, as well as the module parameters and methods to calculate the diffraction angle of each strip.  
It also contains all setup-specific parameters, such as the number of modules or the average_distance of the pixel to the sample. 

.. 
   maybe add a sketch with these parameters 

.. Note:: 
   Due to legacy reasons, MythenDetector can be constructed from either the "historic detector group" DG parameters or the "Best computing" BC parameters. The constructor converts either parameter set to the other one and stores both parameter sets as members.

.. py:currentmodule:: angcal 

.. autoclass:: angcal._angcal.MythenDetector 
    :special-members: __init__
    :show-inheritance:
    :inherited-members:
    :member-order: groupwise
    
