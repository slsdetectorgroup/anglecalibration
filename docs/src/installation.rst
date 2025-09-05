.. _installation: 

Installation
============

To install the Angle Calibration you need to fetch it from GitHub and build it locally on your system. 


.. code-block:: bash 

    git clone https://github.com/slsdetectorgroup/anglecalibration.git
    cd angle_calibration


To build the library locally follow the instructions below: 

.. code-block:: bash

    mkdir build
    cd build
    cmake ../ -DCMAKE_INSTALL_PREFIX=/where/to/put/angle_calibration
    make -j9 
    make install 

If you want to run an example calibration use: 

.. code-block:: bash

    export ANGCAL_TEST_DATA=/path/to/your_data  
    #in build folder
    ./run_example

Use with online visualization: 
------------------------------

If you want to plot and visualize the calibration process build with the following option: 

.. code-block:: bash 

    #in build folder
    cmake ../ -DANGCAL_PLOT=On


Build and use in your C++ project
-----------------------------------

If you want to include it in your own cmake project do:

.. code-block:: bash

    #in build folder
    cmake ../ -DCMAKE_INSTALL_PREFIX=/where/to/put/angle_calibration
    make -j9 
    make install 

In your own cmake project include: 

.. code-block:: cmake

        find_package(angle_calibration REQUIRED)

And build with: 

.. code-block:: bash

    #in build folder
    cmake .. -DCMAKE_PREFIX_PATH=/where/to/put/angle_calibration


Build python module and use in python project: 
-----------------------------------------------

To use the library from python build the python bindings by building with the following option: 

.. code-block:: bash 

    #in build folder
    cmake ../ -DANGCAL_PYTHON_BINDINGS=On 

Append the generated python module to your ``PYTHONPATH``: 
    
.. code-block:: bash 

    export PYTHONPATH=path_to_build_folder:$PYTHONPATH 


Import the module in your python script: 

.. code-block:: python 

    import angcal 

Build the tests: 
-----------------

To build the tests build with the following option: 

.. code-block:: bash 

    #in build folder 
    cmake ../ -DANGCAL_TESTS=On 
    make -j4

Run the tests: 

.. code-block:: bash 

    #in build folder
    ./run_tests 

Build the documentation: 
-------------------------

Building the documentation requires doxygen, sphinx and breathe. To build the documentation build with the following option: 

.. code-block:: bash

    #in build folder 
    cmake ../ -DANGCAL_DOCS=On -DANGCAL_PYTHON_BINDINGS=On 
    make -j4
    make docs 

        
