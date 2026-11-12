# AngleCalibration
C++ library to calibrate PSI Mythen detector. 

## Documentation
More details about the method, installation and the C++ and Python API can be found in [documentation](https://slsdetectorgroup.github.io/anglecalibration/). 

## Build and install

### Developement Install

**Prerequisites:**
- cmake >= 3.14
- C++17 compiler (gcc >= 8)
- python >= 3.10
- HDF5

```bash
git clone https://github.com/slsdetectorgroup/anglecalibration.git
cd angle_calibration 
mkdir build
cd build

#configure using cmake
cmake ../

#build (replace 4 with the number of threads you want to use)
make -j4 
```

#### Build python binaries 

To build the python module build with option ``-DANGCAL_PYTHON_BINDINGS=ON`` e.g. 

```bash 
cmake ../ -DANGCAL_PYTHON_BINDINGS=ON
``` 

**Note:** 

Append the location of your module to your PYTHONPATH such that python can find it during import: 

```bash 
export PYTHONPATH=path_to_build_folder/build:$PYTHONPATH 
```

### Install to a custom location and use in your project

```bash
#build and install angle_calibration 
git clone https://github.com/slsdetectorgroup/anglecalibration.git
cd angle_calibration
mkdir build
cd build

#configure using cmake
cmake ../ -DCMAKE_INSTALL_PREFIX=/where/to/put/angle_calibration

#build (replace 4 with the number of threads you want to use)
make -j4 

#install
make install

#Now configure your project
 cmake .. -DCMAKE_PREFIX_PATH=/where/to/put/angcal
```

### Build Wheel 

If you only want to use the python extension you can use ``python build``. 

We recommend using a specific python environment. 

**Prerequisites:** 

- cmake >= 3.14
- C++17 compiler (gcc >= 8)
- python >= 3.10
- HDF5
- build (https://pypi.org/project/build/)


```bash 
git clone https://github.com/slsdetectorgroup/anglecalibration.git
cd angle_calibration 

#build wheel
python -m build 

pip install dist/angcal-{version}-cp{python-version}-cp{python-version}-linux_x86_64.whl
```
Import in your python project 

```python 
import angcal
```

## Python Example Usage: 

### Calibration

A toy example of a calibration in python can be found in https://github.com/slsdetectorgroup/anglecalibration/blob/main/python/examples/calibration.py. 

The example data can be found here https://zenodo.org/records/20645666 under mythenanglecalibrationtestdata.zip. It contains: 

**Flatfield_EkeV22p0_T11000eV_up_TESTFF1_clean_Jun2025_open_WS.npy**
Contains the already calculated normalized flatfield stored as a 2d numpy array of type ``np.float64``. 
The first column stores the normalized flatfield value for that channel and the second dimension the standard error of the mean (SEM). A value of -1.0 denotes that the channel was not covered by the laser beam at all or not within the soft window or denoted as a bad channel. 

**bc2025_001_RING.npy**
Contains the bad channels stored as a 1d boolean numpy array. 

**angcal_Jul2025_P12_0p0105.csv**
Contains the initial module parameters for calibration. They are saved as "Detector Group (DG) parameters" and thus include the center, conversion and offset of each module. 

The dataset also contains all the acquisition files for different detector positions stored as hdf5 files.

One can set an environment variable to the local path of ``mythenanglecalibrationtestdata``: 

```bash
export ANGCAL_CALIBRATION_DATA=/path/to/mythenanglecalibrationtestdata
```

### Conversion 

An example conversion can be found in https://github.com/slsdetectorgroup/anglecalibration/blob/main/python/examples/conversion.py. 

The example data can be found here https://zenodo.org/records/20645666 under mythenangleconversiontestdata.zip. It contains: 

**Flatfield_E17p5keV_T12500eV_up_AUGCAL2_Sep2023_open_WS_C_X_X.npy**
Contains the already calculated normalized flatfield stored as a 2d numpy array of type ``np.float64``. 
The first column stores the normalized flatfield value for that channel and the second dimension the standard error of the mean (SEM). A value of -1.0 denotes that the channel was not covered by the laser beam at all or not within the soft window or denoted as a bad channel. 

**bad_channels.npy**
Contains the bad channels stored as a 1d boolean numpy array. 

**Angcal_2E_Feb2023_P29.csv**
Contains the already calibrated module parameters. They are saved as "Detector Group (DG) parameters" and thus include the center, conversion and offset of each module. 

The dataset also contains all the acquisition files for different detector positions stored as hdf5 files.

One can set an environment variable to the local path of ``mythenangleconversiontestdata``: 

```bash
export ANGCAL_CONVERSION_DATA=/path/to/mythenangleconversiontestdata
```

### Calculating FlatField

An example of calculating the flatfield can be found in https://github.com/slsdetectorgroup/anglecalibration/blob/main/python/examples/flatfieldcalibration.py. 

The respective example data can be found here https://zenodo.org/records/20645666 under MythenFlatFieldCalibrationData.zip.
It contains: 

**bc2025_001_RING.npy**
Contains the bad channels stored as a 1d boolean numpy array. 

**angcal_Jul2025_P12_0p0105.csv**
Contains the initial module parameters for calibration. They are saved as "Detector Group (DG) parameters" and thus include the center, conversion and offset of each module. 

Additionally it contains all the flatfield acquisitions taken at different detector positions stored as hdf5 files. 

One can set an environment variable to the local path of ``MythenFlatFieldCalibrationData``: 

```bash
export ANGCAL_FLATFIELD_DATA=/path/to/MythenFlatFieldCalibrationData
```

## C++ Example Usage

### Calibration: 

A toy example of a calibration can be found in https://github.com/slsdetectorgroup/anglecalibration/blob/main/examples/example_calibration.cpp 

Information about the example data is given here [here](#calibration)

Run the example as follows: 

```bash 
export ANGCAL_TEST_DATA=/path/to/mythenanglecalibrationtestdata
#in build folder
./run_calibration
```

### Conversion: 

A toy example for a conversion can be found in https://github.com/slsdetectorgroup/anglecalibration/blob/main/examples/example_conversion.cpp 

Information about the example data is given [here](#conversion)

Run the example as follows:

```bash 
export ANGCAL_TEST_DATA=/path/to/mythenangleconversiontestdata
#in build folder
./run_conversion 
```









