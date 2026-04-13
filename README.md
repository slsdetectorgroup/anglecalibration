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

## Example Usage

A toy example of a calibration can be found in https://github.com/slsdetectorgroup/anglecalibration/blob/main/examples/example_calibration.cpp 

The example data is located in https://gitea.psi.ch/angcal/VariaMay2025. It contains: 

- Flatfield_E17p5keV_T8751eV_MIX_Mar2021_open_WS.raw (the already calculated inverse normalized flatfield)

- bcX.txt (the bad channels)

- angcal_Mar2021_P10.off (the initial module parameters to calibrate (saved as "Detector Group (DG) parameters"))

- all the acquisition files for different detector positions (one frame per file) (Note that the file only contains the photon counts for the good channels.)

Run the example as follows: 

```bash 
export ANGCAL_TEST_DATA=/path/to/local/gitea.psi.ch/angcal/VariaMay2025 
#in build folder
./run_calibration
```





