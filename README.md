# AngleCalibration
C++ library to calibrate PSI Mythen detector 

## Build and install

Prerequisites
- cmake >= 3.14
- C++17 compiler (gcc >= 8)
- python >= 3.10
- HDF5
- (boost-iostreams) 

### Developement Install

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
./run_example

```

If you want to add oline visualizations of the calibration compile with the option `ANGCAL_PLOT`
Note this requires `boost-iostreams`. 

```bash 
cd build
cmake ../ -DANGCAL_PLOT=On

#build (replace 4 with the number of threads you want to use)
make -j4
```



