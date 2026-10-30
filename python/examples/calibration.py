
from angcal import MythenDetector, FlatField, AngleCalibration, AngularConversion, EpicsMythenFileReader, MythenFrame, PlotHelper

from angcal.utils import select_base_peak

from pathlib import Path
import os
import numpy as np 
import matplotlib.pyplot as plt

def env_data_path():
    env_value = os.environ.get("ANGCAL_TEST_DATA")
    if not env_value:
        raise RuntimeError("Environment variable ANGCAL_TEST_DATA is not set or is empty")
    
    print("Environment variable ANGCAL_TEST_DATA: ", env_value)

    return Path(env_value)

my_path = env_data_path()

dgparameters = np.loadtxt(my_path / "angcal_Jul2025_P12_0p0105.csv", delimiter=",", dtype=np.float64, skiprows=1) # initial DG parameters for calibration

# setup mythen detector specifications - stores all relevant parameters of the detector setup
mythen_detector = MythenDetector(dgparameters) 

mythen_detector.elastic_correction_factor =  0.0 
mythen_detector.detector_vertical_axis_offset = 0.0
mythen_detector.offset = 0.0 # additional offset to sample detector offset 

# setup flatfield 
flatfield = FlatField(mythen_detector)

flatfield.normalized_flatfield = np.loadtxt(my_path / "Flatfield_EkeV22p0_T11000eV_up_TESTFF1_clean_Jun2025_open_WS.raw", dtype=np.double, usecols=[1,2])


# setup angle calibration - has everything to do conversion and calibration 
anglecalibration = AngleCalibration(mythen_detector, flatfield)

anglecalibration.read_bad_channels_from_file(str(my_path / "bc2025_001_RING.chans"))

# setup mythen data file reader to read EPICS mythen hdf5 files
mythenfilereader = EpicsMythenFileReader(str(my_path / "ang1up_22keV_MIX_0p5mm_48M_a_0000.h5"))


# click through frames to select base peak for calibration
# select_base_peak(anglecalibration, mythenfilereader, module_index=0)

anglecalibration.histogram_bin_width = 0.0036  # in degrees (default value) 

anglecalibration.base_peak_ROI_width = 0.18 # in degrees width of base peak ROI 

anglecalibration.base_peak_angle = 19.0678

plotter = PlotHelper(mythen_detector, anglecalibration.angular_range, anglecalibration.histogram_bin_width)

test_header, test_frame = mythenfilereader.read_frame(209)

angularconversion = AngularConversion(mythen_detector, anglecalibration.histogram_bin_width, anglecalibration.angular_range)

angularconversion.flatfield = flatfield.normalized_flatfield
angularconversion.bad_channels = anglecalibration.bad_channels

first_header, first_frame = mythenfilereader.read_frame(0)
scale_factor = first_header.incident_intensity

angularconversion.scale_factor = scale_factor #set scale factor to get reasonable scales 

# plot one module for fixed angle width bins
module_redistributed_to_fixed_angle_bins = angularconversion.convert(test_frame.photon_counts, test_header.detector_angle, test_header.incident_intensity, test_header.exposure_time)

plotter.plot_diffraction_pattern(module_redistributed_to_fixed_angle_bins[:,0].view(), 1, test_header.detector_angle)

# plot base peak for one module
plotter.plot_base_peak(module_redistributed_to_fixed_angle_bins[:,0].view(), anglecalibration.base_peak_angle, anglecalibration.base_peak_ROI_width, 1)

num_files = 1501

mythenfilereader.seek(0) # go to first frame

headers, frames = mythenfilereader.read_n(num_files)

anglecalibration.calibrate(frames, headers["detector_angle"], anglecalibration.base_peak_angle, module_index = 1, plot_calibration_process=True) 
#anglecalibration.calibrate(frames, headers["detector_angle"], anglecalibration.base_peak_angle)

#calibrated_parameters = anglecalibration.DGparameters

#np.savetxt(my_path / "angcal_Jul2025_P12_0p0105_calibrated.csv", calibrated_parameters, delimiter=",", header="center, conversion, offset", comments="", fmt="%.6f")


print("calibration is done")





