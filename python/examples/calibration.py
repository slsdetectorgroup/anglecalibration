
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

# setup mythen detector specifications - stores all relevant parameters of the detector setup
mythen_detector = MythenDetector() 

mythen_detector.elastic_correction_factor =  0.0 
mythen_detector.detector_vertical_axis_offset = 0.0
mythen_detector.offset = 0.0 # additional offset to sample detector offset 

mythen_detector.DGparameters = np.loadtxt(my_path / "angcal_Jul2025_P12_0p0105.csv", delimiter=",", dtype=np.float64, skiprows=1) # initial DG parameters for calibration

# setup flatfield 
flatfield = FlatField(mythen_detector)

flatfield.normalized_flatfield = np.loadtxt(my_path / "Flatfield_EkeV22p0_T11000eV_up_TESTFF1_clean_Jun2025_open_WS.raw", dtype=np.double, usecols=[1,2])

# setup mythen data file reader to read EPICS mythen hdf5 files
mythenfilereader = EpicsMythenFileReader()

# setup angle calibration - has everything to do conversion and calibration 
anglecalibration = AngleCalibration(mythen_detector, flatfield, mythenfilereader)

anglecalibration.read_bad_channels_from_file(str(my_path / "bc2025_001_RING.chans"))

file_prefix = "ang1up_22keV_MIX_0p5mm_48M_a_"
num_files = 1501

detector_angles = np.empty(shape = (num_files), dtype = np.double)

frames = np.empty(shape=(num_files, mythen_detector.num_strips), dtype=np.uint32)

for frame in range(num_files):
    file_name = str(my_path / file_prefix)+f"{frame:04d}.h5"
    detector_angles[frame] = mythenfilereader.read_detector_angle(file_name)
    mythenfilereader.read_into(frames[frame, :], file_name)


# click through frames to select base peak for calibration
# select_base_peak(anglecalibration, mythenfilereader, file_list, module_index=0)

anglecalibration.histogram_bin_width = 0.0036  # in degrees (default value) 

anglecalibration.base_peak_ROI_width = 0.18 # in degrees width of base peak ROI 

anglecalibration.base_peak_angle = 19.0678

plotter = PlotHelper(mythen_detector, anglecalibration.angular_range, anglecalibration.histogram_bin_width)

test_frame_name = str(my_path / (file_prefix + "0209.h5"))

test_frame = mythenfilereader.read_frame(test_frame_name)

test_frame_angle = mythenfilereader.read_detector_angle(test_frame_name)

angularconversion = AngularConversion(mythen_detector, anglecalibration.histogram_bin_width, anglecalibration.angular_range)

angularconversion.flatfield = flatfield.normalized_flatfield
angularconversion.bad_channels = anglecalibration.bad_channels

first_frame = mythenfilereader.read_frame(str(my_path / file_prefix)+f"{0:04d}.h5")
scale_factor = first_frame.incident_intensity

angularconversion.scale_factor = scale_factor #set scale factor to get reasonable scales 

# plot one module for fixed angle width bins
module_redistributed_to_fixed_angle_bins = angularconversion.convert(test_frame.photon_counts, test_frame_angle, test_frame.incident_intensity, test_frame.exposure_time)

plotter.plot_diffraction_pattern(module_redistributed_to_fixed_angle_bins[:,0].view(), 1, test_frame_angle)

# plot base peak for one module
plotter.plot_base_peak(module_redistributed_to_fixed_angle_bins[:,0].view(), anglecalibration.base_peak_angle, anglecalibration.base_peak_ROI_width, 1)

anglecalibration.calibrate(frames, detector_angles, anglecalibration.base_peak_angle, module_index = 1, plot_calibration_process=True) 
#anglecalibration.calibrate(frames, detector_angles, anglecalibration.base_peak_angle)

#calibrated_parameters = anglecalibration.DGparameters

#np.savetxt(my_path / "angcal_Jul2025_P12_0p0105_calibrated.csv", calibrated_parameters, delimiter=",", header="center, conversion, offset", comments="", fmt="%.6f")


print("calibration is done")





