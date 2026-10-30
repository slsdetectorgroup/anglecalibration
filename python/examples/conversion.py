
from angcal import MythenDetector, FlatField, AngularConversion, EpicsMythenFileReader, MythenFrame
from angcal import PlotHelper

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

dgparameters = np.loadtxt(my_path / "Angcal_2E_Feb2023_P29.csv", delimiter=",", dtype=np.float64, skiprows=1) # initial DG parameters for calibration

# setup mythen detector specifications - stores all relevant parameters of the detector setup
mythen_detector = MythenDetector(dgparameters) 

mythen_detector.elastic_correction_factor =  0.0045 
mythen_detector.detector_vertical_axis_offset = 24.8 
mythen_detector.offset = 0.0 # additional offset to sample detector offset 

print("mythen detector setup")

angular_range = (0.0, 90.5)  # in degrees
bin_width = 0.0036  # in degrees (dfeault value)
angularconversion = AngularConversion(mythen_detector, bin_width, angular_range)

angularconversion.flatfield = np.loadtxt(my_path / "Flatfield_E17p5keV_T12500eV_up_AUGCAL2_Sep2023_open_WS_C_X_X.raw", dtype=np.double, usecols=[1,2])

angularconversion.bad_channels = np.load(str(my_path / "bad_channels.npy"))

# setup mythen data file reader to read EPICS mythen hdf5 files
mythenfilereader = EpicsMythenFileReader()

print("mythen file reader setup")

#set scale factor to get reasonable scales 
frame = mythenfilereader.read_frame(str(my_path / "Fructose_0p2_60_0060.h5"))

incident_intensity_per_second = frame.incident_intensity / frame.exposure_time

scale_factor = 1280.0*round(incident_intensity_per_second/1280.0) # whatever is the same chosen in Antonios code 

print("scale factor: ", scale_factor)

angularconversion.scale_factor = scale_factor

print("scale Factor: " ,angularconversion.scale_factor)

file_list = [str(my_path / f"Fructose_0p2_60_006{i}.h5") for i in range(0,4)] 

num_files = len(file_list)

detector_angles = np.empty(shape = (num_files), dtype = np.double)
incident_intensities = np.empty(shape = (num_files), dtype = np.uint64)

frames = np.empty(shape=(num_files, mythen_detector.num_strips), dtype=np.uint32)

exposure_time = frame.exposure_time

for frame in range(num_files):
    file_name = str(my_path / f"Fructose_0p2_60_006{frame}.h5")
    detector_angles[frame] = mythenfilereader.read_detector_angle(file_name)
    incident_intensities[frame] = mythenfilereader.read_incident_intensity(file_name)
    mythenfilereader.read_into(frames[frame, :], file_name)

print("incident intensities: ", incident_intensities)

redistributed_photon_counts = angularconversion.convert(frames, detector_angles, incident_intensities, exposure_time)

print("conversion is done")

#plot the converted diffraction pattern
plotter = PlotHelper(mythen_detector, angularconversion.angular_range, angularconversion.histogram_bin_width)

plotter.plot_diffraction_pattern(redistributed_photon_counts[:,0].view())




