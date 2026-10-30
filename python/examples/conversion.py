
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
mythenfilereader = EpicsMythenFileReader(str(my_path / "Fructose_0p2_60_0060.h5"))

#set scale factor to get reasonable scales 
header = mythenfilereader.read_header(0)

exposure_time = header.exposure_time

incident_intensity_per_second = header.incident_intensity / exposure_time

scale_factor = 1280.0*round(incident_intensity_per_second/1280.0) # whatever is the same chosen in Antonios code 

print("scale factor: ", scale_factor)

angularconversion.scale_factor = scale_factor

print("scale Factor: " ,angularconversion.scale_factor)

num_files = 4

mythenfilereader.seek(0)
headers, frames = mythenfilereader.read_n(num_files)

redistributed_photon_counts = angularconversion.convert(frames, headers["detector_angle"], headers["incident_intensity"], exposure_time)

print("conversion is done")

#plot the converted diffraction pattern
plotter = PlotHelper(mythen_detector, angularconversion.angular_range, angularconversion.histogram_bin_width)

plotter.plot_diffraction_pattern(redistributed_photon_counts[:,0].view())




