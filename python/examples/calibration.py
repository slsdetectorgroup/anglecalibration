
from angcal import MythenDetectorSpecifications, FlatField, AngleCalibration, EpicsMythenFileReader, MythenFrame

from pathlib import Path
import os
import numpy as np 
import matplotlib.pyplot as plt

def env_data_path():
    env_value = os.environ.get("ANGCAL_TEST_DATA")
    if not env_value:
        raise RuntimeError("Environment variable ANGCAL_TEST_DATA is not set or is empty")

    return Path(env_value)



def plot(array : np.array, x=None): 
    if x is not None:
        plt.plot(x, array)
    else:
        plt.plot(np.arange(0, array.size,1), array)
    plt.show()

def plot_excluding_channels(array: np.array, channels_to_exclude: np.array, x = None):
    good_channels = np.logical_not(channels_to_exclude)
    if x is not None:
        plt.plot(x[good_channels], array[good_channels])
    else:
        plt.plot(np.arange(0, array.size,1)[good_channels], array[good_channels])
    plt.show()


# setup mythen detector specifications - stores all relevant parameters of the detector setup
mythendetectorspecifications = MythenDetectorSpecifications() 

mythendetectorspecifications.elastic_correction_factor =  0.0 
mythendetectorspecifications.detector_vertical_axis_offset = 0.0
mythendetectorspecifications.offset = 0.0 # additional offset to sample detector offset 

# setup flatfield 
flatfield = FlatField(mythendetectorspecifications)

flatfield.normalized_flatfield = np.loadtxt(env_data_path() / "Flatfield_EkeV22p0_T11000eV_up_TESTFF1_clean_Jun2025_open_WS.raw", dtype=np.double, usecols=[1,2])

# setup mythen data file reader to read EPICS mythen hdf5 files
mythenfilereader = EpicsMythenFileReader()

# setup angle calibration - has everything to do conversion and calibration 
anglecalibration = AngleCalibration(mythendetectorspecifications, flatfield, mythenfilereader)

anglecalibration.read_initial_calibration_from_file(str(env_data_path() / "angcal_Jul2025_P12_0p0105.off"))

anglecalibration.read_bad_channels_from_file(str(env_data_path() / "bc2025_001_RING.chans"))

file_prefix = "ang1up_22keV_MIX_0p5mm_48M_a_"
num_files = 1501
file_list = [str(env_data_path() / file_prefix)+f"{i:04d}.h5" for i in range(0,num_files)] 

#set scale factor to get reasonable scales 
first_frame = mythenfilereader.read_frame(file_list[0])
scale_factor = first_frame.incident_intensity

anglecalibration.scale_factor = scale_factor

print("scale Factor: " ,anglecalibration.scale_factor)

anglecalibration.histogram_bin_width = 0.0036  # in degrees (dfeault value) 

anglecalibration.base_peak_ROI_width = 0.18 # in degrees width of base peak ROI 

base_peak_angle = 19.0678

anglecalibration.calibrate(file_list, base_peak_angle, output_filename=str(env_data_path() / "angcal_Jul2025_P12_0p0105_new_calibrated.off"))

print("calibration is done")





