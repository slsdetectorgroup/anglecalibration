
from angcal import MythenDetectorSpecifications, FlatField, AngleCalibration, EpicsMythenFileReader, MythenFrame
from angcal import PlotHelper

from pathlib import Path
import os
import numpy as np 
import matplotlib.pyplot as plt


# Path("/home/mazzol_a/Documents/VariaMay2025/Antonio20250512/AngularConversionTestData/")

def env_data_path():
    env_value = os.environ.get("ANGCAL_TEST_DATA")
    if not env_value:
        raise RuntimeError("Environment variable ANGCAL_TEST_DATA is not set or is empty")

    print("Environment variable ANGCAL_TEST_DATA: ", env_value)
    return Path(env_value)


my_path = env_data_path()

# setup mythen detector specifications - stores all relevant parameters of the detector setup
mythendetectorspecifications = MythenDetectorSpecifications() 

mythendetectorspecifications.elastic_correction_factor =  0.0045 
mythendetectorspecifications.detector_vertical_axis_offset = 24.8 
mythendetectorspecifications.offset = 0.0 # additional offset to sample detector offset 

# setup flatfield 
flatfield = FlatField(mythendetectorspecifications)

flatfield.normalized_flatfield = np.loadtxt(my_path / "Flatfield_E17p5keV_T12500eV_up_AUGCAL2_Sep2023_open_WS_C_X_X.raw", dtype=np.double, usecols=[1,2])

# setup mythen data file reader to read EPICS mythen hdf5 files
mythenfilereader = EpicsMythenFileReader()

# setup angle calibration - has everything to do conversion and calibration 
anglecalibration = AngleCalibration(mythendetectorspecifications, flatfield, mythenfilereader)

anglecalibration.read_initial_calibration_from_file(str(my_path / "Angcal_2E_Feb2023_P29.off"))

anglecalibration.read_bad_channels_from_file(str(my_path / "bc2023_003_RING.chans"))

plt.plot(anglecalibration.bad_channels) 

#set scale factor to get reasonable scales 
frame = mythenfilereader.read_frame(str(my_path / "Fructose_0p2_60_0060.h5"))

incident_intensity_per_second = frame.incident_intensity / frame.exposure_time

scale_factor = 1280.0*round(incident_intensity_per_second/1280.0) # whatever is the same chosen in Antonios code 

anglecalibration.scale_factor = scale_factor

print("scale Factor: " ,anglecalibration.scale_factor)

file_list = [str(my_path / f"Fructose_0p2_60_006{i}.h5") for i in range(0,4)] 

anglecalibration.histogram_bin_width = 0.0036  # in degrees (dfeault value) 

anglecalibration.angular_range = (0.0, 90.5)  # in degrees 

print("Angle Range: ", anglecalibration.angular_range)

redistributed_photon_counts = anglecalibration.convert(file_list)

print("conversion is done")

#plot the converted diffraction pattern
plotter = PlotHelper(anglecalibration)

plotter.plot_diffraction_pattern(redistributed_photon_counts[:,0].view())



