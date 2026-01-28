
from angcal import MythenDetectorSpecifications, FlatField, AngleCalibration, EpicsMythenFileReader, MythenFrame

from pathlib import Path
import os
import numpy as np 
import matplotlib.pyplot as plt

def data_path():
    return Path("/home/mazzol_a/Documents/VariaMay2025/Antonio20250512/AngularConversionTestData/")

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

# setup flatfield 
flatfield = FlatField(mythendetectorspecifications)

flatfield.normalized_flatfield = np.loadtxt(data_path() / "Flatfield_E17p5keV_T12500eV_up_AUGCAL2_Sep2023_open_WS_C_X_X.raw", dtype=np.double, usecols=[1,2])

#plot(flatfield.inverse_normalized_flatfield[:,0])

# setup mythen data file reader to read EPICS mythen hdf5 files
mythenfilereader = EpicsMythenFileReader()

# setup angle calibration - has everything to do conversion and calibration 
anglecalibration = AngleCalibration(mythendetectorspecifications, flatfield, mythenfilereader)

anglecalibration.read_initial_calibration_from_file(str(data_path() / "Angcal_2E_Feb2023_P29.off"))

anglecalibration.read_bad_channels_from_file(str(data_path() / "bc2023_003_RING.chans"))

bad_channels = anglecalibration.bad_channels

#plot(bad_channels)

#set scale factor to get reasonable scales 
frame = mythenfilereader.read_frame(str(data_path() / "Fructose_0p2_60_0060.h5"))
anglecalibration.scale_factor = frame.incident_intensity/10

print("scale Factor: " ,anglecalibration.scale_factor)

file_list = [str(data_path() / f"Fructose_0p2_60_006{i}.h5") for i in range(0,4)] 

redistributed_photon_counts = anglecalibration.convert(file_list)

# plot converted data
bin_indices = np.arange(0, redistributed_photon_counts.size,1)
bin_to_diffraction_angle = lambda bin_index : bin_index * anglecalibration.histogram_bin_width + mythendetectorspecifications.min_angle 

bin_in_degrees = np.apply_along_axis(bin_to_diffraction_angle, 0, bin_indices)

# data stores anything between -180, 180 degrees, we only want to plot region of interest 
zero_channels = redistributed_photon_counts == 0

#plot_excluding_channels(merged_redistributed_photon_counts, zero_channels, bin_in_degrees)

### actual diffraction pattern for comparison 

actual_diffraction_pattern = np.loadtxt(data_path() / "Fructose_0p2_60_m_Alice_WAXS.xye", dtype=np.double)

#plot(actual_diffraction_pattern[:,1], actual_diffraction_pattern[:,0])

plt.plot(bin_in_degrees[~zero_channels],redistributed_photon_counts[~zero_channels], label="my conversion")
plt.plot(actual_diffraction_pattern[:,0], actual_diffraction_pattern[:,1], label="actual diffraction pattern")

plt.xlabel("Diffraction Angle (degrees)")
plt.ylabel("Photon Counts")
plt.legend()
plt.show()




