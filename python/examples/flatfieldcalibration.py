from angcal import FlatField, MythenDetector, EpicsMythenFileReader
import numpy as np
from pathlib import Path
import os
import matplotlib.pyplot as plt

def env_data_path():
    env_value = os.environ.get("ANGCAL_FLATFIELD_DATA")
    if not env_value:
        raise RuntimeError("Environment variable ANGCAL_FLATFIELD_DATA is not set or is empty")

    return Path(env_value)

num_modules = 48

dgparameters = np.loadtxt(env_data_path() / "angcal_Jul2025_P12_0p0105.csv", delimiter=",", dtype=np.float64, skiprows=1) # initial DG parameters for calibration

# setup mythen detector specifications - stores all relevant parameters of the detector setup
mythen_detector = MythenDetector(dgparameters) 

flatfield = FlatField(mythen_detector)

flatfield.bad_channels = np.load(str(env_data_path() / "bc2025_001_RING.npy"))

bad_channels_array = flatfield.bad_channels

mythenfilereader = EpicsMythenFileReader(str(env_data_path() / "opff_up_STRAIGHT_E22p0_T11000_0000.h5"))

incident_intensity_first_frame = mythenfilereader.read_header(0).incident_intensity

flatfield.scale_factor = 78000607.0 #incident_intensity_first_frame

print("incident intensity first frame: ", incident_intensity_first_frame)

flatfield.soft_window = (3.0, 30.0)

for file_idx in range(0, 1601): 
    header, frame = mythenfilereader.read_frame(file_idx)
    flatfield.add_flatfield(frame.photon_counts, header.incident_intensity, header.detector_angle)

flatfield.normalize_flatfield()

flatfield_array = flatfield.normalized_flatfield

plt.plot(flatfield_array[:,0], label="my flatfield")
plt.legend() 
plt.show() 