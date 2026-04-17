from angcal import FlatField, MythenDetectorSpecifications, EpicsMythenFileReader
from angcal._angcal import MythenFileReader
import numpy as np
from pathlib import Path
import os
import matplotlib.pyplot as plt

def env_data_path():
    env_value = os.environ.get("ANGCAL_TEST_DATA")
    if not env_value:
        raise RuntimeError("Environment variable ANGCAL_TEST_DATA is not set or is empty")

    return Path(env_value)

mythen_detector = MythenDetectorSpecifications()

flatfield = FlatField(mythen_detector)

flatfield.read_module_parameters_from_file(env_data_path() / "CALIB_2025/angcal_Jul2025_P12_0p0105.off")

flatfield.read_bad_channels_from_file(str(env_data_path() / "CALIB_2025/bc2025_001_RING.chans"))

bad_channels_array = flatfield.bad_channels

file_reader = EpicsMythenFileReader() 

file_list = [env_data_path() / f"H5DATA/opff_up_STRAIGHT_E22p0_T11000_{i:04d}.h5" for i in range(0,1601)]

incident_intensity_first_frame = file_reader.read_frame(str(file_list[0])).incident_intensity

flatfield.scale_factor = 78000607.0 #incident_intensity_first_frame

print("incident intensity first frame: ", incident_intensity_first_frame)

flatfield.soft_window = (3.0, 30.0)

flatfield.create_normalized_flatfield_from_filelist(file_list, file_reader)

flatfield_array = flatfield.normalized_flatfield

plt.plot(flatfield_array[:,0], label="my flatfield")
plt.legend() 
plt.show() 