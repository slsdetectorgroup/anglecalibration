import matplotlib.pyplot as plt
import h5py    
import numpy as np  
import argparse

import os

parser = argparse.ArgumentParser()
parser.add_argument("--filename")


args = parser.parse_args()

data_path = os.environ["ANGCAL_TEST_DATA"]

with h5py.File(data_path + args.filename, "r+") as f:
    photon_counts = f["entry/data/data"][()]



plt.plot(photon_counts[1286:2553])

plt.show() 



