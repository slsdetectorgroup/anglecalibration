import math
import sys
import pytest 
from angcal import EpicsMythenFileReader, AngleCalibration, MythenDetectorSpecifications, FlatField

from conftest import test_data_path

# Constants (adjust if they come from elsewhere in your code)
one = 1.0
ten_eps_DP = 10.0 * sys.float_info.epsilon

# Example placeholders (replace with your actual values)
betaxmax=0.367879441171442239909878595207253776 # exp(-1) 
beta_ratecorr = [76.08e-9,31.8e-9, 95.e-9] #why are there three? per default use one 
beta_ratecorr = beta_ratecorr[0]  # Use the first value as default


def evalyc(xx):
    """
    Equivalent of Fortran function EVALYC
    Computes inverse Lambert W via fixed-point iteration
    """
    yy0 = xx

    while True:
        yy = xx * math.exp(yy0)     # inverse Lambert W iteration
        tt = abs(one - yy / yy0)
        if tt < ten_eps_DP:
            break
        yy0 = yy

    return yy


def evalcorr(x):
    """
    Equivalent of Fortran function EVALCORR
    """
    evalcorr_val = one

    if x <= one:
        return evalcorr_val

    xx = min(betaxmax, x * beta_ratecorr)
    evalcorr_val = evalyc(xx) / xx

    return evalcorr_val


def test_rate_correction(test_data_path):
    """ Test rate correction calculation"""

    filereader = EpicsMythenFileReader()

    anglecalibration = AngleCalibration(MythenDetectorSpecifications(), FlatField(MythenDetectorSpecifications()), filereader)

    frame = filereader.read_frame(str(test_data_path / "Fructose_0p2_60_0060.h5")) 

    anglecalibration.read_bad_channels_from_file(str(test_data_path / "bc2023_003_RING.chans"))

    photon_counts = frame.photon_counts
    photon_counts_error = photon_counts**2 # Poisson statistics

    for i in range(len(photon_counts)):
        if(anglecalibration.bad_channels[i]):
            continue

        if(photon_counts[i] == 0): 
            rate_corrected_pc = 0.0 
            rate_corrected_pc_errors = 0.0 
        else: 
            rate_corrected_pc, rate_corrected_pc_errors = anglecalibration.rate_correction(photon_counts[i], photon_counts_error[i], frame.exposure_time) 


        expected_rate_corrected_pc = photon_counts[i] * evalcorr(photon_counts[i] / frame.exposure_time)

        #expected_rate_corrected_pc_errors = photon_counts_error[i] * (evalcorr(photon
        assert(rate_corrected_pc == pytest.approx(expected_rate_corrected_pc, rel=1e-6))

    

    