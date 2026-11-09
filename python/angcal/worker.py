"""Runs a calibration in its own process and sends each step to a viewer.

Usage: python -m angcal.worker --data-path /path/to/mythenanglecalibrationtestdata

The default file names are the ones of the example data (see README).
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np

from . import AngleCalibration, CalibrationRecorder, EpicsMythenFileReader, FlatField, MythenDetector
from .steps import DEFAULT_ENDPOINT, StepSender


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-path", default=os.environ.get("ANGCAL_TEST_DATA"), help="folder with the acquisition files (default: $ANGCAL_TEST_DATA)")
    parser.add_argument("--flatfield", default="Flatfield_EkeV22p0_T11000eV_up_TESTFF1_clean_Jun2025_open_WS.raw", help="normalized flatfield file")
    parser.add_argument("--initial-parameters", default="angcal_Jul2025_P12_0p0105.off", help="file with initial DG parameters (.csv or .off)")
    parser.add_argument("--bad-channels", default="bc2025_001_RING.chans", help="bad channels file")
    parser.add_argument("--file-prefix", default="ang1up_22keV_MIX_0p5mm_48M_a_", help="prefix of the acquisition files")
    parser.add_argument("--num-files", type=int, default=1501, help="number of acquisition files")
    parser.add_argument("--base-peak-angle", type=float, default=19.0678, help="angle of base peak center [degree]")
    parser.add_argument("--base-peak-roi-width", type=float, default=0.18, help="width of base peak ROI [degree]")
    parser.add_argument("--histogram-bin-width", type=float, default=0.0036, help="bin width of fixed angle width histogram [degree]")
    parser.add_argument("--module", type=int, default=None, help="only calibrate this module (default: all modules)")
    parser.add_argument("--endpoint", default=DEFAULT_ENDPOINT, help="ZMQ endpoint of the viewer")
    parser.add_argument("--record", default=None, help="hdf5 file to record the calibration steps to")
    parser.add_argument("--output", default=None, help="csv file to write the calibrated DG parameters to (only when calibrating all modules)")
    arguments = parser.parse_args(argv)
    if arguments.data_path is None:
        parser.error("--data-path is required if ANGCAL_TEST_DATA is not set")
    return arguments


def load_dg_parameters(filename):
    """Read DG parameters (center, conversion, offset) of each module from a .csv or .off file."""
    filename = Path(filename)
    if filename.suffix == ".csv":
        return np.loadtxt(filename, delimiter=",", dtype=np.float64, skiprows=1)
    # .off: module <index> center <value> +- <error> conversion <value> +- <error> offset <value> +- <error>
    rows = np.loadtxt(filename, usecols=[1, 3, 7, 11], dtype=np.float64, ndmin=2)
    module_indices = rows[:, 0].astype(int)
    parameters = np.zeros((module_indices.max() + 1, 3), dtype=np.float64)
    parameters[module_indices] = rows[:, 1:]
    return parameters


def run(arguments, sender):
    data_path = Path(arguments.data_path)

    mythen_detector = MythenDetector(load_dg_parameters(data_path / arguments.initial_parameters))

    flatfield = FlatField(mythen_detector)
    flatfield.normalized_flatfield = np.loadtxt(data_path / arguments.flatfield, dtype=np.double, usecols=[1, 2])

    anglecalibration = AngleCalibration(mythen_detector, flatfield)
    anglecalibration.read_bad_channels_from_file(str(data_path / arguments.bad_channels))

    anglecalibration.histogram_bin_width = arguments.histogram_bin_width
    anglecalibration.base_peak_ROI_width = arguments.base_peak_roi_width
    anglecalibration.base_peak_angle = arguments.base_peak_angle

    # the reader takes the first file of the series, frames need to be ordered by detector angle
    mythenfilereader = EpicsMythenFileReader(str(data_path / f"{arguments.file_prefix}{0:04d}.h5"))
    headers, frames = mythenfilereader.read_n(arguments.num_files)
    detector_angles = headers["detector_angle"]

    recorder = CalibrationRecorder(arguments.record) if arguments.record else None

    def observer(step):
        if recorder is not None:
            recorder(step)
        sender(step)

    anglecalibration.set_calibration_observer(observer)

    try:
        if arguments.module is None:
            anglecalibration.calibrate(frames, detector_angles, arguments.base_peak_angle)
            if arguments.output:
                np.savetxt(arguments.output, anglecalibration.DGparameters, delimiter=",", header="center, conversion, offset", comments="", fmt="%.15g")
        else:
            anglecalibration.calibrate(frames, detector_angles, arguments.base_peak_angle, module_index=arguments.module)
    finally:
        if recorder is not None:
            recorder.close()


def main(argv=None):
    arguments = parse_arguments(argv)
    sender = StepSender(arguments.endpoint)
    sender.send_status("started", module=arguments.module, record=arguments.record)
    try:
        run(arguments, sender)
    except Exception as error:
        sender.send_status("error", message=str(error))
        sender.close()
        raise
    sender.send_status("finished", record=arguments.record, num_dropped=sender.num_dropped)
    sender.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
