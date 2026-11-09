"""Calibration steps as plain dictionaries: sent over ZMQ or loaded from a step file.

A step dictionary has the keys ``module_index``, ``phase`` ('coupled' or 'offset'),
``iteration``, ``similarity``, ``parameters`` (list of 3 floats), ``bin_angles``,
``detector_angles``, ``curve_modules`` and ``curves`` (numpy arrays).
"""

import json

import numpy as np

DEFAULT_ENDPOINT = "tcp://127.0.0.1:5556"

PHASES = ("coupled", "offset")

_PHASE_NAMES = {"CoupledParameters": "coupled", "Offset": "offset"}


def step_to_dict(step):
    """Convert an angcal.CalibrationStep to a step dictionary."""
    return {
        "module_index": int(step.module_index),
        "phase": _PHASE_NAMES[step.phase.name],
        "iteration": int(step.iteration),
        "similarity": float(step.similarity),
        "parameters": [float(p) for p in step.parameters],
        "bin_angles": step.bin_angles,
        "detector_angles": step.detector_angles,
        "curve_modules": step.curve_modules,
        "curves": step.base_peak_curves,
    }


def encode_step(step):
    """Encode a step dictionary as a multipart ZMQ message."""
    curves = np.ascontiguousarray(step["curves"], dtype=np.float64)
    header = {
        "module_index": step["module_index"],
        "phase": step["phase"],
        "iteration": step["iteration"],
        "similarity": step["similarity"],
        "parameters": step["parameters"],
        "bin_angles": np.asarray(step["bin_angles"]).tolist(),
        "detector_angles": np.asarray(step["detector_angles"]).tolist(),
        "curve_modules": np.asarray(step["curve_modules"]).tolist(),
        "shape": curves.shape,
    }
    return [b"step", json.dumps(header).encode(), curves.tobytes()]


def encode_status(state, **info):
    """Encode a status message e.g. 'started', 'finished' or 'error'."""
    return [b"status", json.dumps({"state": state, **info}).encode()]


def decode_message(parts):
    """Decode a multipart ZMQ message.

    Returns
    -------
    (kind, content) : (str, dict)
        kind is 'step' or 'status', content the step or status dictionary
    """
    kind = parts[0].decode()
    content = json.loads(parts[1])
    if kind == "step":
        shape = content.pop("shape")
        content["curves"] = np.frombuffer(parts[2], dtype=np.float64).reshape(shape)
        for key in ("bin_angles", "detector_angles"):
            content[key] = np.asarray(content[key], dtype=np.float64)
        content["curve_modules"] = np.asarray(content["curve_modules"], dtype=np.uint64)
    return kind, content


class StepSender:
    """Sends calibration steps to a viewer over ZMQ.

    Can be passed to AngleCalibration.set_calibration_observer. Steps are dropped
    if no viewer is listening.
    """

    def __init__(self, endpoint=DEFAULT_ENDPOINT):
        import zmq

        self._zmq = zmq
        self._socket = zmq.Context.instance().socket(zmq.PUSH)
        self._socket.setsockopt(zmq.LINGER, 2000)  # ms to wait for a viewer on close
        self._socket.connect(endpoint)
        self.num_dropped = 0

    def __call__(self, step):
        self._send(encode_step(step_to_dict(step)))

    def send_status(self, state, **info):
        self._send(encode_status(state, **info))

    def _send(self, parts):
        try:
            self._socket.send_multipart(parts, flags=self._zmq.NOBLOCK)
        except self._zmq.Again:
            self.num_dropped += 1

    def close(self):
        self._socket.close()


def load_step_file(filename):
    """Load a step file written by angcal.CalibrationRecorder.

    Returns
    -------
    list of step dictionaries
    """
    import h5py

    steps = []
    with h5py.File(filename, "r") as file:
        for module_name in sorted(file):
            module_index = int(module_name.split("_")[1])
            for phase in PHASES:
                if phase not in file[module_name]:
                    continue
                group = file[module_name][phase]
                curves = group["curves"][:]
                parameters = group["parameters"][:]
                similarity = group["similarity"][:]
                bin_angles = group["bin_angles"][:]
                detector_angles = group["detector_angles"][:]
                curve_modules = group["curve_modules"][:]
                for iteration in range(curves.shape[0]):
                    steps.append({
                        "module_index": module_index,
                        "phase": phase,
                        "iteration": iteration,
                        "similarity": float(similarity[iteration]),
                        "parameters": parameters[iteration].tolist(),
                        "bin_angles": bin_angles,
                        "detector_angles": detector_angles,
                        "curve_modules": curve_modules,
                        "curves": curves[iteration],
                    })
    return steps
