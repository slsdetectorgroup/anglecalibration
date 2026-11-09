"""Web GUI to run a calibration and to view its steps live or from a step file.

Usage: python -m angcal.gui --data-path /path/to/mythenanglecalibrationtestdata

The calibration runs in its own process (angcal.worker) and sends each step
over ZMQ. Requires nicegui, plotly, pyzmq and h5py.
"""

import argparse
import asyncio
import os
import re
import sys
from collections import deque

import numpy as np
from nicegui import app, background_tasks, ui

from .steps import DEFAULT_ENDPOINT, decode_message, load_step_file

# chart colors
SURFACE = "#fcfcfb"
TEXT_SECONDARY = "#52514e"
TEXT_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
AXIS = "#c3c2b7"
SERIES = ("#2a78d6", "#eb6834")  # module, previous module
SEQUENTIAL = ("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")

PHASE_LABELS = {"coupled": "Module parameters", "offset": "Offset to previous module"}
PARAMETER_LABELS = ("Module center to sample distance", "Angle module center to module normal [degree]", "Angle module center to beam [degree]")

ANSI_ESCAPE = re.compile(r"\x1b\[[0-9;]*m")


class Session:
    """Steps received from the worker or loaded from a file, shared by all browser tabs."""

    def __init__(self):
        self.steps = {}  # (module_index, phase) -> list of step dictionaries
        self.last_key = None
        self.version = 0  # incremented when steps change
        self.status = "Idle"
        self.process = None
        self.log_lines = deque(maxlen=500)
        self.num_log_lines = 0

    def clear(self):
        self.steps = {}
        self.last_key = None
        self.version += 1

    def add_step(self, step):
        key = (step["module_index"], step["phase"])
        self.steps.setdefault(key, []).append(step)
        self.last_key = key
        self.version += 1

    def num_steps(self):
        return sum(len(steps) for steps in self.steps.values())

    def log(self, line):
        self.log_lines.append(line)
        self.num_log_lines += 1

    def handle_message(self, kind, content):
        if kind == "step":
            self.add_step(content)
            self.status = f"Running: module {content['module_index']}, {PHASE_LABELS[content['phase']].lower()}"
        elif content["state"] == "started":
            self.status = "Running"
        elif content["state"] == "finished":
            self.status = f"Finished: {self.num_steps()} steps"
            if content.get("num_dropped"):
                self.status += f" ({content['num_dropped']} not received)"
        elif content["state"] == "error":
            self.status = f"Error: {content['message']}"

    def is_running(self):
        return self.process is not None and self.process.returncode is None

    async def start(self, worker_arguments):
        self.clear()
        self.status = "Starting"
        self.process = await asyncio.create_subprocess_exec(
            sys.executable, "-u", "-m", "angcal.worker", *worker_arguments,
            stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.STDOUT)
        background_tasks.create(self._watch(self.process))

    async def _watch(self, process):
        async for line in process.stdout:
            self.log(ANSI_ESCAPE.sub("", line.decode(errors="replace")).rstrip())
        returncode = await process.wait()
        if returncode != 0 and not self.status.startswith("Error"):
            self.status = "Stopped" if returncode < 0 else f"Error: calibration exited with code {returncode}"

    def stop(self):
        if self.is_running():
            self.process.terminate()

    def load(self, filename):
        steps = load_step_file(filename)
        self.clear()
        for step in steps:
            self.add_step(step)
        self.status = f"Loaded {self.num_steps()} steps from {os.path.basename(filename)}"


session = Session()


async def receive_steps(endpoint):
    import zmq
    import zmq.asyncio

    socket = zmq.asyncio.Context.instance().socket(zmq.PULL)
    socket.bind(endpoint)
    while True:
        session.handle_message(*decode_message(await socket.recv_multipart()))


def to_list(values):
    """numpy array to list with NaN replaced by None (plotted as gap)"""
    return [None if value != value else value for value in np.asarray(values, dtype=float).tolist()]


def layout(x_title, y_title, height):
    axis = {"gridcolor": GRIDLINE, "linecolor": AXIS, "tickcolor": AXIS, "ticks": "outside", "zeroline": False, "automargin": True}
    return {
        "height": height,
        "margin": {"l": 72, "r": 24, "t": 28, "b": 52},
        "paper_bgcolor": SURFACE,
        "plot_bgcolor": SURFACE,
        "font": {"family": 'system-ui, -apple-system, "Segoe UI", sans-serif', "size": 12, "color": TEXT_SECONDARY},
        "xaxis": {**axis, "title": {"text": x_title}},
        "yaxis": {**axis, "title": {"text": y_title}},
        "showlegend": False,
    }


def base_peak_marker(step):
    """vertical line at the base peak angle (center of the ROI)"""
    base_peak_angle = 0.5 * (step["bin_angles"][0] + step["bin_angles"][-1])
    return {
        "shapes": [{"type": "line", "x0": base_peak_angle, "x1": base_peak_angle, "yref": "paper", "y0": 0, "y1": 1,
                    "line": {"color": TEXT_MUTED, "width": 1}}],
        "annotations": [{"x": base_peak_angle, "yref": "paper", "y": 1, "yanchor": "bottom", "showarrow": False,
                         "text": "base peak", "font": {"color": TEXT_MUTED, "size": 11}}],
    }


def covered_curves(step):
    """indices of the curves which cover the base peak ROI"""
    return [i for i in range(step["curves"].shape[0]) if np.isfinite(step["curves"][i]).any()]


def curves_figure(step, max_counts):
    module_index = step["module_index"]
    bin_angles = step["bin_angles"].tolist()
    two_modules = step["phase"] == "offset"
    traces = []
    labeled = set()
    for i in covered_curves(step):
        curve_module = int(step["curve_modules"][i])
        traces.append({
            "type": "scatter", "mode": "lines", "x": bin_angles, "y": to_list(step["curves"][i]),
            "line": {"width": 1.5, "color": SERIES[0] if curve_module == module_index else SERIES[1]},
            "opacity": 0.6,
            "name": f"Module {curve_module}", "legendgroup": str(curve_module),
            "showlegend": two_modules and curve_module not in labeled,
            "hovertemplate": f"detector angle {step['detector_angles'][i]:.2f}°<br>%{{x:.4f}}°: %{{y:,.0f}} counts<extra>Module {curve_module}</extra>",
        })
        labeled.add(curve_module)
    figure_layout = layout("Diffraction angle [degree]", "Photon counts", 420)
    figure_layout["yaxis"]["range"] = [0, 1.05 * max_counts]
    figure_layout["xaxis"]["range"] = [bin_angles[0], bin_angles[-1]]
    figure_layout.update(base_peak_marker(step))
    if two_modules:
        figure_layout["showlegend"] = True
        figure_layout["legend"] = {"orientation": "h", "x": 0, "y": 1.0, "yanchor": "bottom", "itemclick": False, "itemdoubleclick": False}
    return {"data": traces, "layout": figure_layout}


def heatmap_figure(step, max_counts):
    rows = covered_curves(step)
    labels = [f"Module {int(step['curve_modules'][i])}, detector angle {step['detector_angles'][i]:.2f}°" for i in rows]
    trace = {
        "type": "heatmap", "x": step["bin_angles"].tolist(), "y": list(range(len(rows))),
        "z": [to_list(step["curves"][i]) for i in rows],
        "customdata": [[label] * step["curves"].shape[1] for label in labels],
        "zmin": 0, "zmax": max_counts,
        "colorscale": [[i / (len(SEQUENTIAL) - 1), color] for i, color in enumerate(SEQUENTIAL)],
        "colorbar": {"title": {"text": "Photon counts", "side": "right"}, "thickness": 12, "outlinewidth": 0},
        "hovertemplate": "%{customdata}<br>%{x:.4f}°: %{z:,.0f} counts<extra></extra>",
    }
    figure_layout = layout("Diffraction angle [degree]", "Acquisition", 420)
    figure_layout["yaxis"].update({"showgrid": False, "autorange": "reversed"})
    figure_layout["xaxis"]["showgrid"] = False
    figure_layout.update(base_peak_marker(step))
    return {"data": [trace], "layout": figure_layout}


def similarity_figure(steps, index):
    iterations = list(range(len(steps)))
    similarities = [step["similarity"] for step in steps]
    line = {"type": "scatter", "mode": "lines", "x": iterations, "y": similarities, "line": {"width": 2, "color": SERIES[0]},
            "hovertemplate": "step %{x}: %{y:,.0f}<extra></extra>"}
    current = {"type": "scatter", "mode": "markers", "x": [index], "y": [similarities[index]],
               "marker": {"size": 11, "color": SERIES[0], "line": {"width": 2, "color": SURFACE}}, "hoverinfo": "skip"}
    figure_layout = layout("Step", "Similarity criterion", 260)
    figure_layout["xaxis"].update({"dtick": max(1, len(steps) // 10), "range": [-0.5, max(len(steps) - 0.5, 1)]})
    figure_layout["yaxis"]["rangemode"] = "tozero"
    return {"data": [line, current], "layout": figure_layout}


def empty_figure(height):
    figure_layout = layout("", "", height)
    figure_layout["xaxis"]["visible"] = False
    figure_layout["yaxis"]["visible"] = False
    figure_layout["annotations"] = [{"text": "No steps yet", "xref": "paper", "yref": "paper", "x": 0.5, "y": 0.5, "showarrow": False,
                                     "font": {"color": TEXT_MUTED, "size": 14}}]
    return {"data": [], "layout": figure_layout}


def create_page(defaults):

    @ui.page("/")
    def index():
        selection = {"module": None, "phase": "coupled", "index": 0, "version": -1, "log_lines": 0, "busy": False}

        ui.query("body").style("background-color: #f9f9f7")

        with ui.header().classes("items-center justify-between bg-white text-black shadow-sm"):
            ui.label("Angle calibration").classes("text-lg font-medium")
            status = ui.label(session.status).classes("text-sm")

        with ui.row().classes("w-full no-wrap items-start gap-4"):

            with ui.column().classes("w-80 shrink-0 gap-4"):
                with ui.card().classes("w-full"):
                    ui.label("Run calibration").classes("font-medium")
                    data_path = ui.input("Data folder", value=defaults.data_path or "").classes("w-full")
                    base_peak_angle = ui.number("Base peak angle [degree]", value=19.0678, format="%.4f").classes("w-full")
                    roi_width = ui.number("Base peak ROI width [degree]", value=0.18, format="%.4f").classes("w-full")
                    module = ui.input("Module", value="all").props('hint="all, or a module index"').classes("w-full")
                    record = ui.input("Record steps to", value=defaults.record).classes("w-full")
                    with ui.row():
                        start_button = ui.button("Start")
                        stop_button = ui.button("Stop", color="grey").props("outline")
                with ui.card().classes("w-full"):
                    ui.label("Open recorded steps").classes("font-medium")
                    step_file = ui.input("Step file", value=defaults.record).classes("w-full")
                    load_button = ui.button("Open").props("outline")

            with ui.column().classes("grow min-w-0 gap-4"):
                with ui.card().classes("w-full"):
                    with ui.row().classes("w-full items-center gap-6"):
                        module_select = ui.select({}, label="Module").classes("w-32")
                        phase_toggle = ui.toggle(PHASE_LABELS, value="coupled")
                        view_toggle = ui.toggle({"curves": "Curves", "heatmap": "Heatmap"}, value="curves")
                        follow = ui.switch("Follow calibration", value=True)
                    with ui.row().classes("w-full items-center no-wrap gap-4"):
                        step_label = ui.label("Step").classes("w-28 shrink-0 text-sm")
                        slider = ui.slider(min=0, max=1, step=1, value=0).props("label snap markers")
                    curves_plot = ui.plotly(empty_figure(420)).classes("w-full")
                with ui.row().classes("w-full no-wrap items-stretch gap-4"):
                    with ui.card().classes("grow min-w-0"):
                        ui.label("Similarity of the base peaks (lower is better)").classes("text-sm font-medium")
                        similarity_plot = ui.plotly(empty_figure(260)).classes("w-full")
                    with ui.card().classes("w-96 shrink-0"):
                        ui.label("Parameters at this step").classes("text-sm font-medium")
                        parameter_labels = []
                        for label in PARAMETER_LABELS:
                            ui.label(label).classes("text-xs text-gray-600 mt-2")
                            with ui.row().classes("items-baseline gap-3"):
                                value_label = ui.label("-").classes("text-lg")
                                change_label = ui.label("").classes("text-xs text-gray-600")
                            parameter_labels.append((value_label, change_label))
                with ui.expansion("Calibration log").classes("w-full bg-white"):
                    log = ui.log(max_lines=500).classes("w-full h-48 text-xs")

        def render():
            steps = session.steps.get((selection["module"], selection["phase"]), [])
            if not steps:
                step_label.set_text("No steps")
                curves_plot.update_figure(empty_figure(420))
                similarity_plot.update_figure(empty_figure(260))
                for value_label, change_label in parameter_labels:
                    value_label.set_text("-")
                    change_label.set_text("")
                return
            index = min(selection["index"], len(steps) - 1)
            selection["index"] = index
            slider.props["max"] = max(len(steps) - 1, 1)
            slider.set_value(index)
            step_label.set_text(f"Step {index} of {len(steps) - 1}")
            step = steps[index]
            max_counts = max(float(np.nanmax(s["curves"])) for s in steps)
            figure = heatmap_figure if view_toggle.value == "heatmap" else curves_figure
            curves_plot.update_figure(figure(step, max_counts))
            similarity_plot.update_figure(similarity_figure(steps, index))
            for (value_label, change_label), value, initial in zip(parameter_labels, step["parameters"], steps[0]["parameters"]):
                value_label.set_text(f"{value:.8g}")
                change_label.set_text(f"{value - initial:+.3g} since step 0")

        def update(**changes):
            """apply a change of the selection and redraw"""
            if selection["busy"]:
                return
            selection["busy"] = True
            try:
                selection.update(changes)
                if phase_toggle.value != selection["phase"]:
                    phase_toggle.set_value(selection["phase"])
                if module_select.value != selection["module"]:
                    module_select.set_value(selection["module"])
                render()
            finally:
                selection["busy"] = False

        def user_changed(**changes):
            if not selection["busy"]:
                follow.set_value(False)
                update(**changes)

        def refresh():
            status.set_text(session.status)
            start_button.set_enabled(not session.is_running())
            stop_button.set_enabled(session.is_running())
            num_new_lines = min(session.num_log_lines - selection["log_lines"], len(session.log_lines))
            if num_new_lines > 0:
                for line in list(session.log_lines)[-num_new_lines:]:
                    log.push(line)
            selection["log_lines"] = session.num_log_lines
            if session.version == selection["version"]:
                return
            selection["version"] = session.version
            modules = sorted({key[0] for key in session.steps})
            selection["busy"] = True
            try:
                module_select.set_options({m: f"{m}" for m in modules}, value=selection["module"] if selection["module"] in modules else None)
            finally:
                selection["busy"] = False
            if follow.value and session.last_key is not None:
                update(module=session.last_key[0], phase=session.last_key[1], index=len(session.steps[session.last_key]) - 1)
            elif selection["module"] not in modules:
                update(module=modules[0] if modules else None, index=0)
            else:
                update()

        async def start():
            if session.is_running():
                return
            if not os.path.isdir(data_path.value):
                ui.notify(f"Data folder {data_path.value} does not exist", type="negative")
                return
            arguments = ["--data-path", data_path.value, "--endpoint", defaults.endpoint,
                         "--base-peak-angle", str(base_peak_angle.value), "--base-peak-roi-width", str(roi_width.value)]
            if module.value.strip().lower() != "all":
                if not module.value.strip().isdigit():
                    ui.notify("Module must be 'all' or a module index", type="negative")
                    return
                arguments += ["--module", module.value.strip()]
            if record.value.strip():
                arguments += ["--record", record.value.strip()]
                step_file.set_value(record.value.strip())
            follow.set_value(True)
            await session.start(arguments)

        def load():
            if session.is_running():
                ui.notify("Stop the calibration before opening a step file", type="warning")
                return
            try:
                session.load(step_file.value)
            except (OSError, KeyError, ValueError) as error:
                ui.notify(f"Could not open {step_file.value}: {error}", type="negative")
                return
            follow.set_value(False)
            selection["module"] = None

        start_button.on_click(start)
        stop_button.on_click(session.stop)
        load_button.on_click(load)
        module_select.on_value_change(lambda e: user_changed(module=e.value, index=0))
        phase_toggle.on_value_change(lambda e: user_changed(phase=e.value, index=0))
        view_toggle.on_value_change(lambda e: update())
        slider.on_value_change(lambda e: user_changed(index=int(e.value)))

        ui.timer(0.3, refresh)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-path", default=os.environ.get("ANGCAL_TEST_DATA"), help="folder with the acquisition files (default: $ANGCAL_TEST_DATA)")
    parser.add_argument("--record", default="calibration_steps.h5", help="default file to record the calibration steps to")
    parser.add_argument("--endpoint", default=DEFAULT_ENDPOINT, help="ZMQ endpoint the calibration sends its steps to")
    parser.add_argument("--host", default="127.0.0.1", help="address to serve the GUI on")
    parser.add_argument("--port", type=int, default=8080, help="port to serve the GUI on")
    defaults = parser.parse_args(argv)
    defaults.record = os.path.abspath(defaults.record)

    create_page(defaults)
    app.on_startup(lambda: background_tasks.create(receive_steps(defaults.endpoint)))
    app.on_shutdown(session.stop)
    ui.run(host=defaults.host, port=defaults.port, title="Angle calibration", reload=False, show=False)


if __name__ in {"__main__", "__mp_main__"}:
    main()
