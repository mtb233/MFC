#!/usr/bin/env python3
"""Plot single-bubble Lagrangian evolution (D/lag_bubble_evol_*.dat) from a JSON config.

Usage (from the MFC root): python projects/kernel_study/post_process/plot_bubble_evolution.py config.json
Relative paths in the config resolve against the current working directory.
Config format: bubble_plots.schema.json (next to this script).
"""
import argparse
import json
import sys
from pathlib import Path

import jsonschema
import matplotlib.pyplot as plt
import pandas as pd

SCHEMA_PATH = Path(__file__).with_name("bubble_plots.schema.json")

# Defaults mirror generate_cases.py
DEFAULT_REFERENCE = {"x_0": 1.0e-3, "rho_0": 1.0e3, "c_0": 1475.0}

# property -> (dimension, display name, symbol for normalized label)
PROPERTIES = {
    "currentTime": ("time", "Time", "t"),
    "x": ("length", "x", "x"),
    "y": ("length", "y", "y"),
    "z": ("length", "z", "z"),
    "radius": ("length", "Radius", "R"),
    "interfaceVelocity": ("velocity", "Interface velocity", r"\dot{R}"),
    "corePressure": ("pressure", "Core pressure", "p_b"),
    "coreVaporMass": ("mass", "Core vapor mass", "m_v"),
    "coreVaporConcentration": ("dimensionless", "Core vapor concentration", "C_v"),
}

# dimension -> {unit: SI value of one unit}, first entry is the default unit
UNITS = {
    "time": {"us": 1e-6, "s": 1.0, "ms": 1e-3, "ns": 1e-9},
    "length": {"um": 1e-6, "m": 1.0, "mm": 1e-3},
    "velocity": {"m/s": 1.0},
    "pressure": {"kPa": 1e3, "Pa": 1.0, "MPa": 1e6, "atm": 101325.0},
    "mass": {"ng": 1e-12, "kg": 1.0, "pg": 1e-15},
    "dimensionless": {"-": 1.0},
}
UNIT_LABELS = {"us": "μs", "um": "μm"}


def si_scales(ref):
    x_0, rho_0, c_0 = ref["x_0"], ref["rho_0"], ref["c_0"]
    return {
        "time": x_0 / c_0,
        "length": x_0,
        "velocity": c_0,
        "pressure": rho_0 * c_0**2,
        "mass": rho_0 * x_0**3,
        "dimensionless": 1.0,
    }


def find_bubble_file(case_dir):
    files = sorted((case_dir / "D").glob("lag_bubble_evol_*.dat"))
    with_data = []
    for f in files:
        with f.open() as fh:
            if sum(1 for _ in fh) > 1:
                with_data.append(f)
    if len(with_data) != 1:
        sys.exit(f"{case_dir}: expected exactly one D/lag_bubble_evol_*.dat with data, found {len(with_data)} (of {len(files)} files)")
    return with_data[0]


def load_case(case_dir):
    path = find_bubble_file(case_dir)
    data = pd.read_csv(path, sep=r"\s+")
    if data["particleID"].nunique() != 1:
        sys.exit(f"{path}: contains more than one bubble; only single-bubble cases are supported")
    return data


def axis_values(data, axis, scales, time_offset):
    prop = axis["property"]
    dim, name, symbol = PROPERTIES[prop]
    unit = axis.get("unit", next(iter(UNITS[dim])))
    values = data[prop] * scales[dim]
    if unit == "normalized":
        values = values / values.iloc[0]
        default_label = rf"${symbol}/{symbol}_0$"
    elif unit in UNITS[dim]:
        values = values / UNITS[dim][unit]
        if prop == "currentTime":
            values = values + time_offset
        default_label = name if dim == "dimensionless" else f"{name} ({UNIT_LABELS.get(unit, unit)})"
    else:
        sys.exit(f"unit '{unit}' is not valid for '{prop}' ({dim}); choose from {['normalized', *UNITS[dim]]}")
    return values, axis.get("label", default_label)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", type=Path, help="JSON plot configuration")
    args = parser.parse_args()

    config = json.loads(args.config.read_text())
    try:
        jsonschema.validate(config, json.loads(SCHEMA_PATH.read_text()))
    except jsonschema.ValidationError as e:
        sys.exit(f"{args.config}: invalid config at {'/'.join(map(str, e.absolute_path)) or '<root>'}: {e.message}")

    scales = si_scales({**DEFAULT_REFERENCE, **config.get("reference", {})})
    out_dir = Path(config.get("output_dir", "projects/kernel_study/results"))
    out_dir.mkdir(parents=True, exist_ok=True)
    fmt = config.get("format", "png")

    cache = {}
    for plot in config["plots"]:
        fig, ax = plt.subplots(figsize=plot.get("figsize", (5.5, 4.0)))
        x_axis = plot.get("x", {"property": "currentTime"})
        for curve in plot["curves"]:
            case_dir = Path(curve["case"]).resolve()
            if case_dir not in cache:
                cache[case_dir] = load_case(case_dir)
            data = cache[case_dir]
            offset = curve.get("time_offset", 0.0)
            x, xlabel = axis_values(data, x_axis, scales, offset)
            y, ylabel = axis_values(data, plot["y"], scales, offset)
            ax.plot(x, y, label=curve.get("label", case_dir.name), **curve.get("style", {}))

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xscale(x_axis.get("scale", "linear"))
        ax.set_yscale(plot["y"].get("scale", "linear"))
        if "xlim" in plot:
            ax.set_xlim(plot["xlim"])
        if "ylim" in plot:
            ax.set_ylim(plot["ylim"])
        if "title" in plot:
            ax.set_title(plot["title"])
        if plot.get("legend", True):
            ax.legend(fontsize="small")
        fig.tight_layout()
        out = out_dir / f"{plot['name']}.{fmt}"
        fig.savefig(out, dpi=config.get("dpi", 300))
        plt.close(fig)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
