#!/usr/bin/env python3
"""Plot dam-break pressure comparison (experiment vs simulation)."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt


def resolve_exp_path(path: Path) -> Path:
    if path.exists():
        return path

    # Compatible with accidental double extension, e.g. dambreak_pressure.txt.txt
    alt_path = Path(str(path) + ".txt")
    if alt_path.exists():
        return alt_path

    raise FileNotFoundError(f"Experiment file not found: {path}")


def read_experiment_txt(path: Path) -> Tuple[List[float], List[float]]:
    times: List[float] = []
    pressures: List[float] = []

    with path.open("r", encoding="utf-8") as file:
        for line_no, raw_line in enumerate(file, start=1):
            line = raw_line.strip()
            if not line:
                continue

            parts = line.replace(",", " ").split()
            if len(parts) < 2:
                continue

            try:
                t = float(parts[0])
                p = float(parts[1])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid experiment pressure data at line {line_no}: {raw_line.rstrip()}"
                ) from exc

            times.append(t)
            pressures.append(p)

    if not times:
        raise ValueError(f"No valid experiment data found in: {path}")

    return times, pressures


def read_simulation_csv(
    path: Path, time_col: str, pressure_col: str
) -> Tuple[List[float], List[float]]:
    times: List[float] = []
    pressures: List[float] = []

    with path.open("r", encoding="utf-8", newline="") as file:
        reader = csv.DictReader(file)
        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {path}")
        if time_col not in reader.fieldnames:
            raise ValueError(f"Column '{time_col}' not found in {path}")
        if pressure_col not in reader.fieldnames:
            raise ValueError(f"Column '{pressure_col}' not found in {path}")

        for row in reader:
            t_raw = (row.get(time_col) or "").strip()
            p_raw = (row.get(pressure_col) or "").strip()
            if not t_raw or not p_raw:
                continue
            if p_raw.lower() == "nan":
                continue

            try:
                t = float(t_raw)
                p = float(p_raw)
            except ValueError:
                continue

            times.append(t)
            pressures.append(p)

    if not times:
        raise ValueError(f"No valid simulation data found in: {path}")

    return times, pressures


def filter_by_max_time(
    times: List[float], values: List[float], max_time: float
) -> Tuple[List[float], List[float]]:
    filtered_t: List[float] = []
    filtered_v: List[float] = []
    for t, v in zip(times, values):
        if t <= max_time:
            filtered_t.append(t)
            filtered_v.append(v)
    return filtered_t, filtered_v


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Plot comparison of dam-break pressure (experiment vs simulation)."
    )
    parser.add_argument(
        "--exp",
        type=Path,
        default=Path("tools/data/dambreak_pressure.txt"),
        help="Experiment pressure TXT path (col1=time, col2=pressure).",
    )
    parser.add_argument(
        "--sim",
        type=Path,
        required=True,
        help="Simulation pressure CSV path.",
    )
    parser.add_argument(
        "--dt",
        type=float,
        required=True,
        help="Simulation time-step size in seconds.",
    )
    parser.add_argument(
        "--sim-time-col",
        type=str,
        default="step_id",
        help="Simulation step/time column name in CSV.",
    )
    parser.add_argument(
        "--sim-pressure-col",
        type=str,
        default="probe_0_pressure",
        help="Simulation pressure column name in CSV.",
    )
    parser.add_argument(
        "--time-scale",
        type=float,
        default=math.sqrt(9.8 / 0.6),
        help="Scale factor multiplied to simulation time for alignment.",
    )
    parser.add_argument(
        "--max-time",
        type=float,
        default=2.0,
        help="Maximum time shown in figure (seconds).",
    )
    parser.add_argument(
        "--p0",
        type=float,
        default=10000.0,
        help="Reference pressure P0 for normalization (Pa).",
    )
    parser.add_argument(
        "--normalize-exp-with-p0",
        action="store_true",
        help="Normalize experimental pressure by P0 (use only when experiment data is in Pa).",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("output/dambreak_pressure_comparison.png"),
        help="Output figure path.",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="Dam-break pressure comparison",
        help="Figure title.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.dt <= 0.0:
        raise ValueError("--dt must be positive")
    if args.time_scale <= 0.0:
        raise ValueError("--time-scale must be positive")
    if args.max_time <= 0.0:
        raise ValueError("--max-time must be positive")
    if args.p0 <= 0.0:
        raise ValueError("--p0 must be positive")

    exp_path = resolve_exp_path(args.exp)
    exp_time, exp_pressure = read_experiment_txt(exp_path)

    sim_step, sim_pressure = read_simulation_csv(
        args.sim, time_col=args.sim_time_col, pressure_col=args.sim_pressure_col
    )
    sim_time_raw = [step * args.dt for step in sim_step]
    sim_time_raw, sim_pressure = filter_by_max_time(
        sim_time_raw, sim_pressure, args.max_time
    )
    sim_time = [t * args.time_scale for t in sim_time_raw]

    if not exp_time:
        raise ValueError("No experiment data found.")
    if not sim_time:
        raise ValueError("No simulation points left after max-time filtering.")

    if args.normalize_exp_with_p0:
        exp_pressure = [p / args.p0 for p in exp_pressure]
    sim_pressure = [p / args.p0 for p in sim_pressure]

    args.out.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(9.5, 5.2))
    plt.plot(
        exp_time,
        exp_pressure,
        "o-",
        linewidth=1.8,
        markersize=4,
        label="Experiment",
    )
    plt.plot(sim_time, sim_pressure, "-", linewidth=2.0, label="Simulation")
    plt.xlabel("Time (s)")
    plt.ylabel("P/P0")
    plt.title(args.title)
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.out, dpi=220)
    print(f"Saved figure to: {args.out}")
    print(
        "Applied simulation time transform: "
        f"t_plot = step * dt * time_scale = step * {args.dt} * {args.time_scale}"
    )
    print(f"Applied simulation raw-time filter: step * dt <= {args.max_time}")
    print(f"Applied pressure normalization: P/P0, P0 = {args.p0} Pa")


if __name__ == "__main__":
    main()

