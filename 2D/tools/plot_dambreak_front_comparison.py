#!/usr/bin/env python3
"""Plot dam-break front position comparison (experiment vs simulation)."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt


def read_experiment_txt(path: Path) -> Tuple[List[float], List[float]]:
    """Read experiment text file: col1=time, col2=front position."""
    times: List[float] = []
    fronts: List[float] = []

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
                x = float(parts[1])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid experiment data at line {line_no}: {raw_line.rstrip()}"
                ) from exc

            times.append(t)
            fronts.append(x)

    if not times:
        raise ValueError(f"No valid experiment data found in: {path}")

    return times, fronts


def read_simulation_csv(
    path: Path, time_col: str = "time", front_col: str = "front_x"
) -> Tuple[List[float], List[float]]:
    """Read simulation CSV produced by postprocess tool."""
    times: List[float] = []
    fronts: List[float] = []

    with path.open("r", encoding="utf-8", newline="") as file:
        reader = csv.DictReader(file)
        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {path}")
        if time_col not in reader.fieldnames:
            raise ValueError(f"Column '{time_col}' not found in {path}")
        if front_col not in reader.fieldnames:
            raise ValueError(f"Column '{front_col}' not found in {path}")

        for row in reader:
            t_raw = (row.get(time_col) or "").strip()
            x_raw = (row.get(front_col) or "").strip()
            if not t_raw or not x_raw:
                continue
            if x_raw.lower() == "nan":
                continue

            try:
                t = float(t_raw)
                x = float(x_raw)
            except ValueError:
                continue

            times.append(t)
            fronts.append(x)

    if not times:
        raise ValueError(f"No valid simulation data found in: {path}")

    return times, fronts


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Plot comparison of dam-break front position."
    )
    parser.add_argument(
        "--exp",
        type=Path,
        required=True,
        help="Experiment TXT path (col1=time, col2=front_x).",
    )
    parser.add_argument(
        "--sim",
        type=Path,
        required=True,
        help="Simulation CSV path (default columns: time/front_x).",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("dambreak_front_comparison.png"),
        help="Output figure path.",
    )
    parser.add_argument(
        "--sim-time-col",
        type=str,
        default="time",
        help="Time column name in simulation CSV.",
    )
    parser.add_argument(
        "--sim-front-col",
        type=str,
        default="front_x",
        help="Front position column name in simulation CSV.",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="Dam-break front position comparison",
        help="Figure title.",
    )
    parser.add_argument(
        "--max-time",
        type=float,
        default=None,
        help="Only keep data with time <= max-time (seconds).",
    )
    return parser


def filter_by_max_time(
    times: List[float], values: List[float], max_time: float | None
) -> Tuple[List[float], List[float]]:
    if max_time is None:
        return times, values

    filtered_t: List[float] = []
    filtered_v: List[float] = []
    for t, v in zip(times, values):
        if t <= max_time:
            filtered_t.append(t)
            filtered_v.append(v)
    return filtered_t, filtered_v


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    exp_time, exp_front = read_experiment_txt(args.exp)
    sim_time, sim_front = read_simulation_csv(
        args.sim, time_col=args.sim_time_col, front_col=args.sim_front_col
    )

    exp_time, exp_front = filter_by_max_time(exp_time, exp_front, args.max_time)
    sim_time, sim_front = filter_by_max_time(sim_time, sim_front, args.max_time)
    if not exp_time:
      raise ValueError("No experiment points left after max-time filtering.")
    if not sim_time:
      raise ValueError("No simulation points left after max-time filtering.")

    args.out.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(9, 5))
    plt.plot(exp_time, exp_front, "o-", linewidth=1.8, markersize=4, label="Experiment")
    plt.plot(sim_time, sim_front, "-", linewidth=2.0, label="Simulation")
    plt.xlabel("Time (s)")
    plt.ylabel("Front position x (m)")
    plt.title(args.title)
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.out, dpi=200)
    print(f"Saved figure to: {args.out}")


if __name__ == "__main__":
    main()
