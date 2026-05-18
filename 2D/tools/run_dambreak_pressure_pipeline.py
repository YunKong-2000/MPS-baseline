#!/usr/bin/env python3
"""Run dam-break pressure probe + comparison plotting in one command."""

from __future__ import annotations

import argparse
import configparser
import subprocess
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "串联执行 dam-break 测压后处理与实验对比绘图。"
            "默认测压方法为 mps_taylor，仅需提供测压点、邻域半径和绘图最大时间。"
        )
    )
    parser.add_argument("--point-x", type=float, required=True, help="测压点 x 坐标。")
    parser.add_argument("--point-y", type=float, required=True, help="测压点 y 坐标。")
    parser.add_argument("--radius", type=float, required=True, help="测压邻域半径。")
    parser.add_argument(
        "--max-time",
        type=float,
        required=True,
        help="模拟结果绘图最大物理时间（秒），用于筛选模拟数据 time <= max-time。",
    )

    parser.add_argument(
        "--method",
        type=str,
        default="mps_taylor",
        choices=("average", "nearest", "mps_taylor", "lsmps"),
        help="测压方法，默认 mps_taylor（lsmps 为兼容别名）。",
    )
    parser.add_argument(
        "--mps-probe-min-neighbors",
        type=int,
        default=3,
        help="mps_taylor 方法下测压点最少邻域粒子数，默认 3。",
    )
    parser.add_argument(
        "--case-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="2D 算例根目录，默认当前脚本上级目录。",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config.ini"),
        help="配置文件路径（相对 --case-dir 或绝对路径），默认 config.ini。",
    )
    parser.add_argument(
        "--dt",
        type=float,
        default=None,
        help=(
            "手动指定相邻 CSV 行的时间间隔；若不提供则优先读取 "
            "config 的 [Simulation] OutputInterval。"
        ),
    )
    parser.add_argument(
        "--vtk",
        type=Path,
        default=None,
        help="vtk 输入目录/文件；默认使用 --result-dir（不存在时回退读取 config 的 [File] OutputDir）。",
    )
    parser.add_argument(
        "--exp",
        type=Path,
        default=None,
        help="实验压力数据文件（可选；仅在提供该参数时才会绘制实验对比图）。",
    )
    parser.add_argument(
        "--result-dir",
        type=Path,
        default=Path("output"),
        help="结果文件目录（用于 vtk 输入、probe csv 与对比图输出默认路径）。",
    )
    parser.add_argument(
        "--probe-csv",
        type=Path,
        default=Path("output/dambreak_pressure.csv"),
        help="测压后处理输出 CSV 路径。",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("output/dambreak_pressure_comparison.png"),
        help="对比图输出路径。",
    )
    parser.add_argument(
        "--sim-curve-out",
        type=Path,
        default=Path("output/probe_pressure_timeseries.png"),
        help="仅模拟压力曲线输出路径（横轴时间，纵轴压力）。",
    )
    parser.add_argument(
        "--pressure-field",
        type=str,
        default="pressure",
        help="VTK 压力字段名，默认 pressure。",
    )
    parser.add_argument(
        "--tools-build-dir",
        type=Path,
        default=Path("tools/build"),
        help="工具构建目录，默认 tools/build。",
    )
    parser.add_argument(
        "--rebuild-tools",
        action="store_true",
        help="强制重新执行 cmake 配置与编译 pressure_probe_postprocess。",
    )
    return parser


def resolve_path(case_dir: Path, path: Path) -> Path:
    if path.is_absolute():
        return path
    return case_dir / path


def read_timestep_from_config(config_path: Path) -> float:
    config = configparser.ConfigParser()
    if not config.read(config_path):
        raise FileNotFoundError(f"无法读取配置文件: {config_path}")
    if not config.has_option("Simulation", "TimeStep"):
        raise ValueError(f"配置文件缺少 [Simulation] TimeStep: {config_path}")
    dt = config.getfloat("Simulation", "TimeStep")
    if dt <= 0.0:
        raise ValueError(f"配置文件 TimeStep 非法（需为正数）: {dt}")
    return dt


def read_output_interval_from_config(config_path: Path) -> float:
    config = configparser.ConfigParser()
    if not config.read(config_path):
        raise FileNotFoundError(f"无法读取配置文件: {config_path}")
    if not config.has_option("Simulation", "OutputInterval"):
        raise ValueError(f"配置文件缺少 [Simulation] OutputInterval: {config_path}")
    output_interval = config.getfloat("Simulation", "OutputInterval")
    if output_interval <= 0.0:
        raise ValueError(
            f"配置文件 OutputInterval 非法（需为正数）: {output_interval}"
        )
    return output_interval


def read_output_dir_from_config(config_path: Path, case_dir: Path) -> Path:
    config = configparser.ConfigParser()
    if not config.read(config_path):
        raise FileNotFoundError(f"无法读取配置文件: {config_path}")
    if not config.has_option("File", "OutputDir"):
        raise ValueError(f"配置文件缺少 [File] OutputDir: {config_path}")
    output_dir = Path(config.get("File", "OutputDir").strip())
    return resolve_path(case_dir, output_dir)


def run_checked(cmd: list[str], cwd: Path) -> None:
    print("[RUN]", " ".join(str(item) for item in cmd))
    subprocess.run(cmd, cwd=cwd, check=True)


def ensure_pressure_probe_tool(case_dir: Path, tools_build_dir: Path, rebuild: bool) -> Path:
    binary = tools_build_dir / "pressure_probe_postprocess"
    if binary.exists() and not rebuild:
        return binary

    run_checked(
        ["cmake", "-S", str(case_dir / "tools"), "-B", str(tools_build_dir)],
        cwd=case_dir,
    )
    run_checked(
        [
            "cmake",
            "--build",
            str(tools_build_dir),
            "--target",
            "pressure_probe_postprocess",
            "-j",
        ],
        cwd=case_dir,
    )
    if not binary.exists():
        raise FileNotFoundError(f"未找到工具可执行文件: {binary}")
    return binary


def main() -> None:
    args = build_parser().parse_args()
    if args.radius <= 0.0:
        raise ValueError("--radius 必须为正数")
    if args.mps_probe_min_neighbors < 1:
        raise ValueError("--mps-probe-min-neighbors 必须 >= 1")
    if args.max_time <= 0.0:
        raise ValueError("--max-time 必须为正数")

    case_dir = args.case_dir.resolve()
    config_path = resolve_path(case_dir, args.config)
    tools_build_dir = resolve_path(case_dir, args.tools_build_dir)
    tools_build_dir.mkdir(parents=True, exist_ok=True)
    result_dir = resolve_path(case_dir, args.result_dir)
    result_dir.mkdir(parents=True, exist_ok=True)

    if args.dt is not None:
        dt = args.dt
    else:
        try:
            dt = read_output_interval_from_config(config_path)
        except ValueError:
            # 保持向后兼容：若旧配置无 OutputInterval，则退化到 TimeStep。
            dt = read_timestep_from_config(config_path)
    if dt <= 0.0:
        raise ValueError("--dt 必须为正数")

    if args.vtk is not None:
        vtk_path = resolve_path(case_dir, args.vtk)
    else:
        vtk_path = result_dir
        if not vtk_path.exists():
            vtk_path = read_output_dir_from_config(config_path, case_dir)
    if not vtk_path.exists():
        raise FileNotFoundError(f"VTK 输入路径不存在: {vtk_path}")

    exp_path = resolve_path(case_dir, args.exp) if args.exp is not None else None
    has_experiment_data = exp_path is not None
    if has_experiment_data and not exp_path.exists():
        raise FileNotFoundError(f"实验数据文件不存在: {exp_path}")

    default_probe_csv = result_dir / "dambreak_pressure.csv"
    probe_csv = (
        resolve_path(case_dir, args.probe_csv)
        if args.probe_csv != Path("output/dambreak_pressure.csv")
        else default_probe_csv
    )
    probe_csv.parent.mkdir(parents=True, exist_ok=True)
    default_out_png = result_dir / "dambreak_pressure_comparison.png"
    out_png = None
    if has_experiment_data:
        out_png = (
            resolve_path(case_dir, args.out)
            if args.out != Path("output/dambreak_pressure_comparison.png")
            else default_out_png
        )
        out_png.parent.mkdir(parents=True, exist_ok=True)
    default_sim_curve_png = result_dir / "probe_pressure_timeseries.png"
    sim_curve_png = (
        resolve_path(case_dir, args.sim_curve_out)
        if args.sim_curve_out != Path("output/probe_pressure_timeseries.png")
        else default_sim_curve_png
    )
    sim_curve_png.parent.mkdir(parents=True, exist_ok=True)

    probe_bin = ensure_pressure_probe_tool(case_dir, tools_build_dir, args.rebuild_tools)
    plot_script = case_dir / "tools" / "plot_dambreak_pressure_comparison.py"
    sim_curve_plot_script = case_dir / "tools" / "plot_probe_pressure_timeseries.py"
    if not sim_curve_plot_script.exists():
        raise FileNotFoundError(f"未找到模拟压力曲线绘图脚本: {sim_curve_plot_script}")
    if has_experiment_data and not plot_script.exists():
        raise FileNotFoundError(f"未找到绘图脚本: {plot_script}")

    probe_cmd = [
        str(probe_bin),
        "--radius",
        str(args.radius),
        "--point",
        str(args.point_x),
        str(args.point_y),
        "--vtk",
        str(vtk_path),
        "--method",
        "mps_taylor" if args.method == "lsmps" else args.method,
        "--mps-probe-min-neighbors",
        str(args.mps_probe_min_neighbors),
        "--pressure-field",
        args.pressure_field,
        "--output",
        str(probe_csv),
    ]
    run_checked(probe_cmd, cwd=case_dir)

    sim_curve_cmd = [
        sys.executable,
        str(sim_curve_plot_script),
        "--csv",
        str(probe_csv),
        "--dt",
        str(dt),
        "--max-time",
        str(args.max_time),
        "--out",
        str(sim_curve_png),
        "--title",
        "Probe Pressure Time-Series",
    ]
    run_checked(sim_curve_cmd, cwd=case_dir)

    if has_experiment_data:
        plot_cmd = [
            sys.executable,
            str(plot_script),
            "--exp",
            str(exp_path),
            "--sim",
            str(probe_csv),
            "--dt",
            str(dt),
            "--sim-time-col",
            "step_id",
            "--sim-pressure-col",
            "probe_0_pressure",
            "--max-time",
            str(args.max_time),
            "--out",
            str(out_png),
        ]
        run_checked(plot_cmd, cwd=case_dir)
    else:
        print("未检测到实验数据输入，已跳过实验-模拟对比绘图，仅输出模拟测压结果。")

    print("\nPipeline 完成：")
    print(f"- 测压点: ({args.point_x}, {args.point_y})")
    print(f"- 测压半径: {args.radius}")
    print(f"- 方法: {args.method}")
    print(f"- mps_taylor最少邻域粒子数: {args.mps_probe_min_neighbors}")
    print(f"- 时间间隔(dt): {dt}")
    print(f"- 结果目录: {result_dir}")
    if has_experiment_data:
        print(f"- 实验数据文件: {exp_path}")
    else:
        print("- 实验数据文件: 未提供")
    print(f"- VTK 输入: {vtk_path}")
    print(f"- 测压 CSV: {probe_csv}")
    print(f"- 模拟压力曲线: {sim_curve_png}")
    if has_experiment_data:
        print(f"- 对比图: {out_png}")
    else:
        print("- 对比图: 已跳过（未提供实验数据）")


if __name__ == "__main__":
    main()
