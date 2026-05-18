# 算例生成工具

## 生成静水问题算例

### 快速生成（约2万粒子）

```bash
cd /home/amax/mps-baseline/2D
g++ -std=c++17 -O2 tools/gen_hydrostatic_simple.cpp -o tools/gen_hydrostatic_simple
./tools/gen_hydrostatic_simple
```

### 自定义粒子间距

```bash
./tools/gen_hydrostatic_simple 0.015  # 使用0.015m的粒子间距
```

### 生成的算例参数

- **容器宽度**: 2.0 m
- **容器高度**: 2.0 m  
- **水位高度**: 1.0 m
- **粒子间距**: 0.01 m（默认）
- **流体粒子数**: 20,301
- **固体粒子数**: 603
- **总粒子数**: 20,904

### 输出文件

- `data/fluid_particles_2d.txt` - 流体粒子文件
- `data/solid_particles_2d.txt` - 固体粒子文件

### 配置文件

生成算例后，需要确保 `config.ini` 中的粒子参数与生成的算例匹配：

```ini
[Particle]
ParticleSpacing = 0.01
ParticleRadius = 0.005
SmoothingRadius = 0.021
CellSize = 0.042
```

### 运行模拟

```bash
cd build/bin
./MPSBaseline2D
```

注意：如果遇到PETSc MPI版本问题，需要配置正确的MPI环境。

## dam-break 前沿后处理

用于从输出目录下的 `vtk` 文件提取液面前沿位置（`surface_type` 为 `surface` 的自由面粒子中 `x` 最大者）。

### 运行示例

```bash
cd /home/amax/mps-baseline/2D
./tools/build/dambreak_front_postprocess ./output
```

指定输出文件和时间参数：

```bash
./tools/build/dambreak_front_postprocess \
  ./output \
  ./output/dambreak_front_position.csv \
  0.001
```

### 输入参数

- 第1个参数：`vtk` 目录（必填）
- 第2个参数：输出 CSV 文件路径（可选，默认 `<vtk目录>/dambreak_front_position.csv`）
- 第3个参数：时间步大小 `dt`（可选，默认 `1.0`）

### 输出列

- `time`：时间（`time = step_id * dt`，其中 `step_id` 由文件名尾部解析，解析失败则按排序顺序编号）
- `front_x`：前沿 x 坐标

## 液面高度后处理（按观测点 x）

用于从一个或多个 `vtk` 结果中提取每个时刻、各观测点 `x` 处的液面高度（定义为该 `x` 附近流体粒子的最大 `y`）。

### 编译

```bash
cd /home/amax/mps-baseline/2D
cmake -S tools -B tools/build
cmake --build tools/build --target liquid_height_postprocess -j
```

### 运行示例

```bash
./tools/build/liquid_height_postprocess \
  --vtk ./output \
  --x 0.50 \
  --x 0.70 \
  --dt 0.001 \
  --x-tolerance 0.01 \
  --output ./output/liquid_height.csv
```

也可通过文件输入观测点 `x`（每行一个值，支持 `#` 注释）：

```bash
./tools/build/liquid_height_postprocess \
  --vtk ./output \
  --x-file ./tools/data/dambreak/x_probes.txt \
  --dt 0.001 \
  --x-tolerance 0.01 \
  --output ./output/liquid_height.csv
```

### 输入参数

- `--vtk`：vtk 文件或目录，可重复指定（程序仅处理文件名以 `result` 开头的 `.vtk`）
- `--x`：观测点 `x` 坐标，可重复指定
- `--x-file`：观测点 `x` 文件，每行一个 `x`
- `--dt`：相邻输出帧时间间隔（默认 `1.0`）
- `--x-tolerance`：`x` 匹配容差，满足 `|x_i - x_probe| <= tol` 视为命中（默认 `1e-6`）
- `--output`：输出 CSV 路径（默认 `liquid_height.csv`）

### 输出格式

- 第1列：`time`（`time = step_id * dt`，其中 `step_id` 从文件名尾部数字解析，解析失败则按排序顺序）
- 第2列起：`probe_i_height`（第 `i` 个观测点的液面高度，即命中粒子 `y` 的最大值）
- 若某时刻某观测点未命中任何粒子，高度记为 `0`
- 为保证物理意义，高度最小值为 `0`，负值会被截断到 `0`

## 前沿对比绘图程序（Python）

用于将实验前沿数据与模拟前沿数据绘制到同一张图中进行对比。

### 输入数据格式

- 实验数据文件（`--exp`）：文本文件，每行至少两列，分别为 `time front_x`
- 模拟数据文件（`--sim`）：CSV 文件，默认读取列名 `time` 和 `front_x`

### 运行示例

```bash
cd /home/amax/mps-baseline/2D
python3 tools/plot_dambreak_front_comparison.py \
  --exp tools/data/dambreak前沿\ \(2\).txt \
  --sim output/dambreak_front.csv \
  --sim-front-col front_x \
  --out output/dambreak_front_comparison.png \
  --title "Dam-break front comparison"
```

### 常用参数

- `--exp`：实验数据路径（必填）
- `--sim`：模拟数据路径（必填）
- `--sim-time-col`：模拟数据时间列名（默认 `time`）
- `--sim-front-col`：模拟数据前沿列名（默认 `front_x`）
- `--max-time`：仅绘制 `time <= max-time` 的数据（可选）
- `--out`：输出图片路径（默认 `dambreak_front_comparison.png`）

## 压力观测点后处理

用于从一个或多个 `vtk` 结果文件中提取给定观测点的压力时间序列。

### 编译

```bash
cd /home/amax/mps-baseline/2D
cmake -S tools -B tools/build
cmake --build tools/build --target pressure_probe_postprocess dambreak_front_postprocess -j
```

### 运行示例

```bash
./tools/build/pressure_probe_postprocess \
  --radius 0.03 \
  --point 0.50 0.25 \
  --point 0.70 0.25 \
  --vtk ./output \
  --method average \
  --pressure-field pressure \
  --output ./output/probe_pressure.csv
```

### 输入参数

- `--radius`：观测点邻域半径（`average` 与 `nearest` 都使用）
- `--mps-probe-min-neighbors`：`mps_taylor` 下测压点最少邻域粒子数（默认 `3`，少于该值时测压点压力置 `0`）
- `--point x y`：观测点坐标，可重复指定
- `--points-file`：观测点文件，每行 `x y`（支持 `#` 注释）
- `--vtk`：vtk 文件或目录，可重复指定
- `--pressure-field`：压力标量字段名（默认 `pressure`）
- `--method`：测量方法，支持：
  - `average`：半径内流体粒子压力平均；若邻域内无粒子则压力为 `0`
  - `nearest`：先检查半径邻域；若邻域内无粒子则压力为 `0`，否则取邻域内最近流体粒子的压力
- `--output`：输出 CSV 路径（默认 `pressure_probe.csv`）

### 输出格式

输出为 CSV，行表示不同时间步（按 vtk 文件顺序），列表示不同观测点压力：

- `step_id`：时间步序号
- `probe_i_pressure`：第 `i` 个观测点压力

## 压力曲线绘图（Python）

用于将 `pressure_probe_postprocess` 输出的 CSV 绘制成压力-时间曲线图。

### 运行示例

```bash
cd /home/amax/mps-baseline/2D
python3 tools/plot_probe_pressure_timeseries.py \
  --csv output/dambreak_pressure.csv \
  --dt 0.001 \
  --out output/dambreak_pressure_curve.png \
  --title "Dam-break probe pressure"
```

### 输入参数

- `--csv`：压力后处理输出 CSV 路径
- `--dt`：时间步长（秒），程序据此构造时间轴
- `--out`：输出图片路径（默认 `probe_pressure_timeseries.png`）
- `--title`：图标题（可选）
- `--max-time`：模拟数据最大物理时间筛选上限（秒，可选）
- `--show`：保存后显示交互窗口（可选）

### 图像说明

- 横轴：`Time (s)`
- 纵轴：`Pressure (Pa)`
- 每条曲线：一个观测点（对应一个压力列）

## 液面高度曲线绘图（Python）

用于将 `liquid_height_postprocess` 输出的 CSV 绘制成液面高度-时间曲线图。

### 运行示例

```bash
cd /home/amax/mps-baseline/2D
python3 tools/plot_liquid_height_timeseries.py \
  --csv output/liquid_height.csv \
  --out output/liquid_height_timeseries.png \
  --title "Liquid height time-series"
```

与实验数据对比绘图（实验 `txt/csv`）：

```bash
python3 tools/plot_liquid_height_timeseries.py \
  --csv output/liquid_height.csv \
  --exp tools/data/dambreak/H2-height.txt \
  --max-time 2.0 \
  --out output/liquid_height_comparison.png \
  --title "Liquid height comparison"
```

指定列名、实验 CSV 列名与最大时间筛选：

```bash
python3 tools/plot_liquid_height_timeseries.py \
  --csv output/liquid_height.csv \
  --time-col time \
  --height-cols probe_0_height,probe_1_height \
  --exp tools/data/dambreak/height_exp.csv \
  --exp-time-col time \
  --exp-height-col height \
  --max-time 1.0 \
  --out output/liquid_height_timeseries.png
```

### 输入参数

- `--csv`：液面高度后处理输出 CSV 路径（必填）
- `--time-col`：时间列名（默认 `time`）
- `--height-cols`：要绘制的高度列（逗号分隔，默认自动选择所有 `*_height` 列）
- `--exp`：实验数据路径（可选，支持 `txt/csv`；提供后会叠加实验曲线）
- `--exp-time-col`：实验 CSV 时间列名（默认 `time`，仅 CSV 生效）
- `--exp-height-col`：实验 CSV 高度列名（默认 `height`，仅 CSV 生效）
- `--exp-label`：实验曲线图例名（默认 `Experiment`）
- `--max-time`：仅保留 `time <= max-time` 的数据（可选）
- `--out`：输出图片路径（默认 `liquid_height_timeseries.png`）
- `--title`：图标题（可选）
- `--show`：保存后显示交互窗口（可选）

### 图像说明

- 横轴：`Time (s)`
- 纵轴：`Liquid height (m)`
- 每条曲线：一个观测点高度列
- 若提供 `--exp`：同图叠加实验散点用于对比

## 压力对比绘图（实验 vs 模拟）

用于将实验压力数据与模拟压力数据画在同一张图中对比（不做无量纲化）。

- 横轴时间直接使用物理时间：`time = step_id * dt`
- 纵轴压力直接使用原始压力值（Pa）

可选参数 `--max-time` 用于筛选模拟数据 `time <= max-time`（不筛选实验曲线）。

### 运行示例

```bash
cd /home/amax/mps-baseline/2D
python3 tools/plot_dambreak_pressure_comparison.py \
  --exp tools/data/dambreak_pressure.txt \
  --sim output/dambreak_pressure.csv \
  --dt 0.001 \
  --sim-time-col step_id \
  --sim-pressure-col probe_0_pressure \
  --max-time 1.0 \
  --out output/dambreak_pressure_comparison.png
```

### 常用参数

- `--exp`：实验压力 txt 路径（必填）
- `--sim`：模拟压力 csv 路径（必填）
- `--dt`：模拟时间步长（必填）
- `--sim-time-col`：模拟时间/步编号列名（默认 `step_id`）
- `--sim-pressure-col`：模拟压力列名（默认 `probe_0_pressure`）
- `--max-time`：模拟数据物理时间筛选上限（秒，可选）
- `--out`：输出图片路径

## 一键串联：测压到对比图

脚本：`tools/run_dambreak_pressure_pipeline.py`  
功能：自动执行以下两步：

1. `pressure_probe_postprocess`：从 `vtk` 结果提取测压点压力序列
2. `plot_dambreak_pressure_comparison.py`：绘制实验 vs 仿真压力对比图

默认测压方法是 `mps_taylor`（原始MPS梯度 + 最近邻泰勒外推）。

### 最小用法（仅 4 个必填参数）

```bash
cd /home/amax/mps-baseline/2D
python3 tools/run_dambreak_pressure_pipeline.py \
  --point-x 0.50 \
  --point-y 0.12 \
  --radius 0.03 \
  --max-time 2.0
```

### 默认行为

- 自动读取 `config.ini` 中 `[Simulation] OutputInterval` 作为 `dt`（用于 `step_id -> time`）
- 自动读取 `config.ini` 中 `[File] OutputDir` 作为 `vtk` 输入目录
- 默认测压输出：`output/dambreak_pressure.csv`
- 默认模拟压力曲线输出：`output/probe_pressure_timeseries.png`（横轴时间 `time=step_id*dt`，并按 `--max-time` 截断模拟数据）
- 仅在显式提供 `--exp` 时绘制实验-模拟对比图（默认输出 `output/dambreak_pressure_comparison.png`）
- 测压方法默认 `mps_taylor`

### 常用可选参数

- `--config`：指定配置文件（默认 `config.ini`）
- `--dt`：手动指定相邻输出帧时间间隔（覆盖配置读取）
- `--max-time`：模拟结果绘图最大物理时间（秒，对模拟曲线与对比图中的模拟数据生效）
- `--result-dir`：结果目录（默认 `output`，用于 vtk 输入、测压 CSV、对比图默认输出）
- `--vtk`：手动指定 vtk 目录/文件（覆盖默认结果目录）
- `--exp`：手动指定实验数据完整路径（仅提供该参数时才进行实验对比绘图）
- `--method`：`average` / `nearest` / `mps_taylor`（默认 `mps_taylor`，兼容别名 `lsmps`）
- `--mps-probe-min-neighbors`：`mps_taylor` 下测压点最少邻域粒子数（默认 `3`）
- `--out`：指定输出图片路径
- `--sim-curve-out`：指定模拟压力曲线输出图片路径
- `--probe-csv`：指定测压结果 CSV 路径
- `--rebuild-tools`：强制重新编译 `pressure_probe_postprocess`

