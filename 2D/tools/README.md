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

用于从输出目录下的 `vtk` 文件提取液面前沿位置（`surface_type` 为 `surface/splash` 的粒子中 `x` 最大者）。

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
- `--show`：保存后显示交互窗口（可选）

### 图像说明

- 横轴：`Time (s)`
- 纵轴：`Pressure (Pa)`
- 每条曲线：一个观测点（对应一个压力列）

## 压力对比绘图（实验 vs 模拟）

用于将实验压力数据与模拟压力数据画在同一张图中对比。  
其中模拟时间坐标按下面公式对齐实验坐标：

- `t_plot = step_id * dt * sqrt(9.8 / 0.6)`

默认只绘制 `t <= 2s` 的数据。

纵坐标默认绘制为 `P/P0`：

- 模拟压力会自动除以 `P0`
- 实验压力默认认为已经是归一化值，不再除以 `P0`

### 运行示例

```bash
cd /home/amax/mps-baseline/2D
python3 tools/plot_dambreak_pressure_comparison.py \
  --exp tools/data/dambreak_pressure.txt \
  --sim output/dambreak_pressure.csv \
  --dt 0.001 \
  --sim-time-col step_id \
  --sim-pressure-col probe_0_pressure \
  --p0 10000 \
  --max-time 2.0 \
  --out output/dambreak_pressure_comparison.png
```

### 常用参数

- `--exp`：实验压力 txt 路径（默认 `tools/data/dambreak_pressure.txt`）
- `--sim`：模拟压力 csv 路径（必填）
- `--dt`：模拟时间步长（必填）
- `--sim-time-col`：模拟时间/步编号列名（默认 `step_id`）
- `--sim-pressure-col`：模拟压力列名（默认 `probe_0_pressure`）
- `--time-scale`：时间缩放系数（默认 `sqrt(9.8/0.6)`）
- `--max-time`：模拟原始时间筛选上限（`step_id * dt <= max-time`，默认 `2.0` 秒）
- `--p0`：归一化参考压力 `P0`（Pa，默认 `10000`）
- `--normalize-exp-with-p0`：若实验压力原始单位是 Pa，可开启该选项按 `P0` 归一化
- `--out`：输出图片路径

