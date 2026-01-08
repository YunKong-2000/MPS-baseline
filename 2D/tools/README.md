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

