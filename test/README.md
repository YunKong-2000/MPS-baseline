# 测试说明

## 测试文件结构

- `test_particle.cpp`: Particle 类的所有功能测试

## 运行测试

### 使用 CMake 构建

```bash
cd build
cmake ..
make
```

### 运行测试程序

```bash
cd build/bin
./test_particle
```

## 测试内容

1. **FluidParticle 读取功能测试**
   - 从 TXT 文件读取
   - 从 CSV 文件读取
   - 手动调用读取函数

2. **SolidParticle 读取功能测试**
   - 从 TXT 文件读取
   - 从 CSV 文件读取
   - 手动调用读取函数

3. **粒子数据写入功能测试**
   - FluidParticle 写入 TXT/CSV
   - SolidParticle 写入 TXT/CSV（包含法向向量）

4. **粒子数据拷贝功能测试**
   - FluidParticle 拷贝
   - SolidParticle 拷贝

5. **向量写入功能测试（Debug）**
   - double 向量写入
   - double3 向量写入

6. **文件类型自动检测测试**
   - TXT 文件类型识别
   - CSV 文件类型识别

## 测试数据

测试使用的数据文件位于 `data/` 目录：
- `fluid_particles.txt` / `fluid_particles.csv`: 流体粒子数据
- `solid_particles.txt` / `solid_particles.csv`: 固体粒子数据

## 测试输出

测试结果文件会生成在 `build/bin/data/` 目录下：
- `output_fluid.txt` / `output_fluid.csv`
- `output_solid.txt` / `output_solid.csv`
- `debug_double.txt` / `debug_double3.txt` / `debug_double3.csv`

