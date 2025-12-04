# 测试说明

本文档说明项目中所有测试程序的用途和使用方法。

## 测试程序列表

项目包含以下测试程序：

1. **test_particle** - Particle 类和 FileOperator 类功能测试
2. **test_config** - 配置文件读取功能测试

---

## test_particle - 粒子类功能测试

### 用途

测试 Particle 类及其派生类（FluidParticle、SolidParticle）和 FileOperator 类的所有核心功能，包括文件读取、写入、数据拷贝等操作。

### 运行方法

```bash
cd build/bin
./test_particle
```

### 测试内容

#### 1. FluidParticle 读取功能测试
- 从 TXT 文件读取流体粒子数据
- 从 CSV 文件读取流体粒子数据
- 手动调用 `getParticleFromFile()` 函数读取

#### 2. SolidParticle 读取功能测试
- 从 TXT 文件读取固体粒子数据（包含法向向量）
- 从 CSV 文件读取固体粒子数据（包含法向向量）
- 手动调用 `getParticleFromFile()` 函数读取

#### 3. 粒子数据写入功能测试
- FluidParticle 数据写入 TXT 文件
- FluidParticle 数据写入 CSV 文件
- SolidParticle 数据写入 TXT 文件（包含法向向量）
- SolidParticle 数据写入 CSV 文件（包含法向向量）

#### 4. 粒子数据拷贝功能测试
- FluidParticle 数据拷贝（验证粒子数量匹配）
- SolidParticle 数据拷贝（验证粒子数量匹配）

#### 5. 向量写入功能测试（Debug）
- double 向量写入 TXT 文件
- double3 向量写入 TXT 文件
- double3 向量写入 CSV 文件

#### 6. 文件类型自动检测测试
- TXT 文件类型识别和读取
- CSV 文件类型识别和读取

### 测试数据

测试使用的数据文件位于 `data/` 目录：
- `fluid_particles.txt` / `fluid_particles.csv`: 流体粒子数据（每行包含：位置x y z 速度x y z）
- `solid_particles.txt` / `solid_particles.csv`: 固体粒子数据（每行包含：位置x y z 速度x y z 法向向量x y z）

### 测试输出

测试结果文件会生成在 `build/bin/data/` 目录下：
- `output_fluid.txt` / `output_fluid.csv` - FluidParticle 写入结果
- `output_solid.txt` / `output_solid.csv` - SolidParticle 写入结果（包含法向向量）
- `debug_double.txt` - double 向量调试输出
- `debug_double3.txt` / `debug_double3.csv` - double3 向量调试输出

---

## test_config - 配置文件读取功能测试

### 用途

测试配置文件读取功能，包括 SimpleIni 库的基本功能和 MPSConfig 参数类的加载与验证。

### 运行方法

```bash
cd build/bin
./test_config
```

### 测试内容

#### 1. SimpleIni 基本功能测试
- 配置文件加载
- 字符串值读取（`GetValue()`）
- 整数值读取（`GetIntValue()`）
- 浮点数值读取（`GetDoubleValue()`）
- 键存在性检查（`HasKey()`）

#### 2. MPSConfig 参数类测试
- 从配置文件加载参数（`LoadFromConfig()`）
- 参数有效性验证（`Validate()`）
- 文件配置读取（输入/输出目录、文件名等）
- 仿真配置读取（时间步长、总时间、密度、粘度、重力等）
- 粒子配置读取（粒子半径、平滑半径、粒子数量等）

### 测试数据

测试使用的配置文件：
- `config.ini` - MPS 仿真配置文件，包含以下部分：
  - **[File]**: 文件路径配置
  - **[Simulation]**: 仿真参数配置
  - **[Particle]**: 粒子参数配置

### 测试输出

测试程序会输出：
- 配置文件加载状态
- 读取的配置值
- 参数验证结果
- 各配置项的详细内容

---

## 构建所有测试程序

### 使用 CMake 构建

```bash
cd build
cmake ..
make
```

### 运行所有测试

```bash
cd build/bin
./test_particle    # 运行粒子类功能测试
./test_config      # 运行配置文件读取测试
```

---

## 测试数据准备

确保以下文件存在于正确的位置：

### 粒子数据文件（用于 test_particle）
- `data/fluid_particles.txt` / `data/fluid_particles.csv`
- `data/solid_particles.txt` / `data/solid_particles.csv`

### 配置文件（用于 test_config）
- `config.ini`

这些文件会在构建时自动复制到 `build/bin/` 目录下。

---

## 注意事项

1. **文件路径**：测试程序在 `build/bin/` 目录下运行，数据文件和配置文件会自动复制到该目录
2. **文件格式**：
   - 粒子数据文件：每行一个粒子，数据按位置、速度顺序排列
   - SolidParticle 数据还需包含法向向量
   - 支持空行和以 `#` 开头的注释行
3. **错误处理**：如果文件不存在或格式错误，测试程序会输出警告信息并使用默认值
