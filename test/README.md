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

## test_neighbor_list - 邻居列表搜索器功能测试

### 用途

测试 NeighborListSearcher 类的所有核心功能，包括邻居列表构建、粒子排序、网格划分、距离计算等。该测试程序确保邻居列表搜索算法在各种场景下都能正确工作。

### 运行方法

```bash
cd build/bin
./test_neighbor_list
```

### 测试内容

#### 1. 基本功能测试
- 创建简单的 3x3x3 网格场景，每个网格一个粒子
- 验证邻居列表构建功能
- 验证流体-流体和流体-固体邻居搜索
- 检查邻居列表统计信息

#### 2. 空粒子列表测试
- 测试空流体粒子和空固体粒子的边界情况处理
- 验证程序不会因空列表而崩溃
- 确保空列表时返回正确的结果

#### 3. 只有流体粒子测试
- 测试仅包含流体粒子的场景
- 验证流体-流体邻居搜索功能
- 测试线性排列的粒子（5个粒子排成一条线）

#### 4. 密集粒子场景测试
- 创建复杂的密集粒子场景（27个流体粒子 + 8个固体粒子）
- 验证在密集场景下的邻居搜索性能
- 检查最大邻居数量统计
- 验证算法在复杂场景下的正确性

#### 5. 粒子排序验证
- 验证粒子按网格索引排序功能
- 确保排序后粒子的所有属性（位置、速度、密度、压力等）一起移动
- 测试故意打乱顺序的粒子数据
- 验证排序后属性的一致性

#### 6. 边界粒子测试
- 测试粒子在计算域边界上的情况
- 验证分散分布的粒子（距离较远）
- 测试接近边界的粒子邻居搜索
- 验证边界情况下的距离计算

#### 7. 大网格尺寸测试
- 测试网格尺寸大于粒子间距的情况
- 验证大网格下的邻居搜索功能
- 确保算法在不同网格尺寸下都能正常工作

### 验证内容

测试程序会自动验证以下内容：
- **邻居列表正确性**：所有邻居的距离都在 `r_e`（搜索半径）范围内
- **索引有效性**：所有邻居索引都在有效范围内（0 到 particle_num-1）
- **自引用检查**：粒子不会将自己加入邻居列表
- **排序功能**：粒子排序后所有属性保持一致

### 测试参数

测试中使用的关键参数：
- `particle_radius`: 粒子初始间距
- `r_e`: 流体粒子相邻域半径（smoothing_radius）
- `r_cell`: 背景网格尺寸（cell_size）

### 测试输出

测试程序会输出：
- 每个测试用例的详细执行信息
- 邻居列表统计（总邻居数、最大邻居数等）
- 粒子排序前后的位置信息
- 验证结果（✓ 通过 或 ✗ 失败）
- 测试总结

### 测试数据

测试程序使用程序内部生成的测试数据，不需要外部数据文件。测试数据包括：
- 规则网格排列的粒子
- 线性排列的粒子
- 密集分布的粒子
- 边界和分散分布的粒子

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
./test_particle       # 运行粒子类功能测试
./test_config         # 运行配置文件读取测试
./test_neighbor_list  # 运行邻居列表搜索器功能测试
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
