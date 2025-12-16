# 2D 模拟测试程序说明

本文档说明2D子项目中所有测试程序的用途和使用方法。

## 测试程序列表

2D子项目包含以下测试程序：

1. **test_types** - Types.h 中基础类型测试
2. **test_mps_utils** - MPSUtils.h 中工具函数测试
3. **test_particle** - Particle 类功能测试
4. **test_file_operator** - FileOperator 文件读写模块测试
5. **test_neighbor_list** - NeighborListSearcher 邻居列表搜索器测试
6. **test_surface_detection** - SurfaceDetector 自由面检测测试
7. **test_hydrostatic_manual** - LSMPS压力梯度计算参考实现（手动实现）
8. **test_hydrostatic_pressure** - 静水压力梯度测试（使用CorrectiveMatrix接口）
9. **test_pipe_flow** - 管道流动测试（速度梯度、散度、拉普拉斯算子）

---

## test_types - 基础类型测试

### 用途

测试 `Types.h` 中定义的基础类型，包括 `double2`、`int2` 结构体和 `SurfaceType` 枚举类。

### 运行方法

```bash
cd build/bin
./test_types
```

### 测试内容

#### 1. double2 结构体测试
- 基本构造和成员访问
- 向量运算（加法、减法、除法、点积）

#### 2. int2 结构体测试
- 基本构造和成员访问
- 向量运算（加法、减法、点积）

#### 3. SurfaceType 枚举测试
- 枚举值定义和访问
- 枚举向量使用

---

## test_mps_utils - 工具函数测试

### 用途

测试 `MPSUtils.h` 中提供的所有工具函数，包括向量运算、距离计算、权重函数、向量归一化等。

### 运行方法

```bash
cd build/bin
./test_mps_utils
```

### 测试内容

#### 1. double2 基本运算测试
- 向量加法、减法、除法
- 向量点积
- 辅助函数（MakeDouble2、AssignDouble2）
- 流输出操作符

#### 2. int2 基本运算测试
- 向量加法、减法
- 向量点积
- 辅助函数（MakeInt2、AssignInt2）
- 流输出操作符

#### 3. 距离计算测试
- 欧氏距离计算（二维）
- 零距离边界情况
- 距离计算正确性验证

#### 4. 权重函数测试
- MPS标准权重函数计算
- 边界条件验证（r >= r_e 时权重为0）
- 不同距离下的权重值

#### 5. 向量模长计算测试
- 向量模长计算
- 零向量模长
- 模长计算正确性验证

#### 6. 向量归一化测试
- 向量归一化功能
- 归一化后向量模长为1的验证
- 零向量归一化处理

#### 7. 点积计算测试
- 点积计算功能
- 正交向量点积验证（应为0）

---

## test_particle - 粒子类功能测试

### 用途

测试 `Particle` 类及其派生类（`FluidParticle`、`SolidParticle`）的基本功能，包括创建、数据操作、拷贝等。

### 运行方法

```bash
cd build/bin
./test_particle
```

### 测试内容

#### 1. Particle 类创建测试
- 基础 Particle 类初始化
- 成员变量默认值验证

#### 2. Particle 数据操作测试
- 手动设置粒子数据
- 位置和速度向量操作
- `getParticleNum()` 方法测试

#### 3. FluidParticle 类测试
- FluidParticle 创建和初始化
- 密度、压力、表面类型设置
- 邻居列表初始化

#### 4. SolidParticle 类测试
- SolidParticle 创建和初始化
- 法向向量设置和验证
- 法向向量模长检查

#### 5. 粒子数据拷贝测试
- FluidParticle 数据拷贝功能
- SolidParticle 数据拷贝功能
- 拷贝后数据一致性验证

#### 6. 粒子间距离计算测试
- 使用工具函数计算粒子间距离
- 距离计算正确性验证

### 注意事项

⚠️ **重要**: 此测试程序不包含文件I/O功能测试，因为 `FileOperator` 的2D版本尚未完全实现。当前 `FileOperator` 仅提供占位符实现，`getParticleFromFile()` 方法会返回 -1 表示功能未实现。测试主要关注粒子类的内存操作和基本功能。

---

## 构建测试程序

### 使用 CMake 构建

需要在主项目的 `CMakeLists.txt` 中添加2D子项目的测试程序配置，或创建独立的2D子项目CMake配置。

### 运行所有测试

```bash
cd build/bin
./test_types                    # 运行基础类型测试
./test_mps_utils                # 运行工具函数测试
./test_particle                 # 运行粒子类功能测试
./test_file_operator            # 运行文件读写模块测试
./test_neighbor_list            # 运行邻居列表搜索器测试
./test_surface_detection        # 运行自由面检测测试
./test_hydrostatic_manual       # 运行LSMPS压力梯度参考实现测试
./test_hydrostatic_pressure     # 运行静水压力梯度测试
./test_pipe_flow                # 运行管道流动测试
```

---

## 测试输出

所有测试程序都会输出：
- 详细的测试执行信息
- 测试结果验证（✓ 通过 或 ✗ 失败）
- 测试总结

---

## 与3D版本的差异

2D版本的测试程序与3D版本的主要差异：

1. **命名空间**: 使用 `mps2D` 命名空间而不是 `mps`
2. **向量类型**: 使用 `double2` 和 `int2` 而不是 `double3` 和 `int3`
3. **维度**: 所有位置、速度、法向向量只有x和y两个分量
4. **文件I/O**: 暂不包含文件读取/写入测试（FileOperator未实现）

---

## test_file_operator - 文件读写模块测试

### 用途

测试 `FileOperator` 类的所有文件读写功能，包括TXT/CSV文件读取、写入，以及VTK格式输出。包含一个从文件读取均匀分布的2D粒子并写入VTK格式的完整测试用例。

### 运行方法

```bash
cd build/bin
./test_file_operator
```

### 测试内容

#### 1. FluidParticle 读取功能测试
- 从TXT文件读取流体粒子数据
- 从CSV文件读取流体粒子数据
- 验证TXT和CSV读取的数据一致性

#### 2. SolidParticle 读取功能测试
- 从TXT文件读取固体粒子数据（包含法向向量）
- 验证法向向量数据

#### 3. 粒子数据写入功能测试
- FluidParticle 数据写入TXT文件
- FluidParticle 数据写入CSV文件
- SolidParticle 数据写入TXT文件（包含法向向量）
- SolidParticle 数据写入CSV文件（包含法向向量）

#### 4. 均匀分布粒子读取和VTK写入测试（主要测试用例）
- 从文件读取均匀分布的5x5网格粒子（25个粒子）
- 验证粒子分布和位置范围
- 设置测试数据（密度、压力、表面类型）
- 写入VTK基础文件（位置和速度）
- 追加密度标量到VTK文件
- 追加压力标量到VTK文件
- 追加表面类型标量到VTK文件

#### 5. 向量写入功能测试（Debug）
- double向量写入TXT文件
- double2向量写入TXT文件
- double2向量写入CSV文件

#### 6. SolidParticle VTK写入功能测试
- 写入SolidParticle基础VTK文件
- 追加法向向量到VTK文件

### 测试数据

测试使用的数据文件位于 `data/` 目录：
- `fluid_particles_2d.txt` / `fluid_particles_2d.csv`: 2D均匀分布的流体粒子数据（5x5网格，25个粒子）
  - 格式：位置x 位置y 速度x 速度y
- `solid_particles_2d.txt`: 2D均匀分布的固体粒子数据（边界粒子）
  - 格式：位置x 位置y 速度x 速度y 法向向量x 法向向量y

### 测试输出

测试结果文件会生成在 `build/bin/data/` 目录下：
- `output_fluid_2d.txt` / `output_fluid_2d.csv` - FluidParticle写入结果
- `output_solid_2d.txt` / `output_solid_2d.csv` - SolidParticle写入结果（包含法向向量）
- `uniform_fluid_2d.vtk` - 均匀分布粒子VTK文件（主要测试用例输出）
- `solid_particles_2d.vtk` - 固体粒子VTK文件
- `debug_double_2d.txt` / `debug_double2_2d.txt` / `debug_double2_2d.csv` - 调试输出

### VTK文件可视化

生成的VTK文件可以使用以下工具打开查看：
- **ParaView**: 推荐使用，功能强大
- **VisIt**: 科学可视化工具
- **VTK**: 其他支持VTK格式的可视化工具

在ParaView中打开VTK文件后，可以：
- 查看粒子位置分布
- 显示速度向量
- 查看密度、压力等标量场
- 查看表面类型分布

---

## test_neighbor_list - 邻居列表搜索器测试

### 用途

测试 `NeighborListSearcher` 类的所有核心功能，包括邻居列表构建、粒子排序、网格划分、距离计算等。该测试程序确保邻居列表搜索算法在2D场景下都能正确工作。

### 运行方法

```bash
cd build/bin
./test_neighbor_list
```

### 测试内容

#### 1. 基本邻居搜索功能测试
- 创建简单的3x3网格场景，每个网格一个粒子
- 验证邻居列表构建功能
- 验证流体-流体邻居搜索
- 检查邻居列表统计信息
- 验证邻居列表正确性（索引有效性、距离检查、自引用检查）

#### 2. 包含固体粒子的邻居搜索测试
- 测试包含流体粒子和固体粒子的场景
- 验证流体-流体和流体-固体邻居搜索功能
- 测试线性排列的粒子（5个粒子排成一条线）
- 验证边界固体粒子的邻居搜索

#### 3. 均匀分布粒子的邻居搜索测试
- 创建5x5网格的均匀分布粒子场景（25个粒子）
- 验证在均匀分布场景下的邻居搜索性能
- 检查最大/最小/平均邻居数量统计
- 验证中心粒子的邻居数量

#### 4. 空粒子列表测试
- 测试空流体粒子和空固体粒子的边界情况处理
- 验证程序不会因空列表而崩溃
- 确保空列表时返回正确的结果

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
- 邻居列表统计（总邻居数、最大邻居数、平均邻居数等）
- 验证结果（✓ 通过 或 ✗ 失败）
- 测试总结

### 2D版本特点

与3D版本的主要区别：
- **网格维度**：使用2D网格（x, y两个方向）而不是3D网格
- **相邻网格**：每个网格有8个相邻网格（3x3-1=8）而不是26个
- **域计算**：计算4个值的域（min_x, min_y, max_x, max_y）而不是6个值
- **向量类型**：使用 `double2` 和 `int2` 而不是 `double3` 和 `int3`

---

## test_surface_detection - 自由面检测测试

### 用途

测试 `SurfaceDetector` 类的所有核心功能，包括粒子数密度计算、粗筛、细筛（虚拟光源法）、近自由面判定等。该测试程序确保自由面检测算法在2D场景下都能正确工作。

### 运行方法

```bash
cd build/bin
./test_surface_detection
```

### 测试内容

#### 1. 基本自由面检测功能测试
- 创建简单的5x5网格场景
- 验证自由面检测算法
- 统计不同表面类型的粒子数量
- 显示边界粒子的表面类型

#### 2. 正方形壁面包围流体的自由面检测测试
- 创建20x20网格的流体粒子（400个粒子）
- 创建包围流体的正方形壁面（76个固体粒子）
- 构建邻居列表
- 检测自由面
- 计算阴影面积比例
- 输出VTK文件（包含表面类型和阴影面积比例标量）

#### 3. 100x100大尺寸正方形壁面的自由面检测测试
- 创建100x100网格的流体粒子（10,000个粒子）
- 创建包围流体的正方形壁面（396个固体粒子）
- 构建邻居列表
- 检测自由面
- 计算阴影面积比例
- 输出VTK文件（包含表面类型和阴影面积比例标量）

### 算法步骤

测试程序验证以下算法步骤：

1. **计算粒子数密度**：使用权重函数计算每个粒子的粒子数密度
2. **粗筛**：
   - 粒子数密度特别小或邻域粒子数很少的直接判定为飞溅粒子
   - 粒子数密度很大且接近参考粒子数密度的直接判定为内部粒子
   - 其余为待细筛粒子
3. **细筛**：使用虚拟光源法（2D圆形幕布）进行进一步判定
   - 阴影面积大的为内部粒子
   - 阴影面积小的为自由面粒子
4. **近自由面判定**：内部粒子中，1.5倍粒子间距内有自由面粒子的为近自由面粒子

### 验证内容

测试程序会自动验证以下内容：
- **表面类型分布**：统计不同表面类型的粒子数量
- **阴影面积比例**：计算并输出每个粒子的阴影面积比例
- **VTK输出**：验证VTK文件正确生成，包含表面类型和阴影面积比例标量

### 测试输出

测试结果文件会生成在 `build/bin/data/` 目录下：
- `surface_detection_20x20.vtk` - 20x20测试用例的VTK文件
- `surface_detection_100x100.vtk` - 100x100测试用例的VTK文件

每个VTK文件包含：
- 粒子位置和速度
- `surface_type` 标量：表面类型（0=INNER, 1=NEAR_SURFACE, 2=SURFACE, 3=SPLASH）
- `shadow_area_ratio` 标量：阴影面积比例（0.0-1.0）

### VTK文件可视化

在ParaView中打开VTK文件后，可以：
- 使用 `surface_type` 标量着色，查看不同表面类型的分布
- 使用 `shadow_area_ratio` 标量着色，查看阴影面积比例分布
- 分析边界粒子的表面类型是否正确

### 2D版本特点

与3D版本的主要区别：
- **圆形幕布**：使用圆形幕布而不是球面幕布
- **角度网格**：只需要角度theta (0到2π)，不需要theta和phi两个角度
- **阴影计算**：计算圆形幕布上的阴影弧长比例而不是球面阴影面积比例
- **向量类型**：使用 `double2` 而不是 `double3`

---

## test_hydrostatic_manual - LSMPS压力梯度计算参考实现

### 用途

手动实现的LSMPS（Least Square Moving Particle Semi-implicit）压力梯度计算，作为参考实现。该测试程序直接实现了LSMPS算法的核心逻辑，包括corrective matrix的计算和压力梯度的求解，用于验证`CorrectiveMatrix`类的实现是否正确。

### 运行方法

```bash
cd build/bin
./test_hydrostatic_manual
```

### 测试内容

#### 1. 静水压力场设置
- 创建矩形通道中的流体粒子（均匀分布）
- 创建底部和顶部壁面粒子
- 设置静水压力分布：`p(y) = ρ * g * (h - y)`，其中h为通道高度

#### 2. LSMPS算法实现
- **基函数计算**：
  - 流体粒子基函数：`[x/r, y/r, x²/(r·r_e), y²/(r·r_e), x·y/(r·r_e)]`
  - 壁面粒子基函数：`[n_x, n_y, 2·n_x·x/r_e, 2·n_y·y/r_e, (n_x·y + n_y·x)/r_e]`
- **权重函数**：`w(r) = (1 - r/r_e)²`
- **Corrective Matrix构建**：通过加权基函数外积构建系数矩阵
- **压力梯度计算**：使用corrective matrix计算压力梯度

#### 3. 结果验证
- 计算所有流体粒子的压力梯度
- 与理论值对比：理论梯度为 `∇p = (0, -ρ·g)`
- 输出误差统计

#### 4. VTK输出
- 输出所有流体粒子的压力梯度信息到`hydrostatic_pressure_gradient_manual.vtk`
- 包含计算值和理论值的对比

### 算法特点

- **手动实现**：不依赖`CorrectiveMatrix`类，直接实现LSMPS算法
- **参考实现**：作为其他测试程序的参考基准
- **完整流程**：包含从基函数计算到梯度求解的完整流程

### 测试输出

测试结果文件会生成在 `build/bin/` 目录下：
- `hydrostatic_pressure_gradient_manual.vtk` - 压力梯度VTK文件

VTK文件包含：
- 粒子位置和速度
- `computed_gradient` 向量：计算的压力梯度
- `theoretical_gradient` 向量：理论压力梯度 `(0, -ρ·g)`
- `gradient_error` 向量：梯度误差
- `y_coordinate` 标量：Y坐标（用于分析）

---

## test_hydrostatic_pressure - 静水压力梯度测试

### 用途

使用`CorrectiveMatrix`类接口测试LSMPS压力梯度计算。该测试程序验证`CorrectiveMatrix`类的实现是否正确，通过对比`test_hydrostatic_manual`的结果来验证接口的正确性。

### 运行方法

```bash
cd build/bin
./test_hydrostatic_pressure
```

### 测试内容

#### 1. 静水压力场设置
- 创建矩形通道中的流体粒子（均匀分布）
- 创建底部和顶部壁面粒子（带法向量）
- 设置静水压力分布：`p(y) = ρ * g * (h - y)`
- 构建邻居列表

#### 2. Corrective Matrix计算
- 使用`CorrectiveMatrix::ComputeCorrectiveMatrix`接口计算corrective matrix
- 支持边界条件参数（`border_condition`）
- 验证corrective matrix的有效性

#### 3. 压力梯度计算
- 使用corrective matrix的前两行（C1和C2）计算梯度
- 处理流体邻域粒子和固体邻域粒子
- 对于壁面粒子，使用壁面基函数和边界条件

#### 4. 压力拉普拉斯算子计算
- 使用corrective matrix的第3行和第4行（C3和C4）计算拉普拉斯算子
- 计算 `∇²p = ∂²p/∂x² + ∂²p/∂y²`
- 对于静水压力，理论拉普拉斯算子应为0

#### 5. 结果验证
- 与理论梯度对比：`∇p = (0, -ρ·g)`
- 与理论拉普拉斯算子对比：`∇²p = 0`
- 输出误差统计（平均误差、最大误差）

#### 6. VTK输出
- 输出压力、压力梯度、压力拉普拉斯算子到`hydrostatic_pressure.vtk`
- 包含计算值和理论值的对比

### 算法特点

- **接口调用**：使用`CorrectiveMatrix`类的标准接口
- **边界条件**：支持第一类和第二类边界条件
- **完整验证**：验证梯度计算和拉普拉斯算子计算

### 测试输出

测试结果文件会生成在 `build/bin/` 目录下：
- `hydrostatic_pressure.vtk` - 压力场VTK文件

VTK文件包含：
- 粒子位置和速度
- `pressure` 标量：压力值
- `computed_gradient` 向量：计算的压力梯度
- `theoretical_gradient` 向量：理论压力梯度
- `gradient_error` 向量：梯度误差
- `computed_laplacian` 标量：计算的拉普拉斯算子
- `theoretical_laplacian` 标量：理论拉普拉斯算子（0）
- `laplacian_error` 标量：拉普拉斯算子误差

### 与test_hydrostatic_manual的关系

- `test_hydrostatic_manual` 是手动实现的参考版本
- `test_hydrostatic_pressure` 使用标准接口，验证接口实现的正确性
- 两个测试的结果应该一致，用于验证`CorrectiveMatrix`类的正确性

---

## test_pipe_flow - 管道流动测试

### 用途

测试LSMPS方法在管道流动（Poiseuille流动）中的速度梯度、散度和拉普拉斯算子计算。该测试程序验证LSMPS方法在处理向量场（速度场）时的正确性，包括梯度张量、散度和拉普拉斯算子的计算。

### 运行方法

```bash
cd build/bin
./test_pipe_flow
```

### 测试内容

#### 1. Poiseuille流动设置
- 创建矩形通道中的流体粒子（均匀分布）
- 创建底部和顶部壁面粒子（无滑移边界条件，速度为零）
- 设置Poiseuille速度分布：`v_x(y) = v_max * (1 - (2y/h - 1)²)`，`v_y = 0`
- 构建邻居列表

#### 2. Corrective Matrix计算
- 使用`CorrectiveMatrix::ComputeCorrectiveMatrix`接口计算corrective matrix
- 使用第一类边界条件（Dirichlet边界条件，`border_condition = false`）
- 壁面速度为零（无滑移边界条件）

#### 3. 速度梯度计算
- 使用corrective matrix的前两行（C1和C2）计算速度梯度张量
- 计算四个分量：`∂v_x/∂x`, `∂v_x/∂y`, `∂v_y/∂x`, `∂v_y/∂y`
- 处理流体邻域粒子和固体邻域粒子
- 对于壁面粒子，使用标准基函数（第一类边界条件）

#### 4. 速度散度计算
- 使用速度梯度计算散度：`∇·v = ∂v_x/∂x + ∂v_y/∂y`
- 对于不可压缩流动，理论散度应为0

#### 5. 速度拉普拉斯算子计算
- 使用corrective matrix的第3行和第4行（C3和C4）计算拉普拉斯算子
- 计算 `∇²v = (∇²v_x, ∇²v_y)`
- 对于Poiseuille流动：
  - `∇²v_x = -8·v_max/h²`
  - `∇²v_y = 0`

#### 6. 结果验证
- 与理论梯度对比：
  - `∂v_x/∂x = 0`
  - `∂v_x/∂y = -4·v_max·(2y/h - 1)/h`
  - `∂v_y/∂x = 0`
  - `∂v_y/∂y = 0`
- 与理论散度对比：`∇·v = 0`
- 与理论拉普拉斯算子对比：`∇²v_x = -8·v_max/h²`, `∇²v_y = 0`
- 输出误差统计（平均误差、最大误差）

#### 7. VTK输出
- 输出速度、速度梯度、散度、拉普拉斯算子到`pipe_flow.vtk`
- 输出壁面粒子到`pipe_flow_wall_particles.vtk`
- 包含计算值和理论值的对比

### 算法特点

- **向量场处理**：处理速度向量场，计算梯度张量
- **边界条件**：使用第一类边界条件（Dirichlet），壁面速度为零
- **多物理量**：同时计算梯度、散度和拉普拉斯算子
- **完整验证**：验证所有计算的正确性

### 测试参数

测试中使用的关键参数：
- `channel_height`: 通道高度（默认0.5 m）
- `channel_length`: 通道长度（默认1.0 m）
- `particle_spacing`: 粒子间距（默认0.02 m）
- `v_max`: 最大速度（默认1.0 m/s）
- `smoothing_radius`: 平滑半径（默认2.1 * particle_spacing）

### 测试输出

测试结果文件会生成在 `build/bin/` 目录下：
- `pipe_flow.vtk` - 流体粒子VTK文件
- `pipe_flow_wall_particles.vtk` - 壁面粒子VTK文件

`pipe_flow.vtk`文件包含：
- 粒子位置和速度
- `theoretical_velocity` 向量：理论速度
- `computed_velocity` 向量：计算速度
- `computed_gradient` 张量：计算的速度梯度张量
- `theoretical_gradient` 张量：理论速度梯度张量
- `gradient_error` 张量：梯度误差
- `computed_divergence` 标量：计算的散度
- `theoretical_divergence` 标量：理论散度（0）
- `divergence_error` 标量：散度误差
- `computed_velocity_laplacian` 向量：计算的拉普拉斯算子
- `theoretical_velocity_laplacian` 向量：理论拉普拉斯算子
- `laplacian_error` 向量：拉普拉斯算子误差
- `laplacian_error_magnitude` 标量：拉普拉斯算子误差大小
- `y_coordinate` 标量：Y坐标（用于分析）

### VTK文件可视化

在ParaView中打开VTK文件后，可以：
- 查看速度向量分布
- 显示速度梯度张量的各个分量
- 查看散度分布（应该接近0）
- 查看拉普拉斯算子分布
- 分析边界附近的误差分布

### 理论背景

Poiseuille流动是管道中的层流流动，速度分布为：
- `v_x(y) = v_max * (1 - (2y/h - 1)²)`
- `v_y = 0`

其中：
- `v_max` 是管道中心的最大速度
- `h` 是通道高度
- `y` 是从底部壁面开始的距离

理论梯度：
- `∂v_x/∂x = 0`
- `∂v_x/∂y = -4·v_max·(2y/h - 1)/h`
- `∂v_y/∂x = 0`
- `∂v_y/∂y = 0`

理论拉普拉斯算子：
- `∇²v_x = ∂²v_x/∂x² + ∂²v_x/∂y² = 0 + (-8·v_max/h²) = -8·v_max/h²`
- `∇²v_y = 0`

---

## 未来扩展

计划添加的测试功能：
- 配置文件读取测试（2D版本）
- 其他流动场景的测试（如Couette流动、Stokes流动等）

