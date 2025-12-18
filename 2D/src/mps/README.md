# 原始MPS方法模块

## 概述

本模块实现了原始MPS (Moving Particle Semi-implicit) 方法中的算子离散化功能。主要用于处理自由面粒子和近自由面粒子，这些粒子的邻域粒子数较少，无法使用LSMPS方法进行离散化。

## 功能特性

- **梯度计算**：计算标量场的梯度（如压力梯度）
- **散度计算**：计算向量场的散度（如速度散度）
- **拉普拉斯算子计算**：计算标量场的拉普拉斯算子（如压力拉普拉斯）
- **内存安全**：遵循Google C++编程规范，注重内存安全，防止内存泄漏
- **边界条件支持**：支持固体边界条件（壁面粒子）

## 技术细节

### 原始MPS方法公式

对于二维（2D）情况，原始MPS方法的算子离散化公式如下：

#### 1. 梯度计算（标量场）

对于标量场 $\phi$，梯度计算为：

$$
\nabla \phi_i = \frac{d}{n_0} \sum_{j \neq i} \frac{\phi_j - \phi_i}{|\mathbf{r}_j - \mathbf{r}_i|^2} (\mathbf{r}_j - \mathbf{r}_i) w(|\mathbf{r}_j - \mathbf{r}_i|)
$$

其中：
- $d = 2$（2D情况的维度）
- $n_0$ 是参考粒子数密度
- $w(r)$ 是MPS权重函数：$w(r) = (1 - r/r_e)^2$，$r < r_e$
- $r_e$ 是平滑半径

#### 2. 散度计算（向量场）

对于向量场 $\mathbf{u}$，散度计算为：

$$
\nabla \cdot \mathbf{u}_i = \frac{d}{n_0} \sum_{j \neq i} \frac{\mathbf{u}_j - \mathbf{u}_i}{|\mathbf{r}_j - \mathbf{r}_i|^2} \cdot (\mathbf{r}_j - \mathbf{r}_i) w(|\mathbf{r}_j - \mathbf{r}_i|)
$$

#### 3. 拉普拉斯算子计算（标量场）

对于标量场 $\phi$，拉普拉斯算子计算为：

$$
\Delta \phi_i = \frac{2d}{n_0 \lambda} \sum_{j \neq i} (\phi_j - \phi_i) w(|\mathbf{r}_j - \mathbf{r}_i|)
$$

其中 $\lambda$ 是拉普拉斯归一化参数：

**解析公式**：
$$
\lambda = \frac{\int_{V} w(r)r^{2} \,d \nu}{\int_{V} w(r) \,d \nu} = \frac{\int_{0}^{2\pi} \int_{0}^{r_e} (1 - \frac{\rho}{r_e})^{2}\rho^{3} \,d\rho \,d\phi}{\int_{0}^{2\pi} \int_{0}^{r_e} (1 - \frac{\rho}{r_e})^{2}\rho \,d\rho \,d\phi}=\frac{1}{5}r_e^{2}
$$

**离散方法**（实际使用）：
$$
\lambda = \frac{\sum_{j \neq i} w_{ij} |\mathbf{r}_j - \mathbf{r}_i|^2}{\sum_{j \neq i} w_{ij}}
$$

离散方法通过生成均匀分布的粒子网格来计算，使用与算例相同的粒子间距和平滑半径，确保对所有粒子都相同。

### 边界条件处理

- **流体邻域粒子**：使用标准MPS公式
- **固体邻域粒子（壁面）**：
  - 梯度计算：假设边界处标量值为0（Dirichlet边界条件）
  - 散度计算：假设壁面处速度为零（无滑移边界条件）
  - 拉普拉斯计算：假设边界处标量值为0

**注意**：边界条件可以根据实际需求进行修改。

## 使用方法

### 基本用法

```cpp
#include "../src/mps/OriginalMPS.hpp"
#include "../src/core/Particle.hpp"
#include "../src/surface_detection/SurfaceDetector.hpp"

using namespace mps2D;

// 创建原始MPS计算器
OriginalMPS mps_calculator;

// 计算参考粒子数密度（使用OriginalMPS的方法）
double reference_density = mps_calculator.ComputeReferenceDensity(
    fluid_particles, solid_particles, smoothing_radius);

// 计算标量场的梯度
double2 gradient = mps_calculator.ComputeGradient(
    particle_idx,
    pressure_field,  // 压力场
    fluid_particles,
    solid_particles,
    smoothing_radius,
    reference_density);

// 计算向量场的散度
double divergence = mps_calculator.ComputeDivergence(
    particle_idx,
    velocity_field,  // 速度场
    fluid_particles,
    solid_particles,
    smoothing_radius,
    reference_density);

// 计算标量场的拉普拉斯算子（使用离散lambda）
double laplacian = mps_calculator.ComputeLaplacian(
    particle_idx,
    pressure_field,  // 压力场
    fluid_particles,
    solid_particles,
    smoothing_radius,
    reference_density,
    particle_spacing);  // 粒子间距（用于计算离散lambda）
```

### 针对自由面和近自由面粒子的使用

```cpp
// 首先需要识别自由面和近自由面粒子
detector.DetectSurfaceParticles(
    fluid_particles, solid_particles, smoothing_radius, particle_spacing);

// 对自由面和近自由面粒子使用原始MPS方法
for (int i = 0; i < fluid_particles.particle_num; ++i) {
  SurfaceType surface_type = fluid_particles.surface_type[i];
  
  if (surface_type == SurfaceType::SURFACE || 
      surface_type == SurfaceType::NEAR_SURFACE) {
    // 使用原始MPS方法
    double2 grad = mps_calculator.ComputeGradient(
        i, pressure_field, fluid_particles, solid_particles,
        smoothing_radius, reference_density);
    
    double laplacian = mps_calculator.ComputeLaplacian(
        i, pressure_field, fluid_particles, solid_particles,
        smoothing_radius, reference_density, particle_spacing);
  } else {
    // 对内部粒子使用LSMPS方法
    // ...（使用CorrectiveMatrix模块）
  }
}
```

## 参数说明

### 输入参数

- **particle_idx**：粒子索引
- **scalar_field / vector_field**：标量场或向量场值
- **fluid_particles**：流体粒子对象（必须已构建邻域列表）
- **solid_particles**：固体粒子对象
- **smoothing_radius**：平滑半径（$r_e$），建议设为粒子间距的2-3倍
- **reference_density**：参考粒子数密度（$n_0$），可以使用`OriginalMPS::ComputeReferenceDensity`计算

### 输出

- **梯度**：返回`double2`类型，包含x和y方向的梯度分量
- **散度**：返回`double`类型标量值
- **拉普拉斯算子**：返回`double`类型标量值

## 注意事项

1. **邻域列表**：在调用本模块的函数之前，必须确保已构建好邻居列表（`fluid_neighbour_list` 和 `solid_neighbour_list`）。

2. **参考粒子数密度**：参考粒子数密度 $n_0$ 应该使用所有粒子的平均粒子数密度。可以使用`OriginalMPS::ComputeReferenceDensity`方法计算。

3. **平滑半径**：建议平滑半径设为粒子间距的2-3倍，以确保每个粒子有足够的邻域粒子。

4. **边界条件**：当前实现假设固体边界处的标量值为0，壁面处速度为零。如果实际边界条件不同，需要修改代码中的边界值。

5. **数值稳定性**：
   - 自动跳过距离小于 $10^{-10}$ 或大于平滑半径的粒子，避免数值问题
   - 当lambda参数过小时，拉普拉斯算子返回0

6. **适用范围**：本模块主要用于自由面粒子和近自由面粒子，这些粒子的邻域粒子数较少。对于内部粒子，建议使用LSMPS方法（`CorrectiveMatrix`模块）以获得更高的精度。

## 与LSMPS方法的对比

| 特性 | 原始MPS方法 | LSMPS方法 |
|------|------------|-----------|
| 适用范围 | 自由面、近自由面粒子 | 内部粒子 |
| 邻域粒子数要求 | 较少（≥3个） | 较多（≥5个） |
| 计算精度 | 较低 | 较高 |
| 计算复杂度 | 较低 | 较高（需要矩阵求逆） |
| 适用场景 | 邻域粒子不足时 | 邻域粒子充足时 |

## 依赖项

- **core/Particle.hpp**：粒子数据结构
- **core/MPSUtils.h**：MPS工具函数（距离计算、权重函数等）
- **surface_detection/SurfaceDetector.hpp**：自由面检测（可选，用于识别自由面粒子）

## 编译

模块已集成到CMake构建系统中，编译时会自动包含本模块。

## 示例代码

完整的使用示例请参考测试程序（如果存在）。


# 存在的问题
## 梯度计算
1、原始方法对于自由面以及近壁面的压力梯度的计算存在缺陷。对于近壁面粒子，周围的壁面必须具有准确的压力才能准确算出压力梯度；对于自由面粒子，将周围的气体粒子视为和粒子具有相同的压力，这样导致压力梯度累加过程中只累加了流体粒子对其的加权项，但是却同样使用了$n_0$对其归一化，这样就导致分子不变，分母减小，实际计算出的压力梯度小于理论值。例如静水问题中远离壁面的自由面粒子，对于该类粒子，有一半的区域是空的（气体粒子填充），这样实际的pnd应该接近n_0的一半，这就导致算出的压力梯度只有理论值-9800的一半。或许将$n_0$改为$n^*$会有缓解。
## 拉普拉斯算子计算
对于拉普拉斯算子计算，几乎必须要在各个方向具有两层邻域粒子才能计算出准确的结果，这也是原始mps方法的限制所在。对于自由面使用拉普拉斯算子必须补偿气体粒子的影响。但是对近壁面的自由面粒子，如何补偿是需要进一步考虑的
