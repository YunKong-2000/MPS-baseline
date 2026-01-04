# Correction（压力修正）模块

## 概述

Correction模块位于压力求解（PPE）之后，根据PPE求解出的压力计算压力梯度，并根据压力梯度计算粒子的加速度，从而进一步矫正粒子的速度和位置。

## 功能特性

- **压力梯度计算**：使用LSMPS方法计算压力梯度，支持第二类边界条件（Neumann边界条件）
- **加速度计算**：根据压力梯度计算粒子加速度，公式为 `a = -(∇p / ρ)`
- **速度和位置更新**：根据加速度和时间步长更新粒子的速度和位置
- **批量处理**：支持批量计算和更新所有粒子的速度和位置

## 算法步骤

1. **计算压力梯度**：
   - 所有流体粒子使用LSMPS方法来离散压力梯度算子计算压力梯度
   - 壁面边界条件为第二类边界条件：`dp/dn = ρ * n · g`
   - 其中 `n` 是壁面法向量，`g` 是重力加速度向量，`ρ` 是流体密度

2. **计算粒子加速度**：
   - 加速度 = -(压力梯度/密度)
   - 公式：`a = -∇p / ρ`

3. **更新速度和位置**：
   - 速度更新：`v_new = v_old + a * dt`
   - 位置更新：`x_new = x_old + v_new * dt`

## 使用方法

### 基本用法

```cpp
#include "../src/correction/Correction.hpp"
#include "../src/lsmps/CorrectiveMatrix.hpp"
#include "../src/core/Particle.hpp"

using namespace mps2D;

// 创建correction模块实例
Correction correction;

// 假设已经完成PPE求解，压力已存储在fluid_particles.pressure中
// 假设已经计算了corrective matrices（使用第二类边界条件）
std::vector<Eigen::Matrix<double, 5, 5>> corrective_matrices_pressure;

// 批量计算并更新所有粒子的速度和位置
correction.ComputeAndUpdateAllParticles(
    fluid_particles,           // 流体粒子（会被修改）
    solid_particles,            // 固体粒子
    corrective_matrices_pressure,  // corrective matrices（第二类边界条件）
    smoothing_radius,           // 平滑半径
    gravity_x,                   // 重力加速度x分量
    gravity_y,                   // 重力加速度y分量
    density,                     // 流体密度
    time_step                    // 时间步长
);
```

### 单独使用各个功能

```cpp
// 1. 计算单个粒子的压力梯度
double2 pressure_gradient = correction.ComputePressureGradient(
    particle_idx,
    fluid_particles,
    solid_particles,
    corrective_matrix,      // 单个粒子的corrective matrix
    smoothing_radius,
    gravity_x,
    gravity_y,
    density
);

// 2. 计算加速度
double2 acceleration = correction.ComputeAcceleration(
    pressure_gradient,
    density
);

// 3. 更新速度和位置
correction.UpdateVelocityAndPosition(
    particle_idx,
    fluid_particles,
    acceleration,
    time_step
);
```

## 完整示例

```cpp
#include "../src/correction/Correction.hpp"
#include "../src/lsmps/CorrectiveMatrix.hpp"
#include "../src/PPE/PPESolver.hpp"
#include "../src/PPE/PPEMatrixBuilder.hpp"

using namespace mps2D;

// 1. 完成PPE求解（假设已完成）
// ... PPE求解代码 ...

// 2. 计算corrective matrices（第二类边界条件，用于压力梯度计算）
CorrectiveMatrix corrective_matrix_calc;
std::vector<Eigen::Matrix<double, 5, 5>> corrective_matrices_pressure(
    fluid_particles.particle_num);
for (int i = 0; i < fluid_particles.particle_num; ++i) {
    corrective_matrices_pressure[i] = corrective_matrix_calc.ComputeCorrectiveMatrix(
        i, fluid_particles, solid_particles, smoothing_radius, true);  // true = 第二类边界条件
}

// 3. 使用correction模块更新速度和位置
Correction correction;
correction.ComputeAndUpdateAllParticles(
    fluid_particles,
    solid_particles,
    corrective_matrices_pressure,
    smoothing_radius,
    gravity_x,
    gravity_y,
    density,
    time_step
);
```

## 注意事项

1. **Corrective Matrix**：
   - 必须使用第二类边界条件（`border_condition = true`）计算corrective matrix
   - 确保corrective matrix的数量与流体粒子数一致

2. **邻域列表**：
   - 在使用本模块之前，需要确保流体粒子的`fluid_neighbour_list`和`solid_neighbour_list`已经正确构建

3. **压力值**：
   - 确保PPE求解已完成，压力值已存储在`fluid_particles.pressure`中

4. **壁面边界条件**：
   - 壁面处压力梯度的法向分量：`dp/dn = ρ * n · g`
   - 确保`SolidParticle`的`normal_vector`已正确设置

5. **时间步长**：
   - 时间步长应该满足CFL条件，以确保数值稳定性

6. **密度**：
   - 确保密度值大于0，否则会导致数值不稳定

## 技术细节

### 压力梯度计算

对于流体邻域粒子，使用标准基函数：
```
∇p = Σ w_ij * (p_j - p_i) / r_ij * [C1; C2] * P
```

对于壁面邻域粒子，使用壁面基函数（第二类边界条件）：
```
∇p = Σ w_ij * (ρ * n · g) * [C1; C2] * P_wall
```

其中：
- `w_ij` 是权重函数
- `C1, C2` 是corrective matrix的第一行和第二行
- `P` 是基函数向量
- `n` 是壁面法向量
- `g` 是重力加速度向量

### 加速度计算

```
a = -∇p / ρ
```

### 速度和位置更新

使用显式时间积分：
```
v_new = v_old + a * dt
x_new = x_old + v_new * dt
```

## 与其他模块的集成

Correction模块通常与以下模块配合使用：

1. **PPE模块**：在PPE求解完成后调用
2. **LSMPS模块**：用于计算corrective matrix
3. **NeighborList模块**：用于构建邻域列表
4. **SurfaceDetection模块**：用于检测自由面粒子（如果需要）

## 内存安全

- 遵循Google C++编程规范
- 使用智能指针和RAII原则管理资源
- 避免内存泄漏和悬空指针
- 所有数组访问都进行边界检查

## 可扩展性

- 代码设计考虑了未来扩展到CPU/GPU并行架构的可能性
- 粒子级别的计算可以轻松并行化
- 支持OpenMP并行化（如果启用）

