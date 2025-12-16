# 显式力计算模块

## 概述

本模块实现了MPS方法中的显式力计算功能，包括粘性力和重力的计算，以及基于这些力对粒子速度和位置的更新。

## 功能特性

- **速度拉普拉斯算子计算**：使用LSMPS corrective matrix方法计算速度场的拉普拉斯算子
- **粘性力计算**：根据速度拉普拉斯算子和动力学粘性系数计算粘性力加速度
- **重力计算**：直接使用重力加速度
- **速度和位置更新**：使用显式时间积分方法更新粒子的速度和位置
- **批量处理**：支持批量计算和更新所有粒子的速度和位置
- **边界条件**：统一使用第一类边界条件（无滑移边界，壁面速度为零）

## 技术细节

### 粘性力计算

粘性力加速度的计算公式为：

$$
\mathbf{a}_{viscous} = \nu \nabla^2 \mathbf{v}
$$

其中：
- $\nu$ 是动力学粘性系数
- $\nabla^2 \mathbf{v}$ 是速度场的拉普拉斯算子

速度拉普拉斯算子的计算使用LSMPS方法：

$$
\nabla^2 \mathbf{v} = \frac{2}{r_e} \sum_{j} w_{ij} \frac{\mathbf{v}_j - \mathbf{v}_i}{r_{ij}} [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}
$$

其中：
- $r_e$ 是平滑半径
- $w_{ij}$ 是权重函数
- $\mathbf{C}_3$ 和 $\mathbf{C}_4$ 是corrective matrix的第3行和第4行
- $\mathbf{P}$ 是基函数向量

### 重力计算

重力加速度直接使用输入值：

$$
\mathbf{a}_{gravity} = \mathbf{g}
$$

其中 $\mathbf{g} = (g_x, g_y)$ 是重力加速度向量。

### 速度和位置更新

使用显式时间积分方法（前向欧拉法）：

$$
\mathbf{v}^{n+1} = \mathbf{v}^n + \mathbf{a}^{total} \Delta t
$$

$$
\mathbf{x}^{n+1} = \mathbf{x}^n + \mathbf{v}^{n+1} \Delta t
$$

其中：
- $\mathbf{a}^{total} = \mathbf{a}_{viscous} + \mathbf{a}_{gravity}$ 是总加速度
- $\Delta t$ 是时间步长
- 上标 $n$ 表示当前时间步，$n+1$ 表示下一时间步

## 使用方法

### 基本用法

```cpp
#include "../src/explicit_force/ExplicitForce.hpp"
#include "../src/lsmps/CorrectiveMatrix.hpp"
#include "../src/core/Particle.hpp"

using namespace mps2D;

// 创建显式力计算器
ExplicitForce explicit_force;

// 计算corrective matrix（需要预先计算）
CorrectiveMatrix corrective_matrix_calculator;
std::vector<Eigen::Matrix<double, 5, 5>> matrices(fluid_particles.particle_num);
for (int i = 0; i < fluid_particles.particle_num; ++i) {
  matrices[i] = corrective_matrix_calculator.ComputeCorrectiveMatrix(
      i, fluid_particles, solid_particles, smoothing_radius, false);
}

// 批量计算并更新所有粒子的速度和位置
explicit_force.ComputeAndUpdateAllParticles(
    fluid_particles,
    solid_particles,
    matrices,
    smoothing_radius,
    kinematic_viscosity,  // 动力学粘性系数
    gravity_x,            // 重力加速度x分量
    gravity_y,            // 重力加速度y分量
    time_step);           // 时间步长
```

### 单独计算各个分量

如果需要单独计算各个分量，可以使用以下方法：

```cpp
// 计算速度拉普拉斯算子
double2 velocity_laplacian = explicit_force.ComputeVelocityLaplacian(
    particle_idx,
    fluid_particles.velocity,
    fluid_particles,
    solid_particles,
    corrective_matrix,
    smoothing_radius);

// 计算粘性力加速度
double2 viscous_acceleration = explicit_force.ComputeViscousAcceleration(
    velocity_laplacian,
    kinematic_viscosity);

// 计算重力加速度
double2 gravity_acceleration = explicit_force.ComputeGravityAcceleration(
    gravity_x,
    gravity_y);

// 更新速度和位置
explicit_force.UpdateVelocityAndPosition(
    particle_idx,
    fluid_particles,
    viscous_acceleration,
    gravity_acceleration,
    time_step);
```

### 边界条件

模块统一使用第一类边界条件（无滑移边界条件）：
- 壁面速度为零
- 固体邻域粒子使用标准基函数（与流体粒子相同）

## 依赖项

- **Eigen3**：用于矩阵运算（通过CMake自动下载或使用系统安装版本）
- **LSMPS CorrectiveMatrix模块**：用于计算corrective matrix
- **Particle模块**：用于粒子数据结构
- **MPSUtils**：用于工具函数（距离计算、权重函数等）

## 编译

模块已集成到CMake构建系统中，编译时会自动包含。

## 注意事项

1. **Corrective Matrix**：在调用显式力计算之前，必须预先计算所有粒子的corrective matrix。

2. **邻居列表**：在计算corrective matrix之前，必须确保已构建好邻居列表（`fluid_neighbour_list` 和 `solid_neighbour_list`）。

3. **时间步长**：显式时间积分方法对时间步长有稳定性要求。建议使用较小的时间步长以确保数值稳定性。

4. **边界条件**：模块统一使用第一类边界条件（无滑移边界，壁面速度为零），不需要设置边界条件参数。

5. **内存使用**：批量计算需要存储所有粒子的corrective matrix，内存占用为 $N \times 25 \times \text{sizeof(double)}$ 字节，其中 $N$ 是粒子数量。

6. **数值稳定性**：显式时间积分方法可能对时间步长敏感。如果遇到数值不稳定问题，建议减小时间步长或使用隐式时间积分方法。
