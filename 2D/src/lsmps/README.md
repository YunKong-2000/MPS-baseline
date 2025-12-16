# LSMPS Corrective Matrix 模块

## 概述

本模块实现了LSMPS (Least Square Moving Particle Semi-implicit) 方法中的corrective matrix计算功能。Corrective matrix用于恢复MPS方法中离散化的一致性，提高计算精度。

## 功能特性

- **批量矩阵计算**：支持为所有粒子批量计算corrective matrix
- **并行计算**：使用OpenMP实现并行计算，提高性能
- **批量矩阵求逆**：针对5x5矩阵的批量求逆运算，使用Eigen库和并行优化
- **内存安全**：遵循Google C++编程规范，注重内存安全，防止内存泄漏
- **边界条件支持**：支持第一类（Dirichlet）和第二类（Neumann）边界条件

## 技术细节

### 矩阵大小

对于2D情况，corrective matrix是5×5的矩阵，对应5个基函数：

- $x/r_e$（归一化的x坐标）
- $y/r_e$（归一化的y坐标）
- $(x/r_e)^2$（x的二次项）
- $(y/r_e)^2$（y的二次项）
- $(x \cdot y) / r_e^2$（混合二次项）

其中 $r_e$ 是平滑半径（smoothing radius）。

### 计算方法

Corrective matrix的计算基于最小二乘方法：

#### 1. 计算基函数 $\mathbf{P}$

**对于流体邻域粒子**，基函数为：

$$
\mathbf{P} = \left[\frac{x}{r_e}, \frac{y}{r_e}, \left(\frac{x}{r_e}\right)^2, \left(\frac{y}{r_e}\right)^2, \frac{x \cdot y}{r_e^2}\right]
$$

其中 $x = x_j - x_i$，$y = y_j - y_i$ 是相对位置，$r_e$ 是平滑半径。

**对于固体邻域粒子（壁面粒子）**，基函数为：

$$
\mathbf{P} = \left[\mathbf{n}_x, \mathbf{n}_y, \frac{2\mathbf{n}_x x_{ij}}{r_e}, \frac{2\mathbf{n}_y y_{ij}}{r_e}, \frac{\mathbf{n}_x y_{ij} + \mathbf{n}_y x_{ij}}{r_e}\right]
$$

其中 $\mathbf{n}_x, \mathbf{n}_y$ 是壁面法向量的x和y分量，$x_{ij} = x_j - x_i$，$y_{ij} = y_j - y_i$ 是相对位置。

**注意**：当使用第一类边界条件（`border_condition = false`）时，固体邻域粒子也使用标准基函数（与流体粒子相同）。

#### 2. 计算权重 $w$

使用MPS权重函数计算每个邻域粒子的权重：

$$
w(r) = \left(1 - \frac{r}{r_e}\right)^2, \quad r < r_e
$$

#### 3. 构建系数矩阵

$$
\mathbf{C} = \sum_{j} w_{ij} \cdot \mathbf{P}_j \mathbf{P}_j^T
$$

其中对流体粒子和固体粒子分别使用对应的基函数计算方法。

#### 4. 求逆得到corrective matrix

$$
\mathbf{M} = \mathbf{C}^{-1}
$$

### 物理量计算

使用corrective matrix可以计算物理量 $\phi$ 的梯度、散度或拉普拉斯算子。

首先定义差分项：

$$
d_{ij} = \frac{\phi_j - \phi_i}{r_{ij}}
$$

其中 $r_{ij}$ 是粒子 $i$ 和 $j$ 之间的距离。

#### 1. 梯度计算（标量场）

对于标量 $\phi$，梯度计算为：

$$
\nabla \phi = \sum_{j} w_{ij} d_{ij} \begin{bmatrix} \mathbf{C}_1 \\ \mathbf{C}_2 \end{bmatrix} \mathbf{P}
$$

其中 $\mathbf{C}_1$ 和 $\mathbf{C}_2$ 是corrective matrix的第一行和第二行。

#### 2. 散度计算（向量场）

对于向量场 $\boldsymbol{\phi}$，散度计算为：

$$
\nabla \cdot \boldsymbol{\phi} = \sum_{j} w_{ij} \mathbf{d}_{ij} \begin{bmatrix} \mathbf{C}_1 \\ \mathbf{C}_2 \end{bmatrix} \mathbf{P}
$$

其中 $\mathbf{d}_{ij} = [(v_{x,j} - v_{x,i})/r_{ij}, (v_{y,j} - v_{y,i})/r_{ij}]$ 是速度差分的向量形式。

#### 3. 拉普拉斯算子计算（标量场）

对于标量 $\phi$，拉普拉斯算子计算为：

$$
\Delta \phi = \frac{2}{r_e} \sum_{j} w_{ij} d_{ij} [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}
$$

其中 $\mathbf{C}_3$ 和 $\mathbf{C}_4$ 是corrective matrix的第三行和第四行。

**注意**：对于向量场的拉普拉斯算子，需要对每个分量分别计算。

### 并行化策略

- 使用OpenMP并行化粒子级别的计算
- Eigen库内部使用多线程进行矩阵运算
- 支持自定义线程数

## 使用方法

### 基本用法

```cpp
#include "../src/lsmps/CorrectiveMatrix.hpp"
#include "../src/core/Particle.hpp"

using namespace mps2D;

// 创建corrective matrix计算器
CorrectiveMatrix calculator;

// 为单个粒子计算corrective matrix
Eigen::Matrix<double, 5, 5> matrix = calculator.ComputeCorrectiveMatrix(
    particle_idx, fluid_particles, solid_particles, smoothing_radius);

// 如果需要计算所有粒子，可以循环调用
std::vector<Eigen::Matrix<double, 5, 5>> matrices(num_particles);
for (int i = 0; i < num_particles; ++i) {
    matrices[i] = calculator.ComputeCorrectiveMatrix(
        i, fluid_particles, solid_particles, smoothing_radius);
}
```

### 边界条件

`ComputeCorrectiveMatrix` 函数支持边界条件参数：

```cpp
// 第一类边界条件（Dirichlet边界条件，默认）
Eigen::Matrix<double, 5, 5> matrix = calculator.ComputeCorrectiveMatrix(
    particle_idx, fluid_particles, solid_particles, smoothing_radius, false);

// 第二类边界条件（Neumann边界条件）
Eigen::Matrix<double, 5, 5> matrix = calculator.ComputeCorrectiveMatrix(
    particle_idx, fluid_particles, solid_particles, smoothing_radius, true);
```

**边界条件说明**：

- **第一类边界条件** (`border_condition = false`)：
  - 适用于指定场变量值的边界（如速度为零的无滑移边界）
  - 固体邻域粒子使用标准基函数（与流体粒子相同）
  
- **第二类边界条件** (`border_condition = true`)：
  - 适用于指定场变量法向导数的边界（如压力边界）
  - 固体邻域粒子使用壁面基函数（基于法向量）

### 计算梯度示例

```cpp
// 计算标量场的梯度
double2 ComputeGradient(
    int particle_idx,
    const std::vector<double>& scalar_field,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const Eigen::Matrix<double, 5, 5>& corrective_matrix,
    double smoothing_radius) {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  double phi_i = scalar_field[particle_idx];
  
  // 提取corrective matrix的前两行
  Eigen::Matrix<double, 1, 5> C1 = corrective_matrix.row(0);  // x方向
  Eigen::Matrix<double, 1, 5> C2 = corrective_matrix.row(1);  // y方向
  
  double grad_x = 0.0;
  double grad_y = 0.0;
  
  // 处理流体邻域粒子
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = std::sqrt(dx * dx + dy * dy);
    
    if (dist < 1e-10 || dist > smoothing_radius) continue;
    
    double phi_j = scalar_field[j];
    double d_ij = (phi_j - phi_i) / dist;
    double weight = WeightFunction(dist, smoothing_radius);
    
    // 计算基函数
    Eigen::Matrix<double, 5, 1> P;
    P << dx / dist,
         dy / dist,
         dx * dx / (dist * smoothing_radius),
         dy * dy / (dist * smoothing_radius),
         dx * dy / (dist * smoothing_radius);
    
    // 计算梯度贡献
    grad_x += weight * d_ij * (C1 * P)(0, 0);
    grad_y += weight * d_ij * (C2 * P)(0, 0);
  }
  
  // 处理固体邻域粒子（根据边界条件选择基函数）
  // ...（省略固体粒子处理代码）
  
  return {grad_x, grad_y};
}
```

### 计算拉普拉斯算子示例

```cpp
// 计算标量场的拉普拉斯算子
double ComputeLaplacian(
    int particle_idx,
    const std::vector<double>& scalar_field,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const Eigen::Matrix<double, 5, 5>& corrective_matrix,
    double smoothing_radius) {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  double phi_i = scalar_field[particle_idx];
  
  // 提取corrective matrix的第3行和第4行
  Eigen::Matrix<double, 1, 5> C3 = corrective_matrix.row(2);  // x²项
  Eigen::Matrix<double, 1, 5> C4 = corrective_matrix.row(3);  // y²项
  
  double laplacian = 0.0;
  double scalar_factor = 2.0 / smoothing_radius;
  
  // 处理流体邻域粒子
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = std::sqrt(dx * dx + dy * dy);
    
    if (dist < 1e-10 || dist > smoothing_radius) continue;
    
    double phi_j = scalar_field[j];
    double d_ij = (phi_j - phi_i) / dist;
    double weight = WeightFunction(dist, smoothing_radius);
    
    // 计算基函数
    Eigen::Matrix<double, 5, 1> P;
    P << dx / dist,
         dy / dist,
         dx * dx / (dist * smoothing_radius),
         dy * dy / (dist * smoothing_radius),
         dx * dy / (dist * smoothing_radius);
    
    // 计算拉普拉斯算子贡献
    double C3P = (C3 * P)(0, 0);
    double C4P = (C4 * P)(0, 0);
    laplacian += scalar_factor * weight * d_ij * (C3P + C4P);
  }
  
  // 处理固体邻域粒子（根据边界条件选择基函数）
  // ...（省略固体粒子处理代码）
  
  return laplacian;
}
```

**注意事项**：

- 对于流体邻域粒子，使用标准的基函数计算方法
- 对于固体邻域粒子（壁面粒子），根据边界条件选择基函数：
  - 第一类边界条件：使用标准基函数
  - 第二类边界条件：使用壁面基函数，需要确保 `SolidParticle` 的 `normal_vector` 已正确设置

## 依赖项

- **Eigen3**：用于矩阵运算（通过CMake自动下载或使用系统安装版本）
- **OpenMP**：用于并行计算（可选，如果未找到则使用串行版本）

## 编译

模块已集成到CMake构建系统中，编译时会自动：

1. 查找或下载Eigen3库
2. 查找OpenMP（如果可用）
3. 编译corrective matrix模块

## 测试

### 测试程序

运行测试程序：

```bash
cd build/bin
./test_corrective_matrix
```

测试程序会：
- 读取粒子数据
- 构建邻居列表
- 计算corrective matrix（逐个粒子计算）
- 显示前几个粒子的corrective matrix
- 显示性能统计信息和失败原因分析

### 应用测试

本模块还包含以下应用测试程序，展示了corrective matrix在实际问题中的应用：

1. **test_hydrostatic_manual**：手动实现的LSMPS压力梯度计算参考实现
2. **test_hydrostatic_pressure**：使用CorrectiveMatrix接口的静水压力梯度测试
3. **test_pipe_flow**：管道流动测试，计算速度梯度、散度和拉普拉斯算子

详细说明请参考 [测试程序文档](../test/README.md)。

## 注意事项

1. **邻域粒子数要求**：至少需要5个邻域粒子才能计算corrective matrix。如果邻域不足，会返回单位矩阵。

2. **矩阵奇异性**：如果矩阵不可逆（行列式接近0），会返回单位矩阵作为备用方案。

3. **壁面法向量**：
   - 当使用第二类边界条件（`border_condition = true`）时，必须确保 `SolidParticle` 的 `normal_vector` 已正确设置
   - 法向量应该指向流体域（从壁面指向流体）
   - 法向量应该归一化（模长为1）

4. **基函数选择**：
   - **流体邻域粒子**：使用标准基函数（基于相对位置）
   - **固体邻域粒子**：
     - 第一类边界条件（`border_condition = false`）：使用标准基函数
     - 第二类边界条件（`border_condition = true`）：使用壁面基函数（基于法向量和相对位置）

5. **内存使用**：对于N个粒子，需要存储N个5×5矩阵，内存占用为 $N \times 25 \times \text{sizeof(double)}$ 字节。

6. **平滑半径设置**：建议平滑半径设为粒子间距的2-3倍，以确保每个粒子有足够的邻域粒子（至少5个）。

7. **邻居列表**：在计算corrective matrix之前，必须确保已构建好邻居列表（`fluid_neighbour_list` 和 `solid_neighbour_list`）。

8. **距离检查**：在计算基函数时，会自动跳过距离小于 $10^{-10}$ 或大于平滑半径的粒子，避免数值问题。




