# LSMPS Corrective Matrix 模块

## 概述

本模块实现了LSMPS (Least Square Moving Particle Semi-implicit) 方法中的corrective matrix计算功能。Corrective matrix用于恢复MPS方法中离散化的一致性，提高计算精度。

## 功能特性

- **批量矩阵计算**：支持为所有粒子批量计算corrective matrix
- **并行计算**：使用OpenMP实现并行计算，提高性能
- **批量矩阵求逆**：针对5x5矩阵的批量求逆运算，使用Eigen库和并行优化
- **内存安全**：遵循Google C++编程规范，注重内存安全，防止内存泄漏

## 技术细节

### 矩阵大小

对于2D情况，corrective matrix是5x5的矩阵，对应5个基函数：
- x/r_e（归一化的x坐标）
- y/r_e（归一化的y坐标）
- (x/r_e)²（x的二次项）
- (y/r_e)²（y的二次项）
- (x*y / r_e^2) (混合二次项)

其中r_e是平滑半径（smoothing radius）。

### 计算方法

Corrective matrix的计算基于最小二乘方法：

1. **计算基函数$P$**：
   - **对于流体邻域粒子**，基函数为：
     $$P=[\frac{x}{r_e}, \frac{y}{r_e}, (\frac{x}{r_e})^2, (\frac{y}{r_e})^2, (\frac{x}{r_e})*(\frac{y}{r_e})]$$
     其中 $x = x_j - x_i$, $y = y_j - y_i$ 是相对位置，$r_e$ 是平滑半径。
   
   - **对于固体邻域粒子（壁面粒子）**，基函数为：
     $$P=[\mathbf{n}_x, \mathbf{n}_y, \frac{2\mathbf{n}_x x_{ij}}{r_e}, \frac{2\mathbf{n}_y y_{ij}}{r_e}, \frac{\mathbf{n}_x x_{ij} + \mathbf{n}_y y_{ij}}{r_e}]$$
     其中 $\mathbf{n}_x, \mathbf{n}_y$ 是壁面法向量的x和y分量，$x_{ij} = x_j - x_i$, $y_{ij} = y_j - y_i$ 是相对位置。

2. **计算权重$w$**：使用MPS权重函数计算每个邻域粒子的权重
   $$w(r) = (1 - r/r_e)^2, \quad r < r_e$$

3. **构建系数矩阵**：
   $$C=\sum_{j} w_{ij} \cdot \mathbf{P}_j \times \mathbf{P}_j^T$$
   其中对流体粒子和固体粒子分别使用对应的基函数计算方法。

4. **求逆得到corrective matrix**：
   $$M = C^{-1}$$

计算物理量$\phi$的梯度，散度或拉普拉斯算子 
$d_{ij} = \frac{\phi_j - \phi_i}{r_{ij}}$ 
1、梯度，对于标量$\phi$  
$$
\nabla \phi = \sum_{j}  w_{ij}d_{ij} \begin{bmatrix} \mathbf{C_1} \\ \mathbf{C_2} \end{bmatrix} \mathbf{P}
$$
2、散度，对于向量$\mathbf{\phi}$  
$$
\nabla \cdot \phi = \sum_{j}  w_{ij}\mathbf{d}_{ij} \begin{bmatrix} \mathbf{C_1} \\ \mathbf{C_2} \end{bmatrix} \mathbf{P}
$$
3、拉普拉斯算子，对于标量$\phi$ 
$$
\Delta \phi = 2\sum_{j}  w_{ij}d_{ij} \frac{[\mathbf{C_3}+\mathbf{C_4}]\mathbf{P}}{r_e}
$$

### 并行化策略

- 使用OpenMP并行化粒子级别的计算
- Eigen库内部使用多线程进行矩阵运算
- 支持自定义线程数

## 使用方法

### 基本用法

```cpp
#include "../src/lsmps/CorrectiveMatrix.hpp"

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

**注意**：
- 对于流体邻域粒子，使用标准的基函数计算方法
- 对于固体邻域粒子（壁面粒子），自动使用壁面基函数计算方法，需要确保 `SolidParticle` 的 `normal_vector` 已正确设置

## 依赖项

- **Eigen3**：用于矩阵运算（通过CMake自动下载或使用系统安装版本）
- **OpenMP**：用于并行计算（可选，如果未找到则使用串行版本）

## 编译

模块已集成到CMake构建系统中，编译时会自动：

1. 查找或下载Eigen3库
2. 查找OpenMP（如果可用）
3. 编译corrective matrix模块

## 测试

运行测试程序：

```bash
cd 2D/build
make test_corrective_matrix
./bin/test_corrective_matrix
```

测试程序会：
- 读取粒子数据
- 构建邻居列表
- 计算corrective matrix（逐个粒子计算）
- 显示前几个粒子的corrective matrix
- 显示性能统计信息和失败原因分析

## 注意事项

1. **邻域粒子数要求**：至少需要5个邻域粒子才能计算corrective matrix。如果邻域不足，会返回单位矩阵。

2. **矩阵奇异性**：如果矩阵不可逆（行列式接近0），会返回单位矩阵作为备用方案。

3. **壁面法向量**：对于固体粒子（壁面粒子），必须确保 `SolidParticle` 的 `normal_vector` 已正确设置。法向量应该指向流体域（从壁面指向流体）。

4. **基函数选择**：
   - 流体邻域粒子：使用标准基函数（基于相对位置）
   - 固体邻域粒子：使用壁面基函数（基于法向量和相对位置）

5. **内存使用**：对于N个粒子，需要存储N个5x5矩阵，内存占用为 N * 25 * sizeof(double) 字节。

6. **平滑半径设置**：建议平滑半径设为粒子间距的2-3倍，以确保每个粒子有足够的邻域粒子（至少5个）。




