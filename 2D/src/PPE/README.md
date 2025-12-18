# PPE（压力泊松方程）模块

## 概述

本模块实现了压力泊松方程（Pressure Poisson Equation, PPE）的系数矩阵构建和迭代求解功能。基于LSMPS算子离散方法构建PPE系数矩阵，并使用迭代方法求解大规模非对称稀疏线性方程组。

## 功能特性

- **系数矩阵构建**：基于LSMPS算子离散方法构建PPE系数矩阵A和右边项b
- **压缩稀疏矩阵存储**：使用Eigen::SparseMatrix的CSC（Compressed Sparse Column）压缩格式，显著减少内存占用
- **高效构建**：使用三元组格式构建，按列排序后转换为压缩格式，优化内存分配和访问性能
- **迭代求解**：支持BiCGSTAB和GMRES两种迭代求解器，充分利用压缩格式的矩阵-向量乘法性能
- **内存优化**：预统计非零元素个数，精确预分配内存，避免动态扩容
- **内存安全**：遵循Google C++编程规范，注重内存安全，防止内存泄漏

## 技术细节

### PPE方程

压力泊松方程的形式为：

$$
\frac{1}{\rho} \langle \Delta p \rangle_i = \frac{1}{\Delta t} \nabla \cdot \mathbf{u}^k
$$

其中：
- $p$ 是压力
- $\rho$ 是流体密度
- $\Delta t$ 是时间步长
- $\mathbf{u}^k$ 是当前时间步的速度场

### 离散方法

根据 `ppeMatrix.md` 中的方法，对于内部粒子（目前所有流体粒子都按内部粒子处理），PPE方程的离散形式为：

$$
\left(-\frac{2}{r_e \rho} \sum_{j \neq i} \frac{w_{ij}}{r_{ij}} [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij}\right) p_i + \sum_{j \neq i} \left(\frac{2}{r_e \rho} \frac{w_{ij}}{r_{ij}} [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij}\right) p_j = \frac{1}{\Delta t} \sum_{j \neq i} w_{ij} \frac{\mathbf{u}_j - \mathbf{u}_i}{r_{ij}} \begin{bmatrix} \mathbf{C}_1 \\ \mathbf{C}_2 \end{bmatrix} \mathbf{P}_{ij}
$$

其中：
- $r_e$ 是平滑半径
- $w_{ij}$ 是权重函数
- $r_{ij}$ 是粒子$i$和$j$之间的距离
- $\mathbf{C}_1, \mathbf{C}_2, \mathbf{C}_3, \mathbf{C}_4$ 是corrective matrix的行向量
- $\mathbf{P}_{ij}$ 是基函数向量

### 系数矩阵结构

对于远离壁面的内部粒子：
- **对角线元素**：$A_{ii} = -\frac{2}{r_e \rho} \sum_{j \neq i} \frac{w_{ij}}{r_{ij}} [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij}$
- **非对角线元素**：$A_{ij} = \frac{2}{r_e \rho} \frac{w_{ij}}{r_{ij}} [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij}$（当粒子$j$是粒子$i$的流体邻域粒子时）
- **右边项**：$b_i = \frac{1}{\Delta t} \sum_{j \neq i} w_{ij} \frac{\mathbf{u}_j - \mathbf{u}_i}{r_{ij}} \begin{bmatrix} \mathbf{C}_1 \\ \mathbf{C}_2 \end{bmatrix} \mathbf{P}_{ij}$

对于近壁面内部粒子（有壁面邻域粒子时）：
- **对角线元素**：$A_{ii} = -\frac{2}{r_e \rho} \sum_{j \in \mathrm{fluid}} \frac{w_{ij}}{r_{ij}} [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij}$（只考虑流体邻域粒子）
- **非对角线元素**：$A_{ij} = \frac{2}{r_e \rho} \frac{w_{ij}}{r_{ij}} [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij}$（只考虑流体邻域粒子）
- **右边项**：$b_i = \frac{1}{\Delta t} \sum_{j \in \mathrm{fluid}} w_{ij} \frac{\mathbf{u}_j - \mathbf{u}_i}{r_{ij}} \begin{bmatrix} \mathbf{C}_1 \\ \mathbf{C}_2 \end{bmatrix} \mathbf{P}_{ij} + \frac{1}{\Delta t} \sum_{j \in \mathrm{wall}} w_{ij} \frac{\mathbf{u}_{\mathrm{wall}} - \mathbf{u}_i}{r_{ij}} \begin{bmatrix} \mathbf{C}_1 \\ \mathbf{C}_2 \end{bmatrix} \mathbf{P}_{ij} + \frac{2}{r_e \rho} \sum_{j \in \mathrm{wall}} w_{ij} (\rho \mathbf{n}_{\mathrm{wall}} \cdot \mathbf{g}) [\mathbf{C}_3 + \mathbf{C}_4] \mathbf{P}_{ij}$

### 求解器

支持两种迭代求解器：

1. **BiCGSTAB**：双共轭梯度稳定化方法，推荐用于非对称矩阵
2. **GMRES**：广义最小残差方法，适用于非对称矩阵

两种求解器都使用对角预处理器来改善收敛性。

### 稀疏矩阵压缩格式

本模块使用Eigen的CSC（Compressed Sparse Column）压缩格式存储稀疏矩阵：

- **存储格式**：`Eigen::SparseMatrix<double, Eigen::ColMajor>`
- **压缩优势**：
  - 显著减少内存占用（只存储非零元素）
  - 提高矩阵-向量乘法性能（缓存友好的访问模式）
  - 适合迭代求解器的矩阵操作

- **构建过程**：
  1. 使用三元组（Triplet）格式收集所有非零元素
  2. 按列优先排序三元组，提高压缩效率
  3. 调用 `setFromTriplets()` 构建稀疏矩阵
  4. 调用 `makeCompressed()` 转换为压缩格式

- **内存估算**：
  - 对于N个粒子，每个粒子平均有M个邻域粒子
  - 非零元素个数：约 N + N×M
  - 压缩格式内存：约 (N + N×M) × (sizeof(double) + sizeof(int)) + N × sizeof(int)

## 使用方法

### 基本用法

```cpp
#include "PPE/PPEMatrixBuilder.hpp"
#include "PPE/PPESolver.hpp"
#include "lsmps/CorrectiveMatrix.hpp"
#include "core/Particle.hpp"

using namespace mps2D;

// 1. 计算所有粒子的corrective matrix
CorrectiveMatrix corrective_matrix_calculator;
std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>> 
    corrective_matrices(fluid_particles.particle_num);
    
for (int i = 0; i < fluid_particles.particle_num; ++i) {
    corrective_matrices[i] = corrective_matrix_calculator.ComputeCorrectiveMatrix(
        i, fluid_particles, solid_particles, smoothing_radius, false);
}

// 2. 构建PPE系数矩阵和右边项
PPEMatrixBuilder matrix_builder;
// 使用ColMajor格式（CSC压缩格式），适合迭代求解器
Eigen::SparseMatrix<double, Eigen::ColMajor> A;
Eigen::VectorXd b;

// 重力加速度（用于壁面压力边界条件）
double gravity_x = 0.0;
double gravity_y = -9.8;

bool success = matrix_builder.BuildPPEMatrix(
    fluid_particles, solid_particles, corrective_matrices,
    smoothing_radius, density, time_step, gravity_x, gravity_y, A, b);

if (!success) {
    std::cerr << "构建PPE矩阵失败" << std::endl;
    return;
}

// 3. 求解PPE方程
PPESolver::SolverConfig solver_config;
solver_config.solver_type = PPESolver::SolverType::BICGSTAB;
solver_config.max_iterations = 1000;
solver_config.tolerance = 1e-6;

PPESolver solver(solver_config);
Eigen::VectorXd pressure;

success = solver.Solve(A, b, pressure);

if (success) {
    // 将求解结果写入粒子压力数组
    for (int i = 0; i < fluid_particles.particle_num; ++i) {
        fluid_particles.pressure[i] = pressure(i);
    }
    
    std::cout << "求解成功，迭代次数：" << solver.GetLastIterations() 
              << "，残差：" << solver.GetLastResidual() << std::endl;
} else {
    std::cerr << "求解失败" << std::endl;
}
```

### 配置求解器

```cpp
// 使用GMRES求解器
PPESolver::SolverConfig config;
config.solver_type = PPESolver::SolverType::GMRES;
config.max_iterations = 2000;
config.tolerance = 1e-8;
config.restart = 50;  // GMRES重启参数

PPESolver solver(config);
```

## 注意事项

1. **邻域列表**：在使用本模块之前，需要确保流体粒子的`fluid_neighbour_list`和`solid_neighbour_list`已经正确构建。

2. **Corrective Matrix**：需要预先计算所有粒子的corrective matrix，并确保数量与流体粒子数一致。

3. **内存使用**：对于大规模问题，稀疏矩阵的内存占用可能较大。建议在构建矩阵前使用`CountNonZeros()`方法估算内存需求。

4. **收敛性**：如果求解器未收敛，可以尝试：
   - 降低收敛容差（tolerance）
   - 增加最大迭代次数
   - 使用不同的求解器类型
   - 检查系数矩阵的条件数

5. **当前限制**：
   - 所有流体粒子都按内部粒子处理，不区分自由面粒子
   - 暂时不对方程右边项添加其他源项
   - 对于有壁面邻域粒子的流体粒子，按照近壁面内部粒子的公式处理，考虑了壁面速度散度项和壁面压力边界条件项

## 依赖项

- **Eigen3**：用于稀疏矩阵和迭代求解器（通过CMake自动下载或使用系统安装版本）
- **LSMPS模块**：需要corrective matrix计算结果

## 编译

模块已集成到CMake构建系统中，编译时会自动包含PPE模块的源文件。
