#include "../src/PPE/PPEMatrixBuilder.hpp"
#include "../src/PPE/PPESolver.hpp"
#include "../src/lsmps/CorrectiveMatrix.hpp"
#include "../src/neighbour_list/NeighborListSearcher.hpp"
#include "../src/surface_detection/SurfaceDetector.hpp"
#include "../src/core/Particle.hpp"
#include "../src/core/FileOperator.hpp"
#include "core/Types.h"
#include "core/MPSUtils.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>

// PETSc头文件
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>

using namespace mps2D;

// 生成静水压力测试场景
// 参数：
//   container_width: 容器宽度
//   container_height: 容器高度
//   water_height: 水位高度
//   particle_spacing: 粒子间距
//   fluid_particles: 输出流体粒子
//   solid_particles: 输出固体粒子
void GenerateHydrostaticTest(
    double container_width,
    double container_height,
    double water_height,
    double particle_spacing,
    FluidParticle& fluid_particles,
    SolidParticle& solid_particles) {
  
  // 计算粒子数量
  int nx_fluid = static_cast<int>(container_width / particle_spacing) + 1;
  int ny_fluid = static_cast<int>(water_height / particle_spacing) + 1;
  int num_fluid = nx_fluid * ny_fluid;
  
  // 生成流体粒子（水）
  fluid_particles.particle_num = num_fluid;
  fluid_particles.position.resize(num_fluid);
  fluid_particles.velocity.resize(num_fluid);
  fluid_particles.pressure.resize(num_fluid);
  fluid_particles.density.resize(num_fluid);
  fluid_particles.surface_type.resize(num_fluid);
  fluid_particles.fluid_neighbour_list.resize(num_fluid);
  fluid_particles.solid_neighbour_list.resize(num_fluid);
  
  int idx = 0;
  for (int j = 0; j < ny_fluid; ++j) {
    for (int i = 0; i < nx_fluid; ++i) {
      double x = i * particle_spacing;
      double y = j * particle_spacing;
      fluid_particles.position[idx] = {x, y};
      fluid_particles.velocity[idx] = {0.0, 0.0};  // 静止流体
      fluid_particles.pressure[idx] = 0.0;  // 稍后通过PPE求解
      fluid_particles.density[idx] = 1000.0;  // 水的密度
      fluid_particles.surface_type[idx] = SurfaceType::INNER;
      ++idx;
    }
  }
  
  // 生成固体粒子（容器壁面）- 两层壁面粒子
  // 底部壁面
  int nx_bottom = static_cast<int>(container_width / particle_spacing) + 1;
  // 左侧壁面
  int ny_left = static_cast<int>(container_height / particle_spacing) + 1;
  // 右侧壁面
  int ny_right = ny_left;
  
  // 两层壁面：每层都有底部、左侧、右侧
  int num_solid = 2 * (nx_bottom + ny_left + ny_right);
  solid_particles.particle_num = num_solid;
  solid_particles.position.resize(num_solid);
  solid_particles.velocity.resize(num_solid);
  solid_particles.normal_vector.resize(num_solid);
  
  idx = 0;
  
  // 底部壁面 - 第一层（法向量向上）
  for (int i = 0; i < nx_bottom; ++i) {
    double x = i * particle_spacing;
    double y = -particle_spacing;  // 距离底部流体粒子一个粒子间距
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {0.0, 1.0};  // 向上
    ++idx;
  }
  
  // 底部壁面 - 第二层（法向量向上）
  for (int i = 0; i < nx_bottom; ++i) {
    double x = i * particle_spacing;
    double y = -2.0 * particle_spacing;  // 第二层，再向外一个粒子间距
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {0.0, 1.0};  // 向上
    ++idx;
  }
  
  // 左侧壁面 - 第一层（法向量向右）
  for (int j = 0; j < ny_left; ++j) {
    double x = -particle_spacing;
    double y = j * particle_spacing;
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {1.0, 0.0};  // 向右
    ++idx;
  }
  
  // 左侧壁面 - 第二层（法向量向右）
  for (int j = 0; j < ny_left; ++j) {
    double x = -2.0 * particle_spacing;  // 第二层，再向外一个粒子间距
    double y = j * particle_spacing;
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {1.0, 0.0};  // 向右
    ++idx;
  }
  
  // 右侧壁面 - 第一层（法向量向左）
  for (int j = 0; j < ny_right; ++j) {
    double x = container_width;
    double y = j * particle_spacing;
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {-1.0, 0.0};  // 向左
    ++idx;
  }
  
  // 右侧壁面 - 第二层（法向量向左）
  for (int j = 0; j < ny_right; ++j) {
    double x = container_width + particle_spacing;  // 第二层，再向外一个粒子间距
    double y = j * particle_spacing;
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {-1.0, 0.0};  // 向左
    ++idx;
  }
}

int main() {
  std::cout << "=== PPE求解模块静水压力测试 ===" << std::endl;
  
  // 物理参数
  const double rho = 1000.0;  // 水的密度 (kg/m³)
  const double g = 9.8;       // 重力加速度 (m/s²)
  const double gravity_x = 0.0;
  const double gravity_y = -g;
  
  // 几何参数
  const double container_width = 6;   // 容器宽度 (m)
  const double container_height = 6;  // 容器高度 (m)
  const double water_height = 3;      // 水位高度 (m)
  
  // 粒子参数（平滑半径保持为粒子间距的2.1倍）
  const double particle_spacing = M_PI / 64;  // 粒子间距 (m)
  const double particle_radius = particle_spacing / 2.0;
  const double smoothing_radius = 3.1 * particle_spacing;  // 平滑半径是粒子间距的2.1倍
  const double cell_size = 2.0 * smoothing_radius;
  
  // 时间步长
  const double time_step = 0.001;  // 时间步长 (s)
  
  std::cout << "\n物理参数:" << std::endl;
  std::cout << "  流体密度: " << rho << " kg/m³" << std::endl;
  std::cout << "  重力加速度: " << g << " m/s²" << std::endl;
  
  std::cout << "\n几何参数:" << std::endl;
  std::cout << "  容器宽度: " << container_width << " m" << std::endl;
  std::cout << "  容器高度: " << container_height << " m" << std::endl;
  std::cout << "  水位高度: " << water_height << " m" << std::endl;
  
  std::cout << "\n粒子参数:" << std::endl;
  std::cout << "  粒子间距: " << particle_spacing << " m" << std::endl;
  std::cout << "  平滑半径: " << smoothing_radius << " m" << std::endl;
  std::cout << "  网格尺寸: " << cell_size << " m" << std::endl;
    
  // 生成测试场景
  FluidParticle fluid_particles("fluid");
  SolidParticle solid_particles("solid");
  
  std::cout << "\n生成测试场景..." << std::endl;
  GenerateHydrostaticTest(
      container_width, container_height, water_height,
      particle_spacing, fluid_particles, solid_particles);
  
  std::cout << "  流体粒子数: " << fluid_particles.particle_num << std::endl;
  std::cout << "  固体粒子数: " << solid_particles.particle_num << std::endl;
  std::cout << "  总粒子数: " << fluid_particles.particle_num + solid_particles.particle_num << std::endl;
  
  // 构建邻域列表
  std::cout << "\n构建邻域列表..." << std::endl;
  NeighborListSearcher neighbor_searcher;
  neighbor_searcher.BuildNeighborList(
      fluid_particles, solid_particles,
      particle_radius, smoothing_radius, cell_size);
  
  // 统计邻域信息
  int total_fluid_neighbors = 0;
  int total_solid_neighbors = 0;
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    total_fluid_neighbors += fluid_particles.fluid_neighbour_list[i].size();
    total_solid_neighbors += fluid_particles.solid_neighbour_list[i].size();
  }
  std::cout << "  平均流体邻域粒子数: " 
            << static_cast<double>(total_fluid_neighbors) / fluid_particles.particle_num << std::endl;
  std::cout << "  平均固体邻域粒子数: " 
            << static_cast<double>(total_solid_neighbors) / fluid_particles.particle_num << std::endl;
  
  // 检测表面粒子
  std::cout << "\n检测表面粒子..." << std::endl;
  SurfaceDetector surface_detector;
  surface_detector.DetectSurfaceParticles(
      fluid_particles, solid_particles, smoothing_radius, particle_spacing);
  
  // 统计表面粒子数量
  int num_inner = 0, num_near_surface = 0, num_surface = 0, num_splash = 0;
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    switch (fluid_particles.surface_type[i]) {
      case SurfaceType::INNER:
        ++num_inner;
        break;
      case SurfaceType::NEAR_SURFACE:
        ++num_near_surface;
        break;
      case SurfaceType::SURFACE:
        ++num_surface;
        break;
      case SurfaceType::SPLASH:
        ++num_splash;
        break;
    }
  }
  std::cout << "  内部粒子数: " << num_inner << std::endl;
  std::cout << "  近表面粒子数: " << num_near_surface << std::endl;
  std::cout << "  表面粒子数: " << num_surface << std::endl;
  std::cout << "  飞溅粒子数: " << num_splash << std::endl;
  
  // 计算corrective matrix
  std::cout << "\n计算corrective matrix..." << std::endl;
  CorrectiveMatrix corrective_matrix_calc;
  std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>> 
      corrective_matrices(fluid_particles.particle_num);
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    corrective_matrices[i] = corrective_matrix_calc.ComputeCorrectiveMatrix(
        i, fluid_particles, solid_particles, smoothing_radius, false);
  }
  std::cout << "  完成" << std::endl;
  
  // 直接构建PETSc格式的PPE系数矩阵和右边项
  std::cout << "\n构建PPE系数矩阵和右边项（PETSc格式）..." << std::endl;
  PPEMatrixBuilder matrix_builder;
  Mat A_petsc = NULL;
  Vec b_petsc = NULL;
  
  bool success = matrix_builder.BuildPPEMatrixPetsc(
      fluid_particles, solid_particles, corrective_matrices,
      smoothing_radius, rho, time_step, particle_spacing, gravity_x, gravity_y, A_petsc, b_petsc);
  
  if (!success) {
    std::cerr << "错误：构建PPE矩阵失败" << std::endl;
    return 1;
  }
    
  // 获取矩阵信息
  PetscInt m, n;
  MatGetSize(A_petsc, &m, &n);
  MatInfo info;
  MatGetInfo(A_petsc, MAT_GLOBAL_SUM, &info);
  PetscInt nnz = static_cast<PetscInt>(info.nz_used);
  
  std::cout << "  矩阵大小: " << m << " x " << n << std::endl;
  std::cout << "  非零元素数: " << nnz << std::endl;
  std::cout << "  稀疏度: " << (1.0 - static_cast<double>(nnz) / (m * n)) * 100.0 
            << "%" << std::endl;
  
  // 分析矩阵和向量
  std::cout << "\n分析矩阵和向量..." << std::endl;
  
  // 提取对角线元素
  int num_fluid_particles = fluid_particles.particle_num;
  std::vector<double> diagonal_elements(num_fluid_particles);
  std::vector<double> surface_diagonal;
  std::vector<double> inner_diagonal;
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    PetscInt row = static_cast<PetscInt>(i);
    PetscInt col = static_cast<PetscInt>(i);
    PetscScalar value;
    MatGetValues(A_petsc, 1, &row, 1, &col, &value);
    diagonal_elements[i] = value;
    
    if (fluid_particles.surface_type[i] == SurfaceType::SURFACE) {
      surface_diagonal.push_back(value);
    } else {
      inner_diagonal.push_back(value);
    }
  }
  
  // 统计对角线元素
  if (!diagonal_elements.empty()) {
    auto minmax = std::minmax_element(diagonal_elements.begin(), diagonal_elements.end());
    double diag_sum = 0.0;
    for (double val : diagonal_elements) {
      diag_sum += std::abs(val);
    }
    std::cout << "  对角线元素统计:" << std::endl;
    std::cout << "    最小值: " << *minmax.first << std::endl;
    std::cout << "    最大值: " << *minmax.second << std::endl;
    std::cout << "    平均值: " << diag_sum / diagonal_elements.size() << std::endl;
  }
  
  if (!surface_diagonal.empty()) {
    auto minmax = std::minmax_element(surface_diagonal.begin(), surface_diagonal.end());
    std::cout << "  表面粒子对角线元素:" << std::endl;
    std::cout << "    最小值: " << *minmax.first << std::endl;
    std::cout << "    最大值: " << *minmax.second << std::endl;
    std::cout << "    数量: " << surface_diagonal.size() << std::endl;
  }
  
  if (!inner_diagonal.empty()) {
    auto minmax = std::minmax_element(inner_diagonal.begin(), inner_diagonal.end());
    std::cout << "  内部粒子对角线元素:" << std::endl;
    std::cout << "    最小值: " << *minmax.first << std::endl;
    std::cout << "    最大值: " << *minmax.second << std::endl;
    std::cout << "    数量: " << inner_diagonal.size() << std::endl;
  }
  
  // 提取右边项
  std::vector<PetscInt> indices(num_fluid_particles);
  std::vector<PetscScalar> values(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
    indices[i] = static_cast<PetscInt>(i);
  }
  VecGetValues(b_petsc, num_fluid_particles, indices.data(), values.data());
  
  std::vector<double> rhs_values(num_fluid_particles);
  std::vector<double> surface_rhs;
  std::vector<double> inner_rhs;
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    rhs_values[i] = values[i];
    if (fluid_particles.surface_type[i] == SurfaceType::SURFACE) {
      surface_rhs.push_back(values[i]);
    } else {
      inner_rhs.push_back(values[i]);
    }
  }
  
  // 统计右边项
  if (!rhs_values.empty()) {
    auto minmax = std::minmax_element(rhs_values.begin(), rhs_values.end());
    double rhs_sum = 0.0;
    for (double val : rhs_values) {
      rhs_sum += std::abs(val);
    }
    std::cout << "  右边项统计:" << std::endl;
    std::cout << "    最小值: " << *minmax.first << std::endl;
    std::cout << "    最大值: " << *minmax.second << std::endl;
    std::cout << "    平均值: " << rhs_sum / rhs_values.size() << std::endl;
  }
  
  if (!surface_rhs.empty()) {
    auto minmax = std::minmax_element(surface_rhs.begin(), surface_rhs.end());
    std::cout << "  表面粒子右边项:" << std::endl;
    std::cout << "    最小值: " << *minmax.first << std::endl;
    std::cout << "    最大值: " << *minmax.second << std::endl;
  }
  
  if (!inner_rhs.empty()) {
    auto minmax = std::minmax_element(inner_rhs.begin(), inner_rhs.end());
    std::cout << "  内部粒子右边项:" << std::endl;
    std::cout << "    最小值: " << *minmax.first << std::endl;
    std::cout << "    最大值: " << *minmax.second << std::endl;
  }
  
  // 输出调试信息到VTK文件
  std::cout << "\n输出调试信息到VTK文件..." << std::endl;
  std::string debug_filename = "output/ppe_debug_hydrostatic.vtk";
  if (matrix_builder.WriteDebugInfoToVTK(fluid_particles, A_petsc, b_petsc, debug_filename)) {
    std::cout << "  调试信息已保存到: " << debug_filename << std::endl;
  } else {
    std::cerr << "  警告：无法保存调试信息到VTK文件" << std::endl;
  }
  
  std::cout << "\nPPE矩阵构建完成！" << std::endl;
  
  // 求解PPE方程
  std::cout << "\n求解PPE方程..." << std::endl;
  PPESolver::SolverConfig solver_config;
  solver_config.solver_type = PPESolver::SolverType::BICGSTAB;
  solver_config.max_iterations = 1000;
  solver_config.tolerance = 1e-6;
  solver_config.force_iterative = true;  // 强制使用迭代方法
  
  PPESolver ppe_solver(solver_config);
  Vec p_petsc = NULL;  // 压力解向量
  
  if (!ppe_solver.Solve(A_petsc, b_petsc, p_petsc)) {
    std::cerr << "错误：PPE求解失败" << std::endl;
    MatDestroy(&A_petsc);
    VecDestroy(&b_petsc);
    if (p_petsc != NULL) {
      VecDestroy(&p_petsc);
    }
    return 1;
  }
  
  // 输出求解信息
  std::cout << "  求解完成" << std::endl;
  std::cout << "  迭代次数: " << ppe_solver.GetLastIterations() << std::endl;
  std::cout << "  残差: " << ppe_solver.GetLastResidual() << std::endl;
  std::cout << "  是否收敛: " << (ppe_solver.GetLastConverged() ? "是" : "否") << std::endl;
  
  // 提取压力值并更新到粒子
  std::cout << "\n提取压力解..." << std::endl;
  // 重用之前的indices向量
  for (int i = 0; i < num_fluid_particles; ++i) {
    indices[i] = static_cast<PetscInt>(i);
  }
  std::vector<PetscScalar> pressure_values(num_fluid_particles);
  VecGetValues(p_petsc, num_fluid_particles, indices.data(), pressure_values.data());
  
  // 更新粒子压力值
  for (int i = 0; i < num_fluid_particles; ++i) {
    fluid_particles.pressure[i] = pressure_values[i];
  }
  
  // 统计压力值
  std::vector<double> surface_pressure;
  std::vector<double> inner_pressure;
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    if (fluid_particles.surface_type[i] == SurfaceType::SURFACE) {
      surface_pressure.push_back(pressure_values[i]);
    } else {
      inner_pressure.push_back(pressure_values[i]);
    }
  }
  
  if (!pressure_values.empty()) {
    auto minmax = std::minmax_element(pressure_values.begin(), pressure_values.end());
    double pressure_sum = 0.0;
    for (double val : pressure_values) {
      pressure_sum += std::abs(val);
    }
    std::cout << "  压力统计:" << std::endl;
    std::cout << "    最小值: " << *minmax.first << " Pa" << std::endl;
    std::cout << "    最大值: " << *minmax.second << " Pa" << std::endl;
    std::cout << "    平均值: " << pressure_sum / pressure_values.size() << " Pa" << std::endl;
  }
  
  if (!surface_pressure.empty()) {
    auto minmax = std::minmax_element(surface_pressure.begin(), surface_pressure.end());
    std::cout << "  表面粒子压力:" << std::endl;
    std::cout << "    最小值: " << *minmax.first << " Pa" << std::endl;
    std::cout << "    最大值: " << *minmax.second << " Pa" << std::endl;
    std::cout << "    数量: " << surface_pressure.size() << std::endl;
  }
  
  if (!inner_pressure.empty()) {
    auto minmax = std::minmax_element(inner_pressure.begin(), inner_pressure.end());
    std::cout << "  内部粒子压力:" << std::endl;
    std::cout << "    最小值: " << *minmax.first << " Pa" << std::endl;
    std::cout << "    最大值: " << *minmax.second << " Pa" << std::endl;
    std::cout << "    数量: " << inner_pressure.size() << std::endl;
  }
  
  // 计算理论静水压力（p = rho * g * h，其中h是从自由表面的深度）
  // 自由表面在y = water_height，所以深度 h = water_height - y
  std::cout << "\n验证静水压力分布..." << std::endl;
  double max_pressure_error = 0.0;
  double avg_pressure_error = 0.0;
  int error_count = 0;
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    double y = fluid_particles.position[i].y;
    double depth = water_height - y;  // 从自由表面的深度
    double theoretical_pressure = rho * g * depth;  // 理论静水压力
    double computed_pressure = pressure_values[i];
    double error = std::abs(computed_pressure - theoretical_pressure);
    
    max_pressure_error = std::max(max_pressure_error, error);
    avg_pressure_error += error;
    ++error_count;
  }
  
  if (error_count > 0) {
    avg_pressure_error /= error_count;
    std::cout << "  最大压力误差: " << max_pressure_error << " Pa" << std::endl;
    std::cout << "  平均压力误差: " << avg_pressure_error << " Pa" << std::endl;
    std::cout << "  相对误差: " << (avg_pressure_error / (rho * g * water_height)) * 100.0 << "%" << std::endl;
  }
  
  // 使用LSMPS计算压力梯度
  std::cout << "\n使用LSMPS计算压力梯度..." << std::endl;
  std::vector<double2> pressure_gradient(num_fluid_particles);
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    const double2& pos_i = fluid_particles.position[i];
    double p_i = pressure_values[i];
    
    // 提取corrective matrix的前两行（用于梯度计算）
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C1 = corrective_matrices[i].row(0);  // x方向
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C2 = corrective_matrices[i].row(1);  // y方向
    
    double grad_x = 0.0;
    double grad_y = 0.0;
    
    // 处理流体邻域粒子
    for (int j : fluid_particles.fluid_neighbour_list[i]) {
      const double2& pos_j = fluid_particles.position[j];
      double p_j = pressure_values[j];
      
      double dx = pos_j.x - pos_i.x;
      double dy = pos_j.y - pos_i.y;
      double dist = ComputeDistance(pos_i, pos_j);
      
      if (dist < 1e-10 || dist > smoothing_radius) {
        continue;
      }
      
      double d_ij = (p_j - p_i) / dist;
      double weight = WeightFunction(dist, smoothing_radius);
      
      // 计算基函数
      Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis = 
          corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
      
      // 计算梯度贡献
      grad_x += weight * d_ij * (C1 * basis)(0, 0);
      grad_y += weight * d_ij * (C2 * basis)(0, 0);
    }
    
    // 对于压力梯度计算，通常只需要考虑流体邻域粒子
    // 壁面边界条件已经通过corrective matrix的构建得到了考虑
    
    pressure_gradient[i] = {grad_x, grad_y};
  }
  
  // 统计压力梯度
  std::vector<double> grad_magnitude(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
    grad_magnitude[i] = std::sqrt(pressure_gradient[i].x * pressure_gradient[i].x + 
                                   pressure_gradient[i].y * pressure_gradient[i].y);
  }
  
  if (!grad_magnitude.empty()) {
    auto minmax = std::minmax_element(grad_magnitude.begin(), grad_magnitude.end());
    double grad_sum = 0.0;
    for (double val : grad_magnitude) {
      grad_sum += val;
    }
    std::cout << "  压力梯度统计:" << std::endl;
    std::cout << "    最小梯度: " << *minmax.first << " Pa/m" << std::endl;
    std::cout << "    最大梯度: " << *minmax.second << " Pa/m" << std::endl;
    std::cout << "    平均梯度: " << grad_sum / grad_magnitude.size() << " Pa/m" << std::endl;
  }
  
  // 验证理论压力梯度（对于静水压力，dp/dy = -rho*g，dp/dx = 0）
  std::cout << "\n验证压力梯度分布..." << std::endl;
  double max_grad_error = 0.0;
  double avg_grad_error = 0.0;
  int grad_error_count = 0;
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    // 理论压力梯度：dp/dx = 0, dp/dy = -rho*g（负号因为y向下为正）
    double theoretical_grad_x = 0.0;
    double theoretical_grad_y = -rho * g;
    
    double error_x = std::abs(pressure_gradient[i].x - theoretical_grad_x);
    double error_y = std::abs(pressure_gradient[i].y - theoretical_grad_y);
    double error_magnitude = std::sqrt(error_x * error_x + error_y * error_y);
    
    max_grad_error = std::max(max_grad_error, error_magnitude);
    avg_grad_error += error_magnitude;
    ++grad_error_count;
  }
  
  if (grad_error_count > 0) {
    avg_grad_error /= grad_error_count;
    std::cout << "  最大梯度误差: " << max_grad_error << " Pa/m" << std::endl;
    std::cout << "  平均梯度误差: " << avg_grad_error << " Pa/m" << std::endl;
    std::cout << "  理论梯度大小: " << (rho * g) << " Pa/m" << std::endl;
    std::cout << "  相对误差: " << (avg_grad_error / (rho * g)) * 100.0 << "%" << std::endl;
  }
  
  // 输出结果到VTK文件
  std::cout << "\n输出压力结果到VTK文件..." << std::endl;
  std::string result_filename = "output/ppe_result_hydrostatic.vtk";
  FileOperator file_op;
  if (file_op.writeVTKBase(result_filename, fluid_particles)) {
    // 追加压力标量
    std::vector<double> pressure_vec(num_fluid_particles);
    for (int i = 0; i < num_fluid_particles; ++i) {
      pressure_vec[i] = pressure_values[i];
    }
    if (file_op.appendVTKScalar(result_filename, "pressure", pressure_vec)) {
      std::cout << "  压力标量已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加压力标量到VTK文件" << std::endl;
    }
    
    // 追加压力梯度向量
    if (file_op.appendVTKVector(result_filename, "pressure_gradient", pressure_gradient)) {
      std::cout << "  压力梯度向量已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加压力梯度向量到VTK文件" << std::endl;
    }
    
    // 追加梯度大小标量
    if (file_op.appendVTKScalar(result_filename, "gradient_magnitude", grad_magnitude)) {
      std::cout << "  梯度大小标量已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加梯度大小标量到VTK文件" << std::endl;
    }
    
    std::cout << "  结果已保存到: " << result_filename << std::endl;
  } else {
    std::cerr << "  警告：无法创建结果VTK文件" << std::endl;
  }
  
  // 清理PETSc对象
  MatDestroy(&A_petsc);
  VecDestroy(&b_petsc);
  if (p_petsc != NULL) {
    VecDestroy(&p_petsc);
  }
  
  return 0;
}
