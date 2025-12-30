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
#include <iomanip>
#include <map>
#include <chrono>
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
  
  // 计算粒子数量，确保关于中心对称
  // 从中心开始，向两边均匀分布
  double center_x = container_width / 2.0;
  double center_y = 0.0;  // 从底部开始
  
  // 计算x方向的粒子数（确保对称）
  int nx_half = static_cast<int>(std::round(center_x / particle_spacing));
  int nx_fluid = 2 * nx_half + 1;  // 中心粒子 + 左右各nx_half个
  
  // 计算y方向的粒子数
  int ny_fluid = static_cast<int>(std::round(water_height / particle_spacing)) + 1;
  
  int num_fluid = nx_fluid * ny_fluid;
  
  // 生成流体粒子（水）- 关于中心对称
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
      double y = j * particle_spacing;
    // 从中心开始，向两边对称分布
    for (int i = -nx_half; i <= nx_half; ++i) {
      double x = center_x + i * particle_spacing;
      fluid_particles.position[idx] = {x, y};
      fluid_particles.velocity[idx] = {0.0, 0.0};  // 静止流体
      fluid_particles.pressure[idx] = 0.0;  // 稍后通过PPE求解
      fluid_particles.density[idx] = 1000.0;  // 水的密度
      fluid_particles.surface_type[idx] = SurfaceType::INNER;
      ++idx;
    }
  }
  
  // 生成固体粒子（容器壁面）- 单层壁面粒子，完全对称
  // 底部壁面：与流体粒子底部保持一个粒子间距
  int nx_bottom = nx_fluid;  // 与流体粒子x方向数量相同，确保对称
  
  // 左侧和右侧壁面：与流体粒子高度相同
  int ny_wall = ny_fluid;  // 与流体粒子y方向数量相同
  
  // 单层壁面：底部、左侧、右侧
  int num_solid = nx_bottom + 2 * ny_wall;  // 底部 + 左右各一层
  solid_particles.particle_num = num_solid;
  solid_particles.position.resize(num_solid);
  solid_particles.velocity.resize(num_solid);
  solid_particles.normal_vector.resize(num_solid);
  
  idx = 0;
  
  // 底部壁面（法向量向上）
  // 与流体粒子底部保持一个粒子间距，x坐标与流体粒子对齐
  for (int i = -nx_half; i <= nx_half; ++i) {
    double x = center_x + i * particle_spacing;
    double y = -particle_spacing;  // 距离底部流体粒子一个粒子间距
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {0.0, 1.0};  // 向上
    ++idx;
  }
  
  // 左侧壁面（法向量向右）
  // 与流体粒子左侧保持一个粒子间距，y坐标与流体粒子对齐
  double left_x = center_x - (nx_half + 1) * particle_spacing;  // 左侧壁面x坐标
  for (int j = 0; j < ny_wall; ++j) {
    double x = left_x;
    double y = j * particle_spacing;
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {1.0, 0.0};  // 向右
    ++idx;
  }
  
  // 右侧壁面（法向量向左）
  // 与流体粒子右侧保持一个粒子间距，y坐标与流体粒子对齐
  double right_x = center_x + (nx_half + 1) * particle_spacing;  // 右侧壁面x坐标
  for (int j = 0; j < ny_wall; ++j) {
    double x = right_x;
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
  const double particle_spacing = 0.03;  // 粒子间距 (m)
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
  
  // 检测表面粒子
  std::cout << "\n检测表面粒子..." << std::endl;
  SurfaceDetector surface_detector;
  surface_detector.DetectSurfaceParticles(
      fluid_particles, solid_particles, smoothing_radius, particle_spacing);
  
  // 计算corrective matrix
  std::cout << "\n计算corrective matrix..." << std::endl;
  CorrectiveMatrix corrective_matrix_calc;
  // 速度散度使用第一类边界条件（border_condition = false）
  std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>> 
      corrective_matrices_velocity(fluid_particles.particle_num);
  // 压力拉普拉斯算子使用第二类边界条件（border_condition = true）
  std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>> 
      corrective_matrices_pressure(fluid_particles.particle_num);
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    // 速度散度：第一类边界条件
    corrective_matrices_velocity[i] = corrective_matrix_calc.ComputeCorrectiveMatrix(
        i, fluid_particles, solid_particles, smoothing_radius, false);
    // 压力拉普拉斯算子：第二类边界条件
    corrective_matrices_pressure[i] = corrective_matrix_calc.ComputeCorrectiveMatrix(
        i, fluid_particles, solid_particles, smoothing_radius, true);
  }
  std::cout << "  完成（速度：第一类边界条件，压力：第二类边界条件）" << std::endl;
  
  // 直接构建PETSc格式的PPE系数矩阵和右边项
  std::cout << "\n构建PPE系数矩阵和右边项（PETSc格式）..." << std::endl;
  PPEMatrixBuilder matrix_builder;
  Mat A_petsc = NULL;
  Vec b_petsc = NULL;
  
  bool success = matrix_builder.BuildPPEMatrixPetsc(
      fluid_particles, solid_particles, corrective_matrices_velocity, corrective_matrices_pressure,
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
  
  int num_fluid_particles = fluid_particles.particle_num;
  
  // 提取右边项
  std::vector<PetscInt> indices(num_fluid_particles);
  std::vector<PetscScalar> values(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
    indices[i] = static_cast<PetscInt>(i);
  }
  VecGetValues(b_petsc, num_fluid_particles, indices.data(), values.data());
  
  
  std::cout << "\nPPE矩阵构建完成！" << std::endl;
  
  // 使用罚函数方法构建KKT系统
  std::cout << "\n使用罚函数方法构建KKT系统..." << std::endl;
  Mat K_petsc = NULL;
  Vec f_petsc = NULL;
  double penalty_parameter = 1e3;  // 罚函数参数μ
  
  bool penalty_success = matrix_builder.BuildPenaltySystem(
      A_petsc, b_petsc, fluid_particles, penalty_parameter, K_petsc, f_petsc);
  
  if (!penalty_success) {
    std::cerr << "错误：构建罚函数系统失败" << std::endl;
    MatDestroy(&A_petsc);
    VecDestroy(&b_petsc);
    return 1;
  }
  
  // 统计自由面粒子数量
  int num_surface_particles = 0;
  for (int i = 0; i < num_fluid_particles; ++i) {
    if (fluid_particles.surface_type[i] == SurfaceType::SURFACE) {
      num_surface_particles++;
    }
  }
  
  std::cout << "  罚函数参数 μ: " << penalty_parameter << std::endl;
  std::cout << "  自由面粒子数: " << num_surface_particles << std::endl;
  
  // 获取K矩阵信息
  PetscInt K_m, K_n;
  MatGetSize(K_petsc, &K_m, &K_n);
  MatInfo K_info;
  MatGetInfo(K_petsc, MAT_GLOBAL_SUM, &K_info);
  PetscInt K_nnz = static_cast<PetscInt>(K_info.nz_used);
  
  std::cout << "  K矩阵大小: " << K_m << " x " << K_n << std::endl;
  std::cout << "  K矩阵非零元素数: " << K_nnz << std::endl;
  std::cout << "  K矩阵稀疏度: " << (1.0 - static_cast<double>(K_nnz) / (K_m * K_n)) * 100.0 
            << "%" << std::endl;
  
  // 求解PPE方程（使用罚函数系统 K·p = f）
  // 注意：K = A^T A + D 是对称正定矩阵，应使用专门的求解器
  std::cout << "\n========================================" << std::endl;
  std::cout << "对比测试：CG vs GMRES 求解器性能" << std::endl;
  std::cout << "========================================" << std::endl;
  
  // 测试结果结构
  struct SolverResult {
    std::string name;
    int iterations;
    double residual;
    bool converged;
    double time_seconds;
  };
  
  std::vector<SolverResult> results;
  
  // 测试1：CG方法（对称正定矩阵专用）
  {
    std::cout << "\n【测试1】CG方法（适用于对称正定矩阵）..." << std::endl;
    PPESolver::SolverConfig solver_config;
    solver_config.solver_type = PPESolver::SolverType::CG;
    solver_config.max_iterations = 5000;
    solver_config.tolerance = 1e-6;
    solver_config.force_iterative = true;
    solver_config.is_symmetric_positive_definite = true;
    
    PPESolver ppe_solver(solver_config);
    Vec p_petsc_cg = NULL;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    bool success = ppe_solver.Solve(K_petsc, f_petsc, p_petsc_cg);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    double time_seconds = duration.count() / 1e6;
    
    SolverResult result;
    result.name = "CG";
    result.iterations = ppe_solver.GetLastIterations();
    result.residual = ppe_solver.GetLastResidual();
    result.converged = ppe_solver.GetLastConverged();
    result.time_seconds = time_seconds;
    results.push_back(result);
    
    std::cout << "  求解完成" << std::endl;
    std::cout << "  迭代次数: " << result.iterations << std::endl;
    std::cout << "  残差: " << result.residual << std::endl;
    std::cout << "  是否收敛: " << (result.converged ? "是" : "否") << std::endl;
    std::cout << "  求解时间: " << result.time_seconds << " 秒" << std::endl;
    
    if (!success) {
      std::cerr << "错误：CG求解失败" << std::endl;
      if (p_petsc_cg != NULL) {
        VecDestroy(&p_petsc_cg);
      }
    }
  }
  
  // 测试2：GMRES方法（通用方法，也可用于对称正定矩阵）
  {
    std::cout << "\n【测试2】GMRES方法（通用方法）..." << std::endl;
    PPESolver::SolverConfig solver_config;
    solver_config.solver_type = PPESolver::SolverType::GMRES;
    solver_config.max_iterations = 5000;
    solver_config.tolerance = 1e-6;
    solver_config.restart = 30;
    solver_config.force_iterative = true;
    solver_config.is_symmetric_positive_definite = false;  // GMRES不要求对称正定
    
    PPESolver ppe_solver(solver_config);
    Vec p_petsc_gmres = NULL;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    bool success = ppe_solver.Solve(K_petsc, f_petsc, p_petsc_gmres);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    double time_seconds = duration.count() / 1e6;
    
    SolverResult result;
    result.name = "GMRES";
    result.iterations = ppe_solver.GetLastIterations();
    result.residual = ppe_solver.GetLastResidual();
    result.converged = ppe_solver.GetLastConverged();
    result.time_seconds = time_seconds;
    results.push_back(result);
    
    std::cout << "  求解完成" << std::endl;
    std::cout << "  迭代次数: " << result.iterations << std::endl;
    std::cout << "  残差: " << result.residual << std::endl;
    std::cout << "  是否收敛: " << (result.converged ? "是" : "否") << std::endl;
    std::cout << "  求解时间: " << result.time_seconds << " 秒" << std::endl;
    
    if (!success) {
      std::cerr << "错误：GMRES求解失败" << std::endl;
      if (p_petsc_gmres != NULL) {
        VecDestroy(&p_petsc_gmres);
      }
    }
  }
  
  // 输出对比结果
  std::cout << "\n========================================" << std::endl;
  std::cout << "性能对比总结" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << std::left << std::setw(12) << "方法" 
            << std::setw(15) << "迭代次数" 
            << std::setw(15) << "残差" 
            << std::setw(12) << "收敛" 
            << std::setw(15) << "求解时间(秒)" << std::endl;
  std::cout << std::string(70, '-') << std::endl;
  
  for (const auto& result : results) {
    std::cout << std::left << std::setw(12) << result.name
              << std::setw(15) << result.iterations
              << std::setw(15) << std::scientific << std::setprecision(6) << result.residual
              << std::setw(12) << (result.converged ? "是" : "否")
              << std::setw(15) << std::fixed << std::setprecision(4) << result.time_seconds << std::endl;
  }
  
  // 选择最佳方法
  if (results.size() == 2) {
    const auto& cg_result = results[0];
    const auto& gmres_result = results[1];
    
    std::cout << "\n结论：" << std::endl;
    if (cg_result.converged && gmres_result.converged) {
      if (cg_result.time_seconds < gmres_result.time_seconds) {
        std::cout << "  CG方法更快，快 " << (gmres_result.time_seconds / cg_result.time_seconds) 
                  << " 倍" << std::endl;
      } else {
        std::cout << "  GMRES方法更快，快 " << (cg_result.time_seconds / gmres_result.time_seconds) 
                  << " 倍" << std::endl;
      }
      
      if (cg_result.iterations < gmres_result.iterations) {
        std::cout << "  CG方法迭代次数更少（" << cg_result.iterations << " vs " 
                  << gmres_result.iterations << "）" << std::endl;
      } else {
        std::cout << "  GMRES方法迭代次数更少（" << gmres_result.iterations << " vs " 
                  << cg_result.iterations << "）" << std::endl;
      }
    } else if (cg_result.converged) {
      std::cout << "  CG方法收敛，GMRES方法未收敛" << std::endl;
    } else if (gmres_result.converged) {
      std::cout << "  GMRES方法收敛，CG方法未收敛" << std::endl;
  } else {
      std::cout << "  两种方法都未收敛" << std::endl;
    }
  }
  
  // 使用CG的结果作为最终解（因为CG是专门为对称正定矩阵设计的）
  Vec p_petsc = NULL;
  if (results[0].converged) {
    // 使用CG的结果
    std::cout << "\n使用CG方法的解作为最终结果" << std::endl;
    // 注意：这里需要重新求解或复制结果，为了简化，我们重新求解一次
  PPESolver::SolverConfig solver_config;
    solver_config.solver_type = PPESolver::SolverType::CG;
    solver_config.max_iterations = 5000;
  solver_config.tolerance = 1e-6;
    solver_config.force_iterative = true;
    solver_config.is_symmetric_positive_definite = true;
  
  PPESolver ppe_solver(solver_config);
    if (!ppe_solver.Solve(K_petsc, f_petsc, p_petsc)) {
      std::cerr << "错误：最终PPE求解失败" << std::endl;
      MatDestroy(&A_petsc);
      VecDestroy(&b_petsc);
      MatDestroy(&K_petsc);
      VecDestroy(&f_petsc);
      return 1;
    }
  } else {
    std::cerr << "错误：CG求解未收敛，无法继续" << std::endl;
    MatDestroy(&A_petsc);
    VecDestroy(&b_petsc);
    MatDestroy(&K_petsc);
    VecDestroy(&f_petsc);
    return 1;
  }
  
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
  
  // 计算理论静水压力（p = rho * g * h，其中h是从自由表面的深度）
  // 自由表面在y = water_height，所以深度 h = water_height - y
  
  // 创建理论压力向量 X'
  Vec x_theoretical = NULL;
  VecCreate(PETSC_COMM_WORLD, &x_theoretical);
  VecSetSizes(x_theoretical, PETSC_DECIDE, num_fluid_particles);
  VecSetType(x_theoretical, VECSEQ);
  VecSetFromOptions(x_theoretical);
  
  // 设置理论压力值
  std::vector<PetscScalar> theoretical_pressure_vec(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
    double y = fluid_particles.position[i].y;
    double depth = water_height - y;  // 从自由表面的深度
    double theoretical_pressure = rho * g * depth;  // 理论静水压力
    theoretical_pressure_vec[i] = theoretical_pressure;
  }
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    VecSetValue(x_theoretical, i, theoretical_pressure_vec[i], INSERT_VALUES);
  }
  VecAssemblyBegin(x_theoretical);
  VecAssemblyEnd(x_theoretical);
  
  // 计算 b' = A * X'
  Vec b_computed = NULL;
  VecCreate(PETSC_COMM_WORLD, &b_computed);
  VecSetSizes(b_computed, PETSC_DECIDE, num_fluid_particles);
  VecSetType(b_computed, VECSEQ);
  VecSetFromOptions(b_computed);
  
  MatMult(A_petsc, x_theoretical, b_computed);
  
  // 提取 b' 的值
  std::vector<PetscScalar> b_computed_values(num_fluid_particles);
  VecGetValues(b_computed, num_fluid_particles, indices.data(), b_computed_values.data());
  
  // 提取实际右边项 b 的值（之前已经提取过，但为了清晰重新提取）
  std::vector<PetscScalar> b_actual_values(num_fluid_particles);
  VecGetValues(b_petsc, num_fluid_particles, indices.data(), b_actual_values.data());
  
  // 计算差异
  std::vector<double> diff_values(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
    diff_values[i] = std::abs(b_computed_values[i] - b_actual_values[i]);
  }
  
  // 清理临时向量
  VecDestroy(&x_theoretical);
  VecDestroy(&b_computed);
  
  // 使用LSMPS计算压力梯度（第二类边界条件）
  std::cout << "\n使用LSMPS计算压力梯度..." << std::endl;
  std::vector<double2> pressure_gradient(num_fluid_particles);
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    const double2& pos_i = fluid_particles.position[i];
    double p_i = pressure_values[i];
    
    // 提取压力corrective matrix的前两行（用于梯度计算，第二类边界条件）
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C1 = corrective_matrices_pressure[i].row(0);  // x方向
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C2 = corrective_matrices_pressure[i].row(1);  // y方向
    
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
      
      // 计算基函数（标准基函数）
      Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis = 
          corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
      
      // 计算梯度贡献
      grad_x += weight * d_ij * (C1 * basis)(0, 0);
      grad_y += weight * d_ij * (C2 * basis)(0, 0);
    }
    
    // 处理壁面邻域粒子（第二类边界条件）
    // 壁面处压力梯度的法向分量：dp/dn = -rho * g * n_y（对于静水压力）
    for (int j : fluid_particles.solid_neighbour_list[i]) {
      const double2& pos_j = solid_particles.position[j];
      const double2& normal = solid_particles.normal_vector[j];
      
      double dx = pos_j.x - pos_i.x;
      double dy = pos_j.y - pos_i.y;
      double dist = ComputeDistance(pos_i, pos_j);
      
      if (dist < 1e-10 || dist > smoothing_radius) {
        continue;
      }
      
      double weight = WeightFunction(dist, smoothing_radius);
      
      // 使用壁面基函数（第二类边界条件）
      Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis_wall = 
          corrective_matrix_calc.ComputeBasisFunctionsForWall(
              dx, dy, normal.x, normal.y, smoothing_radius);
      
      // 壁面处压力梯度的法向分量：dp/dn = -rho * g * n_y（对于静水压力）
      // 归一化法向量
      double nn = std::sqrt(normal.x * normal.x + normal.y * normal.y);
      double n_y = (nn > 1e-10) ? normal.y / nn : 0.0;
      double d_ij = -rho * g * n_y;
      
      // 计算梯度贡献
      grad_x += weight * d_ij * (C1 * basis_wall)(0, 0);
      grad_y += weight * d_ij * (C2 * basis_wall)(0, 0);
    }
    
    pressure_gradient[i] = {grad_x, grad_y};
  }
  
  // 计算梯度大小
  std::vector<double> grad_magnitude(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
    grad_magnitude[i] = std::sqrt(pressure_gradient[i].x * pressure_gradient[i].x + 
                                   pressure_gradient[i].y * pressure_gradient[i].y);
  }
  
  // 输出所有流体粒子信息到同一个VTK文件
  std::cout << "\n输出所有流体粒子信息到VTK文件..." << std::endl;
  std::string result_filename = "output/ppe_result_hydrostatic.vtk";
  FileOperator file_op;
  
  // 输出流体粒子基础信息
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
    
    // 追加理论压力
    std::vector<double> theoretical_pressure_double(num_fluid_particles);
    for (int i = 0; i < num_fluid_particles; ++i) {
      theoretical_pressure_double[i] = static_cast<double>(theoretical_pressure_vec[i]);
    }
    if (file_op.appendVTKScalar(result_filename, "theoretical_pressure", theoretical_pressure_double)) {
      std::cout << "  理论压力已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加理论压力到VTK文件" << std::endl;
    }
    
    // 追加 b' = A * X' 的计算结果
    std::vector<double> b_computed_double(num_fluid_particles);
    for (int i = 0; i < num_fluid_particles; ++i) {
      b_computed_double[i] = static_cast<double>(b_computed_values[i]);
    }
    if (file_op.appendVTKScalar(result_filename, "b_computed_AX", b_computed_double)) {
      std::cout << "  b' = A * X' 已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加 b' 到 VTK 文件" << std::endl;
    }
    
    // 追加差异 |b' - b|
    if (file_op.appendVTKScalar(result_filename, "diff_b_computed_actual", diff_values)) {
      std::cout << "  差异 |b' - b| 已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加差异到 VTK 文件" << std::endl;
    }
    
    // 追加自由面类型信息
    std::vector<int> surface_type_values(num_fluid_particles);
    for (int i = 0; i < num_fluid_particles; ++i) {
      // 将SurfaceType枚举转换为整数
      // INNER = 0, NEAR_SURFACE = 1, SURFACE = 2, SPLASH = 3
      surface_type_values[i] = static_cast<int>(fluid_particles.surface_type[i]);
    }
    if (file_op.appendVTKScalar(result_filename, "surface_type", surface_type_values)) {
      std::cout << "  自由面类型已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加自由面类型到VTK文件" << std::endl;
    }
    
    // 追加调试信息（对角线元素和右边项）
    if (matrix_builder.WriteDebugInfoToVTK(fluid_particles, A_petsc, b_petsc, result_filename)) {
      std::cout << "  调试信息（对角线元素、右边项等）已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加调试信息到VTK文件" << std::endl;
    }
    
    std::cout << "  所有流体粒子信息已保存到: " << result_filename << std::endl;
  } else {
    std::cerr << "  警告：无法创建流体粒子VTK文件" << std::endl;
  }
  
  // 单独输出壁面粒子到新的VTK文件
  std::string wall_filename = "output/ppe_result_hydrostatic_wall.vtk";
  std::cout << "\n输出壁面粒子到VTK文件..." << std::endl;
  std::ofstream wall_file(wall_filename);
  if (!wall_file.is_open()) {
    std::cerr << "  警告：无法创建壁面粒子VTK文件" << std::endl;
  } else {
    int num_solid = solid_particles.particle_num;
    
    wall_file << std::fixed << std::setprecision(15);
    
    // VTK文件头
    wall_file << "# vtk DataFile Version 3.0\n";
    wall_file << "MPS Particle Data 2D - Wall Particles\n";
    wall_file << "ASCII\n";
    wall_file << "DATASET POLYDATA\n";
    
    // 写入点坐标
    wall_file << "POINTS " << num_solid << " float\n";
    for (int i = 0; i < num_solid; ++i) {
      const auto& pos = solid_particles.position[i];
      wall_file << pos.x << " " << pos.y << " 0.0\n";
    }
    
    // 写入顶点
    wall_file << "VERTICES " << num_solid << " " << (num_solid * 2) << "\n";
    for (int i = 0; i < num_solid; ++i) {
      wall_file << "1 " << i << "\n";
    }
    
    // 写入点数据
    wall_file << "POINT_DATA " << num_solid << "\n";
    
    // 写入速度向量
    wall_file << "VECTORS velocity float\n";
    for (int i = 0; i < num_solid; ++i) {
      const auto& vel = solid_particles.velocity[i];
      wall_file << vel.x << " " << vel.y << " 0.0\n";
    }
    
    // 写入法向量
    wall_file << "VECTORS normal_vector float\n";
    for (int i = 0; i < num_solid; ++i) {
      const auto& normal = solid_particles.normal_vector[i];
      wall_file << normal.x << " " << normal.y << " 0.0\n";
    }
    
    // 写入法向量X分量
    wall_file << "SCALARS normal_vector_x float\n";
    wall_file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < num_solid; ++i) {
      wall_file << solid_particles.normal_vector[i].x << "\n";
    }
    
    // 写入法向量Y分量
    wall_file << "SCALARS normal_vector_y float\n";
    wall_file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < num_solid; ++i) {
      wall_file << solid_particles.normal_vector[i].y << "\n";
    }
    
    wall_file.close();
    std::cout << "  壁面粒子结果已保存到: " << wall_filename << " (包含 " 
              << num_solid << " 个壁面粒子)" << std::endl;
  }
  
  // 验证自由面粒子的压力（应该接近0）
  if (num_surface_particles > 0) {
    std::cout << "\n验证自由面粒子压力约束..." << std::endl;
    double max_surface_pressure = 0.0;
    double sum_surface_pressure = 0.0;
    int surface_count = 0;
    
    for (int i = 0; i < num_fluid_particles; ++i) {
      if (fluid_particles.surface_type[i] == SurfaceType::SURFACE) {
        double p = pressure_values[i];
        double abs_p = std::abs(p);
        if (abs_p > max_surface_pressure) {
          max_surface_pressure = abs_p;
        }
        sum_surface_pressure += abs_p;
        surface_count++;
      }
    }
    
    double avg_surface_pressure = (surface_count > 0) ? sum_surface_pressure / surface_count : 0.0;
    std::cout << "  自由面粒子数: " << num_surface_particles << std::endl;
    std::cout << "  自由面粒子最大压力绝对值: " << max_surface_pressure << std::endl;
    std::cout << "  自由面粒子平均压力绝对值: " << avg_surface_pressure << std::endl;
    std::cout << "  罚函数参数 μ: " << penalty_parameter << std::endl;
    
    if (max_surface_pressure < 1.0 / penalty_parameter * 10.0) {
      std::cout << "  ✓ 自由面压力约束满足（压力接近0）" << std::endl;
    } else {
      std::cout << "  ⚠ 警告：自由面压力约束可能未完全满足" << std::endl;
    }
  } else {
    std::cout << "\n注意：当前测试场景中没有自由面粒子（完全封闭容器）" << std::endl;
  }
  
  // 清理PETSc对象
  MatDestroy(&A_petsc);
  VecDestroy(&b_petsc);
  MatDestroy(&K_petsc);
  VecDestroy(&f_petsc);
  if (p_petsc != NULL) {
    VecDestroy(&p_petsc);
  }
  
  return 0;
}
