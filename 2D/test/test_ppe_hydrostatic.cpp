#include "../src/PPE/PPEMatrixBuilder.hpp"
#include "../src/PPE/PPESolver.hpp"
#include "../src/lsmps/CorrectiveMatrix.hpp"
#include "../src/neighbour_list/NeighborListSearcher.hpp"
#include "../src/surface_detection/SurfaceDetector.hpp"
#include "../src/core/Particle.hpp"
#include "../src/core/FileOperator.hpp"
#include "../src/explicit_force/ExplicitForce.hpp"
#include "core/Types.h"
#include "core/MPSUtils.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <map>
#include <limits>
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
    double /* container_height */,  // 未使用，保留用于接口兼容性
    double water_height,
    double particle_spacing,
    FluidParticle& fluid_particles,
    SolidParticle& solid_particles) {
  
  // 计算粒子数量，确保关于中心对称
  // 从中心开始，向两边均匀分布
  double center_x = container_width / 2.0;
  // center_y = 0.0 从底部开始（未使用，保留注释）
  
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

int main(int argc, char** argv) {
  // 在程序入口统一初始化PETSc（包括MPI）
  PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr);
  if (ierr != 0) {
    std::cerr << "错误：PETSc初始化失败" << std::endl;
    return 1;
  }
  std::cout << "=== PPE求解模块静水压力测试 ===" << std::endl;
  
  // 物理参数
  const double rho = 1000.0;  // 水的密度 (kg/m³)
  const double g = 9.8;       // 重力加速度 (m/s²)
  const double gravity_x = 0.0;
  const double gravity_y = -g;
  
  // 几何参数
  const double container_width = 1.0;   // 容器宽度 (m)
  const double container_height = 2.0;   // 容器高度 (m)
  const double water_height = 1.0;      // 水位高度 (m)
  
  // 粒子参数（平滑半径保持为粒子间距的3.1倍）
  const double particle_spacing = 0.02;  // 粒子间距 (m)
  const double particle_radius = particle_spacing / 2.0;
  const double smoothing_radius = 3.1 * particle_spacing;  // 平滑半径是粒子间距的3.1倍
  const double cell_size = 2.0 * smoothing_radius;
  
  // 时间步长
  const double time_step = 0.01;  // 时间步长 (s)
  
  // 粘性参数（用于显式力计算）
  const double kinematic_viscosity = 1.0e-6;  // 运动粘性系数 (m²/s)，水的典型值
  
  std::cout << "\n物理参数:" << std::endl;
  std::cout << "  流体密度: " << rho << " kg/m³" << std::endl;
  std::cout << "  重力加速度: " << g << " m/s²" << std::endl;
  std::cout << "  运动粘性系数: " << kinematic_viscosity << " m²/s" << std::endl;
  
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
  
  // 计算corrective matrix（用于显式力计算，第一类边界条件）
  std::cout << "\n计算corrective matrix（显式力用，第一类边界条件）..." << std::endl;
  CorrectiveMatrix corrective_matrix_calc;
  std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>> 
      corrective_matrices_explicit(fluid_particles.particle_num);
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    // 显式力计算使用第一类边界条件
    corrective_matrices_explicit[i] = corrective_matrix_calc.ComputeCorrectiveMatrix(
        i, fluid_particles, solid_particles, smoothing_radius, false);
  }
  std::cout << "  完成（第一类边界条件）" << std::endl;
  
  // ========== 显式更新模块：计算粘性力和重力，更新临时速度 ==========
  std::cout << "\n========== 显式更新模块：计算粘性力和重力，更新临时速度 ==========" << std::endl;
  ExplicitForce explicit_force;
  std::vector<double2> viscous_acceleration(fluid_particles.particle_num);

  // 保存初始速度（用于对比）
  std::vector<double2> initial_velocity(fluid_particles.particle_num);
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    initial_velocity[i] = fluid_particles.velocity[i];
  }

  // 计算粘性力和重力，更新流体粒子的临时速度，并返回粘性加速度
  explicit_force.ComputeAndUpdateVelocity(
      fluid_particles, solid_particles,
      corrective_matrices_explicit,
      smoothing_radius,
      kinematic_viscosity,
      gravity_x, gravity_y,
      time_step,
      viscous_acceleration);

  // 保存显式更新后的速度（临时速度）
  std::vector<double2> velocity_explicit(fluid_particles.particle_num);
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    velocity_explicit[i] = fluid_particles.velocity[i];
  }

  // 壁面速度在 PPEMatrixBuilder::BuildPPEMatrixPetsc 内会临时加上 Δt*g 用于速度散度计算，计算完后恢复，此处无需预设

  std::cout << "  显式更新完成（速度已更新为临时速度，位置未改变）" << std::endl;
  
  // ========== 重新计算corrective matrix（用于PPE，使用更新后的速度）==========
  // 注意：显式更新只更新速度，不更新位置，因此用于速度散度的corrective matrix
  // 与显式力计算时使用的相同（考虑壁面粒子，第一类边界条件），可以直接复用corrective_matrices_explicit
  // 只需要计算用于压力梯度的corrective matrix（第二类边界条件）
  std::cout << "\n计算corrective matrix（PPE用，压力梯度，第二类边界条件）..." << std::endl;
  std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>> 
      corrective_matrices_pressure(fluid_particles.particle_num);
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    // 压力拉普拉斯算子：第二类边界条件
    corrective_matrices_pressure[i] = corrective_matrix_calc.ComputeCorrectiveMatrix(
        i, fluid_particles, solid_particles, smoothing_radius, true);
  }
  std::cout << "  完成（第二类边界条件）" << std::endl;
  
  // 复用显式力计算时的corrective matrix作为速度散度用（第一类边界条件，考虑壁面粒子）
  // 因为显式更新不改变粒子位置，所以可以直接复用
  const auto& corrective_matrices_velocity = corrective_matrices_explicit;
  
  // 直接构建PETSc格式的PPE系数矩阵和右边项（使用显式更新后的速度）
  std::cout << "\n构建PPE系数矩阵和右边项（PETSc格式，使用显式更新后的速度）..." << std::endl;
  PPEMatrixBuilder matrix_builder;
  Mat A_petsc = NULL;
  Vec b_petsc = NULL;
  const double penalty_mu = 0.0;  // 使用默认罚参数尺度
  
  bool success = matrix_builder.BuildPPEMatrixPetsc(
      fluid_particles, solid_particles, corrective_matrices_velocity, corrective_matrices_pressure,
      smoothing_radius, rho, time_step, particle_spacing, penalty_mu,
      gravity_x, gravity_y, A_petsc, b_petsc);
  
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
  
#if 0  // 诊断统计代码（已禁用）
  // ========== 诊断：检查矩阵和右边项的状态 ==========
  std::cout << "\n[诊断] 检查矩阵和右边项的状态..." << std::endl;
  
  // 检查矩阵的对角线元素
  double min_diag = std::numeric_limits<double>::max();
  double max_diag = std::numeric_limits<double>::lowest();
  int zero_diag_count = 0;
  
  for (PetscInt i = 0; i < m && i < 1000; ++i) {  // 只检查前1000行以节省时间
    PetscScalar diag_val;
    MatGetValue(A_petsc, i, i, &diag_val);
    double diag = std::abs(static_cast<double>(diag_val));
    min_diag = std::min(min_diag, diag);
    max_diag = std::max(max_diag, diag);
    if (diag < 1e-15) {
      zero_diag_count++;
    }
  }
  
  // 检查右边项的范数
  PetscReal b_norm;
  VecNorm(b_petsc, NORM_2, &b_norm);
  
  std::cout << "  矩阵A规模: " << m << " x " << n << std::endl;
  std::cout << "  矩阵A对角线元素范围（前1000行）: [" << std::scientific << std::setprecision(3) 
            << min_diag << ", " << max_diag << "]" << std::fixed << std::endl;
  std::cout << "  零对角线元素数量（前1000行）: " << zero_diag_count << std::endl;
  std::cout << "  右边项b的范数: " << std::scientific << std::setprecision(3) 
            << b_norm << std::fixed << std::endl;
  
  if (zero_diag_count > 0) {
    std::cerr << "  ⚠️  警告：矩阵A有零对角线元素，可能导致求解困难！" << std::endl;
  }
  if (b_norm < 1e-15) {
    std::cerr << "  ⚠️  警告：右边项b的范数极小（接近零）！" << std::endl;
  }
  // ========== 诊断代码结束 ==========
#endif
  
  // ========== 提取对角线元素（用于检查压力拉普拉斯算子的离散）==========
  std::cout << "\n[分析] 提取对角线元素（检查压力拉普拉斯算子离散）..." << std::endl;
  std::vector<double> diagonal_elements(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
    PetscInt row = static_cast<PetscInt>(i);
    PetscInt col = static_cast<PetscInt>(i);
    PetscScalar diag_val;
    PetscErrorCode ierr = MatGetValue(A_petsc, row, col, &diag_val);
    if (ierr != 0) {
      std::cerr << "警告：无法获取矩阵对角线元素 " << i << std::endl;
      diagonal_elements[i] = 0.0;
    } else {
      diagonal_elements[i] = static_cast<double>(diag_val);
    }
  }
  
  // 对角线元素统计代码已移除，仅保留对角线本身用于后处理
  
  // ========== 计算右边项中速度散度部分（用于后处理）==========
  std::vector<double> divergence_contribution(num_fluid_particles, 0.0);
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    const double2& pos_i = fluid_particles.position[i];
    const double2& vel_i = fluid_particles.velocity[i];
    
    // 提取速度corrective matrix的行向量（用于速度散度计算，第一类边界条件）
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C1_velocity = corrective_matrices_velocity[i].row(0);
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C2_velocity = corrective_matrices_velocity[i].row(1);
    
    // 提取压力corrective matrix的行向量（用于压力拉普拉斯算子计算，第二类边界条件）
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> M2_pressure = corrective_matrices_pressure[i].row(2);
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> M3_pressure = corrective_matrices_pressure[i].row(3);
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> M2_plus_M3_pressure = M2_pressure + M3_pressure;
    
    double divergence_sum = 0.0;
    
    // 遍历流体邻域粒子
    for (int j : fluid_particles.fluid_neighbour_list[i]) {
      const double2& pos_j = fluid_particles.position[j];
      const double2& vel_j = fluid_particles.velocity[j];
      
      double dx = pos_j.x - pos_i.x;
      double dy = pos_j.y - pos_i.y;
      double dist = ComputeDistance(pos_i, pos_j);
      
      if (dist < 1e-10 || dist > smoothing_radius) {
        continue;
      }
      
      double weight = WeightFunction(dist, smoothing_radius);
      
      // 计算基函数
      Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis = 
          corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
      
      // 计算速度散度项（使用速度corrective matrix，第一类边界条件）
      double dvx = vel_j.x - vel_i.x;
      double dvy = vel_j.y - vel_i.y;
      double C1P = (C1_velocity * basis)(0, 0);
      double C2P = (C2_velocity * basis)(0, 0);
      divergence_sum += weight * (C1P * dvx + C2P * dvy);
    }
    
    // 处理壁面邻域粒子
    for (int j : fluid_particles.solid_neighbour_list[i]) {
      const double2& pos_j = solid_particles.position[j];
      const double2& normal = solid_particles.normal_vector[j];
      // 在速度散度计算中，壁面粒子的速度临时设为 delta_t * g
      double2 vel_wall_effective = {gravity_x * time_step, gravity_y * time_step};
      
      double dx = pos_j.x - pos_i.x;
      double dy = pos_j.y - pos_i.y;
      double dist = ComputeDistance(pos_i, pos_j);
      
      if (dist < 1e-10 || dist > smoothing_radius) {
        continue;
      }
      
      double weight = WeightFunction(dist, smoothing_radius);
      
      // 对于速度散度：使用第一类边界条件的基函数（标准基函数）
      Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis_velocity = 
          corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
      
      // 计算速度散度项（使用速度corrective matrix，第一类边界条件）
      double dvx = vel_wall_effective.x - vel_i.x;
      double dvy = vel_wall_effective.y - vel_i.y;
      double C1P = (C1_velocity * basis_velocity)(0, 0);
      double C2P = (C2_velocity * basis_velocity)(0, 0);
      divergence_sum += weight * (C1P * dvx + C2P * dvy);
      
    }
    
    // 计算贡献
    divergence_contribution[i] = (1.0 / (smoothing_radius * time_step)) * divergence_sum;
  }
  
  // ========== 求解PPE方程（直接求解 A·p = b）==========
  std::cout << "\n求解PPE方程（直接求解 A·p = b）..." << std::endl;
  PPESolver::SolverConfig solver_config;
  solver_config.solver_type = PPESolver::SolverType::CG;
  solver_config.max_iterations = 10000;  // 最大迭代次数
  solver_config.tolerance = 1e-6;        // 容差
  solver_config.force_iterative = true;  // 强制使用迭代求解器
  solver_config.is_symmetric_positive_definite = true;
  solver_config.restart = 30;
  
  PPESolver ppe_solver(solver_config);
  Vec p_petsc = NULL;  // 压力解向量
  
  if (!ppe_solver.Solve(A_petsc, b_petsc, p_petsc)) {
    std::cerr << "错误：PPE求解失败" << std::endl;
    if (p_petsc != NULL) {
      VecDestroy(&p_petsc);
    }
    if (b_petsc != NULL) {
      VecDestroy(&b_petsc);
    }
    if (A_petsc != NULL) {
      MatDestroy(&A_petsc);
    }
    PetscFinalize();
    return 1;
  }
  
  std::cout << "  求解完成" << std::endl;
  std::cout << "  迭代次数: " << ppe_solver.GetLastIterations() << std::endl;
  std::cout << "  残差: " << ppe_solver.GetLastResidual() << std::endl;
  std::cout << "  是否收敛: " << (ppe_solver.GetLastConverged() ? "是" : "否") << std::endl;
  
  // 提取压力解并写回粒子
  std::vector<PetscScalar> pressure_values(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
    indices[i] = static_cast<PetscInt>(i);
  }
  VecGetValues(p_petsc, num_fluid_particles, indices.data(), pressure_values.data());
  for (int i = 0; i < num_fluid_particles; ++i) {
    fluid_particles.pressure[i] = static_cast<double>(pressure_values[i]);
  }
  
  // ========== 使用LSMPS计算压力梯度 ==========
  std::cout << "\n使用LSMPS计算压力梯度..." << std::endl;
  std::vector<double2> pressure_gradient(num_fluid_particles);
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    const double2& pos_i = fluid_particles.position[i];
    double p_i = static_cast<double>(pressure_values[i]);
    
    // 提取 [M_{i,0}, M_{i,1}]
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> M0 = corrective_matrices_pressure[i].row(0);  // x方向
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> M1 = corrective_matrices_pressure[i].row(1);  // y方向
    
    double sum_x = 0.0;
    double sum_y = 0.0;
    
    // 流体邻域
    for (int j : fluid_particles.fluid_neighbour_list[i]) {
      const double2& pos_j = fluid_particles.position[j];
      double p_j = static_cast<double>(pressure_values[j]);
      
      double dx = pos_j.x - pos_i.x;
      double dy = pos_j.y - pos_i.y;
      double dist = ComputeDistance(pos_i, pos_j);
      if (dist < 1e-10 || dist > smoothing_radius) {
        continue;
      }
      
      double dp = p_j - p_i;
      double weight = WeightFunction(dist, smoothing_radius);
      
      Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis =
          corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
      
      double M0P = (M0 * basis)(0, 0);
      double M1P = (M1 * basis)(0, 0);
      sum_x += weight * dp * M0P;
      sum_y += weight * dp * M1P;
    }
    
    // 壁面邻域（第二类边界条件）
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
      Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis_wall =
          corrective_matrix_calc.ComputeBasisFunctionsForWall(
              dx, dy, normal.x, normal.y, smoothing_radius);
      
      double nn = std::sqrt(normal.x * normal.x + normal.y * normal.y);
      if (nn < 1e-10) {
        continue;
      }
      double n_x = normal.x / nn;
      double n_y = normal.y / nn;
      double dp_dn = rho * (n_x * gravity_x + n_y * gravity_y);
      
      double M0Q = (M0 * basis_wall)(0, 0);
      double M1Q = (M1 * basis_wall)(0, 0);
      sum_x += weight * smoothing_radius * dp_dn * M0Q;
      sum_y += weight * smoothing_radius * dp_dn * M1Q;
    }
    
    double inv_re = 1.0 / smoothing_radius;
    pressure_gradient[i] = {inv_re * sum_x, inv_re * sum_y};
  }
  
  // 梯度模长
  std::vector<double> grad_magnitude(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
    grad_magnitude[i] = std::sqrt(pressure_gradient[i].x * pressure_gradient[i].x +
                                   pressure_gradient[i].y * pressure_gradient[i].y);
  }
  
  // ========== 输出结果到VTK文件 ==========
  std::cout << "\n输出所有流体粒子信息到VTK文件..." << std::endl;
  std::string result_filename = "output/ppe_result_hydrostatic.vtk";
  FileOperator file_op;
  
  if (file_op.writeVTKBase(result_filename, fluid_particles)) {
    // 位置和速度由 writeVTKBase 写入
    
    // 压力
    std::vector<double> pressure_out(num_fluid_particles);
    for (int i = 0; i < num_fluid_particles; ++i) {
      pressure_out[i] = static_cast<double>(pressure_values[i]);
    }
    file_op.appendVTKScalar(result_filename, "pressure", pressure_out);
    
    // 压力梯度
    file_op.appendVTKVector(result_filename, "pressure_gradient", pressure_gradient);
    file_op.appendVTKScalar(result_filename, "gradient_magnitude", grad_magnitude);
    
    // 速度散度
    file_op.appendVTKScalar(result_filename, "velocity_divergence", divergence_contribution);
    
    // 对角线元素大小 |A_ii|
    std::vector<double> diag_abs(num_fluid_particles);
    for (int i = 0; i < num_fluid_particles; ++i) {
      diag_abs[i] = std::abs(diagonal_elements[i]);
    }
    file_op.appendVTKScalar(result_filename, "diag_abs", diag_abs);
    
    // 右边项大小 |b_i|
    std::vector<double> rhs_abs(num_fluid_particles);
    for (int i = 0; i < num_fluid_particles; ++i) {
      rhs_abs[i] = std::abs(static_cast<double>(values[i]));
    }
    file_op.appendVTKScalar(result_filename, "rhs_abs", rhs_abs);
    
    // 自由面类型
    std::vector<int> surface_type_values(num_fluid_particles);
    for (int i = 0; i < num_fluid_particles; ++i) {
      surface_type_values[i] = static_cast<int>(fluid_particles.surface_type[i]);
    }
    file_op.appendVTKScalar(result_filename, "surface_type", surface_type_values);
    
    std::cout << "  结果文件已保存到: " << result_filename << std::endl;
  } else {
    std::cerr << "  错误：无法创建结果VTK文件" << std::endl;
  }
  
  // 清理PETSc对象
  if (p_petsc != NULL) {
    VecDestroy(&p_petsc);
  }
  if (b_petsc != NULL) {
    VecDestroy(&b_petsc);
  }
  if (A_petsc != NULL) {
    MatDestroy(&A_petsc);
  }
  
  // 在程序结束前统一结束PETSc（包括MPI）
  PetscFinalize();
  
  return 0;
}
