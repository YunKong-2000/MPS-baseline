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

int main() {
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
  
  // 使用重载方法，直接获取explicit_force模块计算出的粘性力加速度
  explicit_force.ComputeAndUpdateVelocity(
      fluid_particles, solid_particles,
      corrective_matrices_explicit,
      smoothing_radius,
      kinematic_viscosity,
      gravity_x, gravity_y,
      time_step,
      viscous_acceleration);
  
  // 保存显式更新后的速度
  std::vector<double2> velocity_explicit(fluid_particles.particle_num);
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    velocity_explicit[i] = fluid_particles.velocity[i];
  }
  
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
  
  // 分析对角线元素的统计信息
  double min_diag_all = std::numeric_limits<double>::max();
  double max_diag_all = std::numeric_limits<double>::lowest();
  double sum_abs_diag = 0.0;
  int negative_diag_count = 0;
  int zero_diag_count_all = 0;
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    double diag = diagonal_elements[i];
    double abs_diag = std::abs(diag);
    if (diag < min_diag_all) min_diag_all = diag;
    if (diag > max_diag_all) max_diag_all = diag;
    sum_abs_diag += abs_diag;
    if (diag < 0.0) negative_diag_count++;
    if (abs_diag < 1e-15) zero_diag_count_all++;
  }
  
  double avg_abs_diag = sum_abs_diag / num_fluid_particles;
  
  std::cout << "  对角线元素统计（所有粒子）:" << std::endl;
  std::cout << "    最小值: " << std::scientific << std::setprecision(3) << min_diag_all << std::fixed << std::endl;
  std::cout << "    最大值: " << std::scientific << std::setprecision(3) << max_diag_all << std::fixed << std::endl;
  std::cout << "    平均值（绝对值）: " << std::scientific << std::setprecision(3) << avg_abs_diag << std::fixed << std::endl;
  std::cout << "    负值数量: " << negative_diag_count << " / " << num_fluid_particles << std::endl;
  std::cout << "    零值数量: " << zero_diag_count_all << " / " << num_fluid_particles << std::endl;
  
  // 理论值：对于压力拉普拉斯算子，对角线元素应该是负值
  // 当前实现与 PPEMatrixBuilder 一致：
  // 对角线系数 = - (2/(r_e^2 * ρ)) * Σ_{j∈fluid} w_ij * ([M_{i,2}+M_{i,3}] P_ij)
  if (negative_diag_count == 0) {
    std::cerr << "  ⚠️  警告：所有对角线元素都是非负值，可能不符合压力拉普拉斯算子的离散形式！" << std::endl;
  }
  if (zero_diag_count_all > 0) {
    std::cerr << "  ⚠️  警告：有 " << zero_diag_count_all << " 个零对角线元素，可能导致矩阵奇异！" << std::endl;
  }
  // ========== 对角线元素分析结束 ==========
  
  // ========== 计算右边项的详细贡献 ==========
  std::cout << "\n[分析] 计算右边项的详细贡献..." << std::endl;
  
  // 系数因子（与BuildPPEMatrixPetsc中使用的相同）
  // 2 / (r_e^2 * ρ)
  double coeff_factor = 2.0 / (smoothing_radius * smoothing_radius * rho);
  
  // 计算每个粒子的速度散度贡献和压力边界条件贡献
  std::vector<double> divergence_contribution(num_fluid_particles, 0.0);
  std::vector<double> wall_pressure_contribution(num_fluid_particles, 0.0);
  std::vector<double> b_total_recomputed(num_fluid_particles, 0.0);
  
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
    double wall_pressure_term = 0.0;
    
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
      const double2& vel_wall = solid_particles.velocity[j];
      const double2& normal = solid_particles.normal_vector[j];
      
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
      
      // 对于壁面压力边界条件项：使用第二类边界条件的基函数（壁面基函数）
      Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis_pressure = 
          corrective_matrix_calc.ComputeBasisFunctionsForWall(
              dx, dy, normal.x, normal.y, smoothing_radius);
      
      // 计算速度散度项（使用速度corrective matrix，第一类边界条件）
      double dvx = vel_wall.x - vel_i.x;
      double dvy = vel_wall.y - vel_i.y;
      double C1P = (C1_velocity * basis_velocity)(0, 0);
      double C2P = (C2_velocity * basis_velocity)(0, 0);
      divergence_sum += weight * (C1P * dvx + C2P * dvy);
      
      // 计算壁面压力边界条件项（使用压力corrective matrix，第二类边界条件）
      double n_dot_g = normal.x * gravity_x + normal.y * gravity_y;
      // 文档壁面项：w_ij * r_e * ρ (n·g) * ([M2+M3] Q_ij)
      double wall_pressure_coeff = weight * smoothing_radius * rho * n_dot_g *
          (M2_plus_M3_pressure * basis_pressure)(0, 0);
      wall_pressure_term += wall_pressure_coeff;
    }
    
    // 计算贡献
    divergence_contribution[i] = (1.0 / (smoothing_radius * time_step)) * divergence_sum;
    wall_pressure_contribution[i] = coeff_factor * wall_pressure_term;
    b_total_recomputed[i] = divergence_contribution[i] + wall_pressure_contribution[i];
  }
  
  // 统计信息
  double max_div_contrib = 0.0;
  double max_wall_contrib = 0.0;
  double sum_abs_div = 0.0;
  double sum_abs_wall = 0.0;
  for (int i = 0; i < num_fluid_particles; ++i) {
    double abs_div = std::abs(divergence_contribution[i]);
    double abs_wall = std::abs(wall_pressure_contribution[i]);
    if (abs_div > max_div_contrib) max_div_contrib = abs_div;
    if (abs_wall > max_wall_contrib) max_wall_contrib = abs_wall;
    sum_abs_div += abs_div;
    sum_abs_wall += abs_wall;
  }
  double avg_abs_div = sum_abs_div / num_fluid_particles;
  double avg_abs_wall = sum_abs_wall / num_fluid_particles;
  
  std::cout << "  速度散度贡献统计:" << std::endl;
  std::cout << "    最大值: " << std::scientific << std::setprecision(3) << max_div_contrib << std::fixed << std::endl;
  std::cout << "    平均值: " << std::scientific << std::setprecision(3) << avg_abs_div << std::fixed << std::endl;
  std::cout << "  压力边界条件贡献统计:" << std::endl;
  std::cout << "    最大值: " << std::scientific << std::setprecision(3) << max_wall_contrib << std::fixed << std::endl;
  std::cout << "    平均值: " << std::scientific << std::setprecision(3) << avg_abs_wall << std::fixed << std::endl;
  
  // 验证重新计算的b与实际b的差异
  double max_b_diff = 0.0;
  double sum_b_diff = 0.0;
  for (int i = 0; i < num_fluid_particles; ++i) {
    double b_actual = static_cast<double>(values[i]);
    double b_recomputed = b_total_recomputed[i];
    double diff = std::abs(b_actual - b_recomputed);
    if (diff > max_b_diff) max_b_diff = diff;
    sum_b_diff += diff;
  }
  double avg_b_diff = sum_b_diff / num_fluid_particles;
  std::cout << "  重新计算的b与实际b的差异:" << std::endl;
  std::cout << "    最大差异: " << std::scientific << std::setprecision(3) << max_b_diff << std::fixed << std::endl;
  std::cout << "    平均差异: " << std::scientific << std::setprecision(3) << avg_b_diff << std::fixed << std::endl;
  
  // ========== 分析：为什么右边项b的改变会影响收敛性 ==========
  std::cout << "\n[分析] 为什么右边项b的改变会影响迭代求解器的收敛性？" << std::endl;
  std::cout << "  虽然系数矩阵A不变，但右边项b的改变会影响收敛性的原因：" << std::endl;
  
  // 1. 右边项范数的影响
  PetscReal b_norm_actual;
  VecNorm(b_petsc, NORM_2, &b_norm_actual);
  double b_norm_recomputed = 0.0;
  for (int i = 0; i < num_fluid_particles; ++i) {
    b_norm_recomputed += b_total_recomputed[i] * b_total_recomputed[i];
  }
  b_norm_recomputed = std::sqrt(b_norm_recomputed);
  
  std::cout << "\n  1. 相对残差的影响：" << std::endl;
  std::cout << "    右边项b的范数: " << std::scientific << std::setprecision(3) 
            << b_norm_actual << std::fixed << std::endl;
  std::cout << "    迭代求解器使用相对残差: ||r|| / ||b||" << std::endl;
  std::cout << "    如果b的范数变化，即使绝对残差相同，相对残差也会不同" << std::endl;
  const double target_tolerance = 1e-6;  // 与后面配置的容差一致
  std::cout << "    目标相对残差: " << target_tolerance << std::endl;
  std::cout << "    对应的绝对残差目标: " << std::scientific << std::setprecision(3) 
            << (target_tolerance * b_norm_actual) << std::fixed << std::endl;
  
  // 2. 初始残差的影响
  // 初始残差 r0 = b - A*p0，其中p0 = 0（零初始猜测）
  // 所以 r0 = b，初始残差范数 = ||b||
  std::cout << "\n  2. 初始残差的影响：" << std::endl;
  std::cout << "    初始猜测: p0 = 0（零向量）" << std::endl;
  std::cout << "    初始残差: r0 = b - A*p0 = b" << std::endl;
  std::cout << "    初始残差范数: ||r0|| = ||b|| = " << std::scientific << std::setprecision(3) 
            << b_norm_actual << std::fixed << std::endl;
  std::cout << "    如果b的范数很大，初始残差也很大，需要更多迭代才能收敛" << std::endl;
  std::cout << "    如果b的范数很小，初始残差也很小，但可能接近机器精度，导致数值问题" << std::endl;
  
  // 3. 数值精度的影响
  double max_abs_b = 0.0;
  double min_abs_b = std::numeric_limits<double>::max();
  for (int i = 0; i < num_fluid_particles; ++i) {
    double abs_b = std::abs(static_cast<double>(values[i]));
    if (abs_b > max_abs_b) max_abs_b = abs_b;
    if (abs_b < min_abs_b && abs_b > 1e-20) min_abs_b = abs_b;  // 忽略真正的零
  }
  double b_range = (min_abs_b < std::numeric_limits<double>::max()) ? (max_abs_b / min_abs_b) : 0.0;
  
  std::cout << "\n  3. 数值精度的影响：" << std::endl;
  std::cout << "    右边项b的元素范围: [" << std::scientific << std::setprecision(3) 
            << min_abs_b << ", " << max_abs_b << "]" << std::fixed << std::endl;
  std::cout << "    范围比值: " << std::scientific << std::setprecision(3) 
            << b_range << std::fixed << std::endl;
  if (b_range > 1e10) {
    std::cerr << "    ⚠️  警告：右边项b的元素范围很大，可能导致数值精度问题！" << std::endl;
    std::cerr << "    在矩阵-向量乘法中，大值和小值的混合可能导致舍入误差累积" << std::endl;
  }
  
  // 4. 速度散度贡献的分布
  double max_abs_div = 0.0;
  double max_abs_wall = 0.0;
  for (int i = 0; i < num_fluid_particles; ++i) {
    double abs_div = std::abs(divergence_contribution[i]);
    double abs_wall = std::abs(wall_pressure_contribution[i]);
    if (abs_div > max_abs_div) max_abs_div = abs_div;
    if (abs_wall > max_abs_wall) max_abs_wall = abs_wall;
  }
  
  std::cout << "\n  4. 右边项组成的影响：" << std::endl;
  std::cout << "    速度散度贡献最大值: " << std::scientific << std::setprecision(3) 
            << max_abs_div << std::fixed << std::endl;
  std::cout << "    压力边界条件贡献最大值: " << std::scientific << std::setprecision(3) 
            << max_abs_wall << std::fixed << std::endl;
  if (max_abs_div > max_abs_wall * 10.0) {
    std::cout << "    速度散度贡献占主导，显式更新速度会显著改变右边项" << std::endl;
  } else if (max_abs_wall > max_abs_div * 10.0) {
    std::cout << "    压力边界条件贡献占主导，显式更新速度影响较小" << std::endl;
  } else {
    std::cout << "    两种贡献相当，显式更新速度会改变右边项的平衡" << std::endl;
  }
  
  // 5. 矩阵条件数的影响（虽然A不变，但不同的b可能暴露不同的数值特性）
  std::cout << "\n  5. 矩阵条件数和数值稳定性：" << std::endl;
  std::cout << "    虽然矩阵A不变，但不同的b可能：" << std::endl;
  std::cout << "    - 在矩阵-向量乘法中产生不同的舍入误差累积" << std::endl;
  std::cout << "    - 影响预处理器的效果（预处理器可能对不同的b有不同的表现）" << std::endl;
  std::cout << "    - 改变迭代求解器的收敛路径" << std::endl;
  
  // 6. 对比分析：计算不使用显式更新时的右边项（用于对比）
  std::cout << "\n  6. 对比分析：计算不使用显式更新时的右边项..." << std::endl;
  
  // 临时恢复初始速度，计算不使用显式更新时的右边项
  std::vector<double2> velocity_backup(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
    velocity_backup[i] = fluid_particles.velocity[i];  // 保存当前速度
    fluid_particles.velocity[i] = initial_velocity[i];   // 恢复初始速度（零速度）
  }
  
  // 计算不使用显式更新时的速度散度贡献
  std::vector<double> divergence_contribution_no_explicit(num_fluid_particles, 0.0);
  std::vector<double> wall_pressure_contribution_no_explicit(num_fluid_particles, 0.0);
  std::vector<double> b_total_no_explicit(num_fluid_particles, 0.0);
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    const double2& pos_i = fluid_particles.position[i];
    const double2& vel_i = fluid_particles.velocity[i];  // 现在是初始速度（零）
    
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C1_velocity = corrective_matrices_velocity[i].row(0);
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C2_velocity = corrective_matrices_velocity[i].row(1);
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> M2_pressure = corrective_matrices_pressure[i].row(2);
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> M3_pressure = corrective_matrices_pressure[i].row(3);
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> M2_plus_M3_pressure = M2_pressure + M3_pressure;
    
    double divergence_sum = 0.0;
    double wall_pressure_term = 0.0;
    
    // 遍历流体邻域粒子
    for (int j : fluid_particles.fluid_neighbour_list[i]) {
      const double2& pos_j = fluid_particles.position[j];
      const double2& vel_j = fluid_particles.velocity[j];  // 现在是初始速度（零）
      
      double dx = pos_j.x - pos_i.x;
      double dy = pos_j.y - pos_i.y;
      double dist = ComputeDistance(pos_i, pos_j);
      
      if (dist < 1e-10 || dist > smoothing_radius) {
        continue;
      }
      
      double weight = WeightFunction(dist, smoothing_radius);
      Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis = 
          corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
      
      double dvx = vel_j.x - vel_i.x;
      double dvy = vel_j.y - vel_i.y;
      double C1P = (C1_velocity * basis)(0, 0);
      double C2P = (C2_velocity * basis)(0, 0);
      divergence_sum += weight * (C1P * dvx + C2P * dvy);
    }
    
    // 处理壁面邻域粒子
    for (int j : fluid_particles.solid_neighbour_list[i]) {
      const double2& pos_j = solid_particles.position[j];
      const double2& vel_wall = solid_particles.velocity[j];
      const double2& normal = solid_particles.normal_vector[j];
      
      double dx = pos_j.x - pos_i.x;
      double dy = pos_j.y - pos_i.y;
      double dist = ComputeDistance(pos_i, pos_j);
      
      if (dist < 1e-10 || dist > smoothing_radius) {
        continue;
      }
      
      double weight = WeightFunction(dist, smoothing_radius);
      Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis_velocity = 
          corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
      Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis_pressure = 
          corrective_matrix_calc.ComputeBasisFunctionsForWall(
              dx, dy, normal.x, normal.y, smoothing_radius);
      
      double dvx = vel_wall.x - vel_i.x;
      double dvy = vel_wall.y - vel_i.y;
      double C1P = (C1_velocity * basis_velocity)(0, 0);
      double C2P = (C2_velocity * basis_velocity)(0, 0);
      divergence_sum += weight * (C1P * dvx + C2P * dvy);
      
      double n_dot_g = normal.x * gravity_x + normal.y * gravity_y;
      double wall_pressure_coeff = weight * smoothing_radius * rho * n_dot_g *
          (M2_plus_M3_pressure * basis_pressure)(0, 0);
      wall_pressure_term += wall_pressure_coeff;
    }
    
    divergence_contribution_no_explicit[i] = (1.0 / (smoothing_radius * time_step)) * divergence_sum;
    wall_pressure_contribution_no_explicit[i] = coeff_factor * wall_pressure_term;
    b_total_no_explicit[i] = divergence_contribution_no_explicit[i] + wall_pressure_contribution_no_explicit[i];
  }
  
  // 恢复显式更新后的速度
  for (int i = 0; i < num_fluid_particles; ++i) {
    fluid_particles.velocity[i] = velocity_backup[i];
  }
  
  // 对比分析：统计总体右边项的最小值和最大值
  double b_norm_no_explicit = 0.0;
  double max_abs_b_no_explicit = 0.0;
  double min_abs_b_no_explicit = std::numeric_limits<double>::max();
  double min_abs_b_nonzero_no_explicit = std::numeric_limits<double>::max();  // 非零值的最小值
  int min_b_particle_idx = -1;  // 产生最小值的粒子索引
  int zero_count = 0;  // 零值粒子数量
  int near_wall_count = 0;  // 靠近壁面的粒子数
  int far_wall_count = 0;   // 远离壁面的粒子数
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    b_norm_no_explicit += b_total_no_explicit[i] * b_total_no_explicit[i];
    double abs_b = std::abs(b_total_no_explicit[i]);
    if (abs_b > max_abs_b_no_explicit) max_abs_b_no_explicit = abs_b;
    
    // 统计所有值的最小值（包括0）
    if (abs_b < min_abs_b_no_explicit) {
      min_abs_b_no_explicit = abs_b;
    }
    
    // 统计非零值的最小值
    if (abs_b < 1e-10) {
      zero_count++;
    } else {
      if (abs_b < min_abs_b_nonzero_no_explicit) {
        min_abs_b_nonzero_no_explicit = abs_b;
        min_b_particle_idx = i;
      }
    }
    
    // 统计靠近壁面和远离壁面的粒子
    if (!fluid_particles.solid_neighbour_list[i].empty()) {
      near_wall_count++;
    } else {
      far_wall_count++;
      // 检查远离壁面的粒子的右边项是否真的为0
      if (abs_b > 1e-10) {
        std::cerr << "    ⚠️  警告：远离壁面的粒子 " << i << " 的右边项不为0: " 
                  << std::scientific << abs_b << std::fixed << std::endl;
        std::cerr << "      速度散度贡献: " << std::scientific 
                  << divergence_contribution_no_explicit[i] << std::fixed << std::endl;
        std::cerr << "      压力边界条件贡献: " << std::scientific 
                  << wall_pressure_contribution_no_explicit[i] << std::fixed << std::endl;
      }
    }
  }
  b_norm_no_explicit = std::sqrt(b_norm_no_explicit);
  // 计算范围比值（如果最小值为0，使用非零最小值）
  double min_for_range = (min_abs_b_no_explicit < 1e-10) ? 
                         ((min_abs_b_nonzero_no_explicit < std::numeric_limits<double>::max()) ? 
                          min_abs_b_nonzero_no_explicit : min_abs_b_no_explicit) : 
                         min_abs_b_no_explicit;
  double b_range_no_explicit = (min_for_range > 1e-20) ? 
                               (max_abs_b_no_explicit / min_for_range) : 0.0;
  
  std::cout << "\n  对比结果：" << std::endl;
  std::cout << "    不使用显式更新：" << std::endl;
  std::cout << "      右边项b的范数: " << std::scientific << std::setprecision(3) 
            << b_norm_no_explicit << std::fixed << std::endl;
  std::cout << "      右边项b的最小值（包括0）: " << std::scientific << std::setprecision(3) 
            << min_abs_b_no_explicit << std::fixed << std::endl;
  std::cout << "      右边项b的最小值（非零）: " << std::scientific << std::setprecision(3) 
            << ((min_abs_b_nonzero_no_explicit < std::numeric_limits<double>::max()) ? min_abs_b_nonzero_no_explicit : 0.0) 
            << std::fixed << std::endl;
  std::cout << "      右边项b的最大值: " << std::scientific << std::setprecision(3) 
            << max_abs_b_no_explicit << std::fixed << std::endl;
  std::cout << "      范围比值: " << std::scientific << std::setprecision(3) 
            << b_range_no_explicit << std::fixed << std::endl;
  std::cout << "      零值粒子数: " << zero_count << " / " << num_fluid_particles << std::endl;
  std::cout << "      靠近壁面的粒子数: " << near_wall_count << std::endl;
  std::cout << "      远离壁面的粒子数: " << far_wall_count << std::endl;
  
  if (zero_count > 0) {
    std::cout << "      ✓ 有 " << zero_count << " 个粒子的右边项为0（远离壁面的粒子）" << std::endl;
  }
  if (min_abs_b_no_explicit < 1e-10) {
    std::cout << "      ✓ 最小值确实为0（来自远离壁面的粒子）" << std::endl;
  } else {
    std::cout << "      ⚠️  最小值不为0，可能所有粒子都靠近壁面，或者有数值误差" << std::endl;
  }
  
  // 详细分析产生最小值的粒子
  if (min_b_particle_idx >= 0) {
    std::cout << "\n      最小值来源分析（粒子 " << min_b_particle_idx << "）：" << std::endl;
    std::cout << "        粒子位置: (" << std::fixed << std::setprecision(6) 
              << fluid_particles.position[min_b_particle_idx].x << ", " 
              << fluid_particles.position[min_b_particle_idx].y << ")" << std::endl;
    std::cout << "        速度散度贡献: " << std::scientific << std::setprecision(6) 
              << divergence_contribution_no_explicit[min_b_particle_idx] << std::fixed << std::endl;
    std::cout << "        压力边界条件贡献: " << std::scientific << std::setprecision(6) 
              << wall_pressure_contribution_no_explicit[min_b_particle_idx] << std::fixed << std::endl;
    std::cout << "        右边项b总值: " << std::scientific << std::setprecision(6) 
              << b_total_no_explicit[min_b_particle_idx] << std::fixed << std::endl;
    std::cout << "        壁面邻域粒子数: " << fluid_particles.solid_neighbour_list[min_b_particle_idx].size() << std::endl;
    
    // 计算参数
    std::cout << "\n        计算参数：" << std::endl;
    std::cout << "          coeff_factor = 2.0 / (smoothing_radius^2 * rho) = " 
              << std::scientific << std::setprecision(6) << coeff_factor << std::fixed << std::endl;
    std::cout << "          smoothing_radius = " << smoothing_radius << " m" << std::endl;
    std::cout << "          rho = " << rho << " kg/m³" << std::endl;
    std::cout << "          gravity_y = " << gravity_y << " m/s²" << std::endl;
    
    // 如果有壁面邻域粒子，显示详细计算
    if (!fluid_particles.solid_neighbour_list[min_b_particle_idx].empty()) {
      std::cout << "\n        壁面压力边界条件贡献的组成：" << std::endl;
      std::cout << "          wall_pressure_term = sum(-weight * rho * n_dot_g * (C3+C4) * basis_pressure)" << std::endl;
      std::cout << "          wall_pressure_contribution = coeff_factor * wall_pressure_term" << std::endl;
      std::cout << "          = " << std::scientific << std::setprecision(6) << coeff_factor << " * " 
                << (wall_pressure_contribution_no_explicit[min_b_particle_idx] / coeff_factor) << std::fixed << std::endl;
      std::cout << "          = " << std::scientific << std::setprecision(6) 
                << wall_pressure_contribution_no_explicit[min_b_particle_idx] << std::fixed << std::endl;
    }
  }
  
  std::cout << "    使用显式更新：" << std::endl;
  std::cout << "      右边项b的范数: " << std::scientific << std::setprecision(3) 
            << b_norm_actual << std::fixed << std::endl;
  std::cout << "      右边项b的最小值: " << std::scientific << std::setprecision(3) 
            << ((min_abs_b < std::numeric_limits<double>::max()) ? min_abs_b : 0.0) 
            << std::fixed << std::endl;
  std::cout << "      右边项b的最大值: " << std::scientific << std::setprecision(3) 
            << max_abs_b << std::fixed << std::endl;
  std::cout << "      范围比值: " << std::scientific << std::setprecision(3) 
            << b_range << std::fixed << std::endl;
  
  std::cout << "\n  关键发现：" << std::endl;
  if (b_norm_no_explicit < b_norm_actual * 0.1) {
    std::cout << "    ✓ 不使用显式更新时，右边项b的范数显著更小（" 
              << (b_norm_actual / b_norm_no_explicit) << "倍）" << std::endl;
    std::cout << "    ✓ 即使范围比值很大，但整体范数小，相对残差更容易达到" << std::endl;
    std::cout << "    ✓ 目标绝对残差: " << std::scientific << std::setprecision(3) 
              << (target_tolerance * b_norm_no_explicit) << " (不使用显式更新)" << std::endl;
    std::cout << "    ✓ 目标绝对残差: " << std::scientific << std::setprecision(3) 
              << (target_tolerance * b_norm_actual) << " (使用显式更新)" << std::endl;
    std::cout << "    ✓ 使用显式更新需要达到的绝对残差是 " 
              << (b_norm_actual / b_norm_no_explicit) << " 倍！" << std::endl;
  }
  
  if (b_range_no_explicit > 1e10 && b_norm_no_explicit < 1.0) {
    std::cout << "\n    ✓ 重要发现：即使范围比值很大，如果右边项范数较小，仍能收敛！" << std::endl;
    std::cout << "      原因：" << std::endl;
    std::cout << "      1. 相对残差 = ||r|| / ||b||，如果||b||小，即使||r||较大，相对残差也可能满足" << std::endl;
    std::cout << "      2. 数值精度：如果整体范数小，矩阵-向量乘法中的舍入误差也相对较小" << std::endl;
    std::cout << "      3. 迭代求解器的收敛性主要取决于相对残差，而不是绝对残差" << std::endl;
    std::cout << "      4. 范围比值大但范数小，说明大部分元素都很小，只有少数元素较大" << std::endl;
    std::cout << "      5. 这种分布可能对迭代求解器更友好（主要能量集中在少数元素）" << std::endl;
  }
  
  // ========== 分析结束 ==========
  
  // ========== 提前输出VTK文件（求解前）==========
  std::cout << "\n[输出] 提前输出VTK文件（求解前）..." << std::endl;
  std::string result_filename_before = "output/ppe_result_hydrostatic_before_solve.vtk";
  FileOperator file_op_before;
  
  // 输出流体粒子基础信息
  if (file_op_before.writeVTKBase(result_filename_before, fluid_particles)) {
    // 追加显式更新后的速度向量
    if (file_op_before.appendVTKVector(result_filename_before, "velocity_explicit", velocity_explicit)) {
      std::cout << "  显式更新后的速度向量已追加" << std::endl;
    }
    
    // 追加初始速度向量
    if (file_op_before.appendVTKVector(result_filename_before, "velocity_initial", initial_velocity)) {
      std::cout << "  初始速度向量已追加" << std::endl;
    }
    
    // 追加粘性力加速度向量
    if (file_op_before.appendVTKVector(result_filename_before, "viscous_acceleration", viscous_acceleration)) {
      std::cout << "  粘性力加速度向量已追加" << std::endl;
    }
    
    // 追加右边项b（从PETSc向量提取）
    std::vector<double> b_values_double(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
      b_values_double[i] = static_cast<double>(values[i]);
    }
    if (file_op_before.appendVTKScalar(result_filename_before, "b_total", b_values_double)) {
      std::cout << "  右边项b（总）已追加" << std::endl;
    }
    
    // 追加速度散度贡献
    if (file_op_before.appendVTKScalar(result_filename_before, "b_divergence_contribution", divergence_contribution)) {
      std::cout << "  速度散度贡献已追加" << std::endl;
    }
    
    // 追加压力边界条件贡献
    if (file_op_before.appendVTKScalar(result_filename_before, "b_wall_pressure_contribution", wall_pressure_contribution)) {
      std::cout << "  压力边界条件贡献已追加" << std::endl;
    }
    
    // 追加重新计算的b（用于验证）
    if (file_op_before.appendVTKScalar(result_filename_before, "b_total_recomputed", b_total_recomputed)) {
      std::cout << "  重新计算的b（验证）已追加" << std::endl;
    }
    
    // 追加对角线元素（用于检查压力拉普拉斯算子离散）
    if (file_op_before.appendVTKScalar(result_filename_before, "diagonal_elements", diagonal_elements)) {
      std::cout << "  对角线元素已追加" << std::endl;
    }
    
    // 追加自由面类型信息
    std::vector<int> surface_type_values(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
      surface_type_values[i] = static_cast<int>(fluid_particles.surface_type[i]);
    }
    if (file_op_before.appendVTKScalar(result_filename_before, "surface_type", surface_type_values)) {
      std::cout << "  自由面类型已追加" << std::endl;
    }
    
    std::cout << "  求解前的VTK文件已保存到: " << result_filename_before << std::endl;
    } else {
    std::cerr << "  警告：无法创建求解前的VTK文件" << std::endl;
  }
  // ========== 提前输出结束 ==========
  
  // ========== 求解PPE方程（直接求解 A·p = b）==========
  std::cout << "\n求解PPE方程（直接求解 A·p = b）..." << std::endl;
  PPESolver::SolverConfig solver_config;
  solver_config.solver_type = PPESolver::SolverType::BICGSTAB;  // 推荐用于非对称矩阵
  solver_config.max_iterations = 10000;  // 最大迭代次数
  solver_config.tolerance = 1e-6;  // 容差
  solver_config.force_iterative = true;  // 强制使用迭代求解器
  solver_config.is_symmetric_positive_definite = false;
  solver_config.restart = 30;
  
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
  
  // 使用LSMPS计算压力梯度（第二类边界条件，按diffuse_derivative_constraint.md）
  std::cout << "\n使用LSMPS计算压力梯度..." << std::endl;
  std::vector<double2> pressure_gradient(num_fluid_particles);
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    const double2& pos_i = fluid_particles.position[i];
    double p_i = pressure_values[i];
    
    // 提取 [M_{i,0}, M_{i,1}]
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> M0 = corrective_matrices_pressure[i].row(0);  // x方向
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> M1 = corrective_matrices_pressure[i].row(1);  // y方向
    
    double sum_x = 0.0;
    double sum_y = 0.0;
    
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
      
      double dp = p_j - p_i;
      double weight = WeightFunction(dist, smoothing_radius);
      
      // 计算基函数（标准基函数）
      Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis = 
          corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
      
      double M0P = (M0 * basis)(0, 0);
      double M1P = (M1 * basis)(0, 0);
      sum_x += weight * dp * M0P;
      sum_y += weight * dp * M1P;
    }
    
    // 处理壁面邻域粒子（第二类边界条件）
    // 壁面处压力梯度的法向分量：dp/dn = ρ (n · g)
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
      
      // 归一化法向量
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
  
  // 计算梯度大小
  std::vector<double> grad_magnitude(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
    grad_magnitude[i] = std::sqrt(pressure_gradient[i].x * pressure_gradient[i].x + 
                                   pressure_gradient[i].y * pressure_gradient[i].y);
  }
  
  // 输出所有流体粒子信息到同一个VTK文件
  std::cout << "\n输出所有流体粒子信息到VTK文件..." << std::endl;
  std::string result_filename = "output/ppe_result_hydrostatic.vtk";
  // 复用之前创建的file_op_before，或者创建新的file_op
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
    
    // 追加显式更新后的速度向量（vector格式）
    if (file_op.appendVTKVector(result_filename, "velocity_explicit", velocity_explicit)) {
      std::cout << "  显式更新后的速度向量已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加显式更新后的速度向量到VTK文件" << std::endl;
    }
    
    // 追加初始速度向量（vector格式）
    if (file_op.appendVTKVector(result_filename, "velocity_initial", initial_velocity)) {
      std::cout << "  初始速度向量已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加初始速度向量到VTK文件" << std::endl;
    }
    
    // 追加粘性力加速度向量（vector格式）
    if (file_op.appendVTKVector(result_filename, "viscous_acceleration", viscous_acceleration)) {
      std::cout << "  粘性力加速度向量已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加粘性力加速度向量到VTK文件" << std::endl;
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
    
    // 追加调试信息（原始矩阵A的对角线元素和右边项b）
    if (matrix_builder.WriteDebugInfoToVTK(fluid_particles, A_petsc, b_petsc, result_filename)) {
      std::cout << "  调试信息（原始矩阵A的对角线元素、右边项b等）已追加" << std::endl;
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
  
  // 清理PETSc对象（按照创建顺序的逆序销毁）
  if (p_petsc != NULL) {
    VecDestroy(&p_petsc);
  }
  if (b_petsc != NULL) {
    VecDestroy(&b_petsc);
  }
  if (A_petsc != NULL) {
    MatDestroy(&A_petsc);
  }
  
  return 0;
}
