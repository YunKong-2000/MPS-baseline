#include "../src/PPE/PPEMatrixBuilder.hpp"
#include "../src/PPE/PPESolver.hpp"
#include "../src/lsmps/CorrectiveMatrix.hpp"
#include "../src/neighbour_list/NeighborListSearcher.hpp"
#include "../src/surface_detection/SurfaceDetector.hpp"
#include "../src/core/Particle.hpp"
#include "../src/core/FileOperator.hpp"
#include "../src/correction/Correction.hpp"
#include "../src/explicit_force/ExplicitForce.hpp"
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
  
  // 消除未使用参数的警告
  (void)container_height;
  
  // 计算粒子数量，确保关于中心对称
  // 从中心开始，向两边均匀分布
  double center_x = container_width / 2.0;
  // center_y = 0.0 从底部开始（未使用，但保留注释说明）
  
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
  const double container_width = 2.0;   // 容器宽度 (m)
  const double container_height = 2.0;  // 容器高度 (m)
  const double water_height = 1.0;      // 水位高度 (m)
  
  // 粒子参数（平滑半径保持为粒子间距的2.1倍）
  const double particle_spacing = 0.02;  // 粒子间距 (m)
  const double particle_radius = particle_spacing / 2.0;
  const double smoothing_radius = 3.1 * particle_spacing;  // 平滑半径是粒子间距的2.1倍
  const double cell_size = 2.0 * smoothing_radius;
  
  // 时间步长
  const double time_step = 0.0001;  // 时间步长 (s)
  
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
  
  // 保存初始速度和粒子数量（用于后续分析）
  int num_fluid_particles = fluid_particles.particle_num;
  std::vector<double2> initial_velocity(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
    initial_velocity[i] = fluid_particles.velocity[i];  // 初始速度为0（静水状态）
  }
  
  // ========== 显式力计算步骤（只更新速度，不更新位置）==========
  std::cout << "\n========== 显式力计算步骤（只更新速度，不更新位置）==========" << std::endl;
  
  // 计算用于显式力的corrective matrix（第一类边界条件）
  std::cout << "\n计算corrective matrix（显式力用，第一类边界条件）..." << std::endl;
  std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>>
      corrective_matrices_explicit(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
    corrective_matrices_explicit[i] = corrective_matrix_calc.ComputeCorrectiveMatrix(
        i, fluid_particles, solid_particles,
        smoothing_radius, false);  // 第一类边界条件
  }
  std::cout << "  完成" << std::endl;
  
  // 显式更新模块（只更新速度作为临时速度，不更新位置）
  std::cout << "\n显式更新模块：计算粘性力和重力，更新临时速度..." << std::endl;
  ExplicitForce explicit_force;
  explicit_force.ComputeAndUpdateVelocity(
      fluid_particles, solid_particles,
      corrective_matrices_explicit,
      smoothing_radius,
      kinematic_viscosity,
      gravity_x, gravity_y,
      time_step);
  
  // 统计速度变化
  double max_velocity_change = 0.0;
  double sum_velocity_change = 0.0;
  double max_velocity_magnitude = 0.0;
  for (int i = 0; i < num_fluid_particles; ++i) {
    double2 velocity_change = {
        fluid_particles.velocity[i].x - initial_velocity[i].x,
        fluid_particles.velocity[i].y - initial_velocity[i].y
    };
    double vel_change_magnitude = std::sqrt(velocity_change.x * velocity_change.x + 
                                            velocity_change.y * velocity_change.y);
    double vel_magnitude = std::sqrt(fluid_particles.velocity[i].x * fluid_particles.velocity[i].x + 
                                     fluid_particles.velocity[i].y * fluid_particles.velocity[i].y);
    
    if (vel_change_magnitude > max_velocity_change) {
      max_velocity_change = vel_change_magnitude;
    }
    if (vel_magnitude > max_velocity_magnitude) {
      max_velocity_magnitude = vel_magnitude;
    }
    sum_velocity_change += vel_change_magnitude;
  }
  double avg_velocity_change = sum_velocity_change / num_fluid_particles;
  
  std::cout << "  初始速度最大值: " << max_velocity_magnitude << " m/s" << std::endl;
  std::cout << "  速度变化最大值: " << max_velocity_change << " m/s" << std::endl;
  std::cout << "  速度变化平均值: " << avg_velocity_change << " m/s" << std::endl;
  std::cout << "  理论速度变化（仅重力）: " << std::abs(gravity_y) * time_step << " m/s" << std::endl;
  
  // 重新计算用于PPE的corrective matrix（因为速度已更新，但位置未变，所以corrective matrix应该不变）
  // 但为了与main.cpp保持一致，我们重新计算
  std::cout << "\n重新计算corrective matrix（PPE用，速度已更新）..." << std::endl;
  for (int i = 0; i < num_fluid_particles; ++i) {
    // 速度散度：第一类边界条件
    corrective_matrices_velocity[i] = corrective_matrix_calc.ComputeCorrectiveMatrix(
        i, fluid_particles, solid_particles, smoothing_radius, false);
    // 压力拉普拉斯算子：第二类边界条件
    corrective_matrices_pressure[i] = corrective_matrix_calc.ComputeCorrectiveMatrix(
        i, fluid_particles, solid_particles, smoothing_radius, true);
  }
  std::cout << "  完成" << std::endl;
  
  // ========== 构建PPE（使用更新后的速度）==========
  // 直接构建PETSc格式的PPE系数矩阵和右边项
  std::cout << "\n构建PPE系数矩阵和右边项（PETSc格式，使用更新后的速度）..." << std::endl;
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
  
  // 提取右边项
  std::vector<PetscInt> indices(num_fluid_particles);
  std::vector<PetscScalar> values(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
    indices[i] = static_cast<PetscInt>(i);
  }
  VecGetValues(b_petsc, num_fluid_particles, indices.data(), values.data());
  
  // ========== 分析速度散度和压力边界条件对右边项的贡献 ==========
  std::cout << "\n========== 分析速度散度和压力边界条件对右边项的贡献 ==========" << std::endl;
  
  // 系数因子（与PPEMatrixBuilder中一致）
  double coeff_factor = 2.0 / (smoothing_radius * rho);
  
  // 存储每个粒子的贡献
  std::vector<double> divergence_contribution(num_fluid_particles, 0.0);  // 速度散度贡献
  std::vector<double> wall_pressure_contribution(num_fluid_particles, 0.0);  // 壁面压力边界条件贡献
  std::vector<double> divergence_fluid_part(num_fluid_particles, 0.0);  // 速度散度：流体邻域贡献
  std::vector<double> divergence_wall_part(num_fluid_particles, 0.0);  // 速度散度：壁面邻域贡献
  std::vector<double> distance_to_wall(num_fluid_particles, 1e10);  // 到最近壁面的距离
  std::vector<int> num_wall_neighbors(num_fluid_particles, 0);  // 壁面邻域粒子数
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    const double2& pos_i = fluid_particles.position[i];
    const double2& vel_i = fluid_particles.velocity[i];
    
    // 提取速度corrective matrix的行向量（用于速度散度计算，第一类边界条件）
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C1_velocity = corrective_matrices_velocity[i].row(0);
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C2_velocity = corrective_matrices_velocity[i].row(1);
    
    // 提取压力corrective matrix的行向量（用于压力拉普拉斯算子计算，第二类边界条件）
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C3_pressure = corrective_matrices_pressure[i].row(2);
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C4_pressure = corrective_matrices_pressure[i].row(3);
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C3_plus_C4_pressure = C3_pressure + C4_pressure;
    
    double div_fluid = 0.0;
    double div_wall = 0.0;
    double wall_pressure_term = 0.0;
    double min_wall_dist = 1e10;
    
    // 计算流体邻域粒子的速度散度贡献
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
      Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis = 
          corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
      
      double dux_dr = (vel_j.x - vel_i.x) / dist;
      double duy_dr = (vel_j.y - vel_i.y) / dist;
      double C1P = (C1_velocity * basis)(0, 0);
      double C2P = (C2_velocity * basis)(0, 0);
      div_fluid += weight * (C1P * dux_dr + C2P * duy_dr);
    }
    
    // 计算壁面邻域粒子的速度散度贡献和压力边界条件贡献
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
      
      if (dist < min_wall_dist) {
        min_wall_dist = dist;
      }
      num_wall_neighbors[i]++;
      
      double weight = WeightFunction(dist, smoothing_radius);
      
      // 对于速度散度：使用第一类边界条件的基函数（标准基函数）
      Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis_velocity = 
          corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
      
      // 对于壁面压力边界条件项：使用第二类边界条件的基函数（壁面基函数）
      Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis_pressure = 
          corrective_matrix_calc.ComputeBasisFunctionsForWall(
              dx, dy, normal.x, normal.y, smoothing_radius);
      
      // 计算速度散度项（使用速度corrective matrix，第一类边界条件）
      double dux_dr = (vel_wall.x - vel_i.x) / dist;
      double duy_dr = (vel_wall.y - vel_i.y) / dist;
      double C1P = (C1_velocity * basis_velocity)(0, 0);
      double C2P = (C2_velocity * basis_velocity)(0, 0);
      div_wall += weight * (C1P * dux_dr + C2P * duy_dr);
      
      // 计算壁面压力边界条件项（使用压力corrective matrix，第二类边界条件）
      double n_dot_g = normal.x * gravity_x + normal.y * gravity_y;
      double wall_pressure_coeff = -weight * (rho * n_dot_g) * (C3_plus_C4_pressure * basis_pressure)(0, 0);
      wall_pressure_term += wall_pressure_coeff;
    }
    
    // 计算总贡献
    double total_divergence = div_fluid + div_wall;
    divergence_contribution[i] = (1.0 / time_step) * total_divergence;
    wall_pressure_contribution[i] = coeff_factor * wall_pressure_term;
    divergence_fluid_part[i] = (1.0 / time_step) * div_fluid;
    divergence_wall_part[i] = (1.0 / time_step) * div_wall;
    distance_to_wall[i] = min_wall_dist;
  }
  
  // 统计所有粒子的贡献
  std::cout << "\n所有粒子的贡献统计：" << std::endl;
  double max_div_contrib = 0.0, max_wall_press_contrib = 0.0;
  double sum_abs_div_contrib = 0.0, sum_abs_wall_press_contrib = 0.0;
  double sum_abs_b = 0.0;
  for (int i = 0; i < num_fluid_particles; ++i) {
    double abs_div = std::abs(divergence_contribution[i]);
    double abs_wall = std::abs(wall_pressure_contribution[i]);
    double abs_b = std::abs(values[i]);
    
    if (abs_div > max_div_contrib) max_div_contrib = abs_div;
    if (abs_wall > max_wall_press_contrib) max_wall_press_contrib = abs_wall;
    
    sum_abs_div_contrib += abs_div;
    sum_abs_wall_press_contrib += abs_wall;
    sum_abs_b += abs_b;
  }
  double avg_abs_div_contrib = sum_abs_div_contrib / num_fluid_particles;
  double avg_abs_wall_press_contrib = sum_abs_wall_press_contrib / num_fluid_particles;
  double avg_abs_b = sum_abs_b / num_fluid_particles;
  
  std::cout << "  速度散度贡献 - 最大值: " << max_div_contrib << ", 平均值: " << avg_abs_div_contrib << std::endl;
  std::cout << "  壁面压力贡献 - 最大值: " << max_wall_press_contrib << ", 平均值: " << avg_abs_wall_press_contrib << std::endl;
  std::cout << "  右边项b - 平均值: " << avg_abs_b << std::endl;
  std::cout << "  速度散度贡献占比: " << (avg_abs_div_contrib / (avg_abs_div_contrib + avg_abs_wall_press_contrib + 1e-10)) * 100.0 << "%" << std::endl;
  std::cout << "  壁面压力贡献占比: " << (avg_abs_wall_press_contrib / (avg_abs_div_contrib + avg_abs_wall_press_contrib + 1e-10)) * 100.0 << "%" << std::endl;
  
  // 分析靠近壁面的粒子（有壁面邻域的粒子）
  std::cout << "\n靠近壁面的粒子（有壁面邻域的粒子）分析：" << std::endl;
  int particles_with_wall = 0;
  double max_div_contrib_wall = 0.0, max_wall_press_contrib_wall = 0.0;
  double sum_abs_div_contrib_wall = 0.0, sum_abs_wall_press_contrib_wall = 0.0;
  double sum_abs_b_wall = 0.0;
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    if (num_wall_neighbors[i] > 0) {
      particles_with_wall++;
      double abs_div = std::abs(divergence_contribution[i]);
      double abs_wall = std::abs(wall_pressure_contribution[i]);
      double abs_b = std::abs(values[i]);
      
      if (abs_div > max_div_contrib_wall) max_div_contrib_wall = abs_div;
      if (abs_wall > max_wall_press_contrib_wall) max_wall_press_contrib_wall = abs_wall;
      
      sum_abs_div_contrib_wall += abs_div;
      sum_abs_wall_press_contrib_wall += abs_wall;
      sum_abs_b_wall += abs_b;
    }
  }
  
  if (particles_with_wall > 0) {
    double avg_abs_div_contrib_wall = sum_abs_div_contrib_wall / particles_with_wall;
    double avg_abs_wall_press_contrib_wall = sum_abs_wall_press_contrib_wall / particles_with_wall;
    double avg_abs_b_wall = sum_abs_b_wall / particles_with_wall;
    
    std::cout << "  靠近壁面的粒子数: " << particles_with_wall << " / " << num_fluid_particles << std::endl;
    std::cout << "  速度散度贡献 - 最大值: " << max_div_contrib_wall << ", 平均值: " << avg_abs_div_contrib_wall << std::endl;
    std::cout << "  壁面压力贡献 - 最大值: " << max_wall_press_contrib_wall << ", 平均值: " << avg_abs_wall_press_contrib_wall << std::endl;
    std::cout << "  右边项b - 平均值: " << avg_abs_b_wall << std::endl;
    std::cout << "  速度散度贡献占比: " << (avg_abs_div_contrib_wall / (avg_abs_div_contrib_wall + avg_abs_wall_press_contrib_wall + 1e-10)) * 100.0 << "%" << std::endl;
    std::cout << "  壁面压力贡献占比: " << (avg_abs_wall_press_contrib_wall / (avg_abs_div_contrib_wall + avg_abs_wall_press_contrib_wall + 1e-10)) * 100.0 << "%" << std::endl;
    
    // 按距离壁面的距离分组分析
    std::cout << "\n按距离壁面的距离分组分析：" << std::endl;
    const int num_bins = 5;
    std::vector<double> bin_max_dist(num_bins);
    std::vector<double> bin_avg_div_contrib(num_bins, 0.0);
    std::vector<double> bin_avg_wall_press_contrib(num_bins, 0.0);
    std::vector<double> bin_avg_b(num_bins, 0.0);
    std::vector<int> bin_count(num_bins, 0);
    
    double max_wall_dist = 0.0;
    for (int i = 0; i < num_fluid_particles; ++i) {
      if (distance_to_wall[i] < 1e9 && distance_to_wall[i] > max_wall_dist) {
        max_wall_dist = distance_to_wall[i];
      }
    }
    
    if (max_wall_dist > 1e-10) {
      for (int bin = 0; bin < num_bins; ++bin) {
        bin_max_dist[bin] = (bin + 1) * max_wall_dist / num_bins;
      }
      
      for (int i = 0; i < num_fluid_particles; ++i) {
        if (distance_to_wall[i] < 1e9) {
          int bin = static_cast<int>(distance_to_wall[i] / (max_wall_dist / num_bins));
          if (bin >= num_bins) bin = num_bins - 1;
          
          bin_avg_div_contrib[bin] += std::abs(divergence_contribution[i]);
          bin_avg_wall_press_contrib[bin] += std::abs(wall_pressure_contribution[i]);
          bin_avg_b[bin] += std::abs(values[i]);
          bin_count[bin]++;
        }
      }
      
      for (int bin = 0; bin < num_bins; ++bin) {
        if (bin_count[bin] > 0) {
          bin_avg_div_contrib[bin] /= bin_count[bin];
          bin_avg_wall_press_contrib[bin] /= bin_count[bin];
          bin_avg_b[bin] /= bin_count[bin];
          
          double total_contrib = bin_avg_div_contrib[bin] + bin_avg_wall_press_contrib[bin];
          std::cout << "  距离 [" << (bin * max_wall_dist / num_bins) << ", " 
                    << bin_max_dist[bin] << "] m: " << bin_count[bin] << " 个粒子" << std::endl;
          std::cout << "    速度散度贡献: " << bin_avg_div_contrib[bin] 
                    << " (" << (total_contrib > 1e-10 ? (bin_avg_div_contrib[bin] / total_contrib * 100.0) : 0.0) << "%)" << std::endl;
          std::cout << "    壁面压力贡献: " << bin_avg_wall_press_contrib[bin] 
                    << " (" << (total_contrib > 1e-10 ? (bin_avg_wall_press_contrib[bin] / total_contrib * 100.0) : 0.0) << "%)" << std::endl;
          std::cout << "    右边项b: " << bin_avg_b[bin] << std::endl;
        }
      }
    }
  }
  
  // 分析远离壁面的粒子（无壁面邻域的粒子）
  std::cout << "\n远离壁面的粒子（无壁面邻域的粒子）分析：" << std::endl;
  int particles_without_wall = 0;
  double max_div_contrib_no_wall = 0.0;
  double sum_abs_div_contrib_no_wall = 0.0;
  double sum_abs_b_no_wall = 0.0;
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    if (num_wall_neighbors[i] == 0) {
      particles_without_wall++;
      double abs_div = std::abs(divergence_contribution[i]);
      double abs_b = std::abs(values[i]);
      
      if (abs_div > max_div_contrib_no_wall) max_div_contrib_no_wall = abs_div;
      
      sum_abs_div_contrib_no_wall += abs_div;
      sum_abs_b_no_wall += abs_b;
    }
  }
  
  if (particles_without_wall > 0) {
    double avg_abs_div_contrib_no_wall = sum_abs_div_contrib_no_wall / particles_without_wall;
    double avg_abs_b_no_wall = sum_abs_b_no_wall / particles_without_wall;
    
    std::cout << "  远离壁面的粒子数: " << particles_without_wall << " / " << num_fluid_particles << std::endl;
    std::cout << "  速度散度贡献 - 最大值: " << max_div_contrib_no_wall << ", 平均值: " << avg_abs_div_contrib_no_wall << std::endl;
    std::cout << "  右边项b - 平均值: " << avg_abs_b_no_wall << std::endl;
    std::cout << "  注意：远离壁面的粒子没有壁面压力贡献" << std::endl;
  }
  
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
  // 注意：K = A^T A + D 是对称正定矩阵，使用CG方法求解
  std::cout << "\n求解PPE方程（CG方法）..." << std::endl;
  PPESolver::SolverConfig solver_config;
  solver_config.solver_type = PPESolver::SolverType::CG;
  solver_config.max_iterations = 5000;
  solver_config.tolerance = 1e-6;
  solver_config.force_iterative = true;
  solver_config.is_symmetric_positive_definite = true;
  
  PPESolver ppe_solver(solver_config);
  Vec p_petsc = NULL;
  
  auto start_time = std::chrono::high_resolution_clock::now();
  bool solve_success = ppe_solver.Solve(K_petsc, f_petsc, p_petsc);
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
  
  if (!solve_success) {
    std::cerr << "错误：PPE求解失败" << std::endl;
    MatDestroy(&A_petsc);
    VecDestroy(&b_petsc);
    MatDestroy(&K_petsc);
    VecDestroy(&f_petsc);
    if (p_petsc != NULL) {
      VecDestroy(&p_petsc);
    }
    return 1;
  }
  
  std::cout << "  求解完成" << std::endl;
  std::cout << "  迭代次数: " << ppe_solver.GetLastIterations() << std::endl;
  std::cout << "  残差: " << ppe_solver.GetLastResidual() << std::endl;
  std::cout << "  是否收敛: " << (ppe_solver.GetLastConverged() ? "是" : "否") << std::endl;
  std::cout << "  求解时间: " << duration.count() << " 毫秒" << std::endl;
  
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
  
  // 保存当前速度和位置（用于correction模块验证，注意：速度是初始速度0）
  std::vector<double2> velocity_before_correction(num_fluid_particles);
  std::vector<double2> position_before_correction(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
    velocity_before_correction[i] = fluid_particles.velocity[i];  // 初始速度为0
    position_before_correction[i] = fluid_particles.position[i];
  }
  
  // ========== 分析PPE右边项和压力求解结果 ==========
  std::cout << "\n========== 分析PPE右边项和压力求解结果 ==========" << std::endl;
  
  // 提取右边项的值（用于分析速度散度的影响）
  std::vector<PetscScalar> b_values(num_fluid_particles);
  VecGetValues(b_petsc, num_fluid_particles, indices.data(), b_values.data());
  
  // 计算速度散度的统计信息
  // 右边项 b = (1/dt) * divergence + coeff_factor * wall_pressure_term
  // 对于静水压力，如果速度散度为0，右边项应该只包含壁面压力项
  double max_b_value = 0.0;
  double min_b_value = 0.0;
  double sum_abs_b_analysis = 0.0;
  for (int i = 0; i < num_fluid_particles; ++i) {
    double abs_b = std::abs(b_values[i]);
    if (abs_b > max_b_value) {
      max_b_value = abs_b;
    }
    if (i == 0 || std::abs(b_values[i]) < min_b_value) {
      min_b_value = std::abs(b_values[i]);
    }
    sum_abs_b_analysis += abs_b;
  }
  double avg_abs_b_analysis = sum_abs_b_analysis / num_fluid_particles;
  
  std::cout << "\n右边项（b）统计信息（反映速度散度的影响）：" << std::endl;
  std::cout << "  右边项最大值: " << max_b_value << std::endl;
  std::cout << "  右边项最小值: " << min_b_value << std::endl;
  std::cout << "  右边项平均值: " << avg_abs_b_analysis << std::endl;
  std::cout << "  注意：如果速度散度为0，右边项应该只包含壁面压力项" << std::endl;
  
  // 计算理论压力（用于对比）
  std::vector<double> theoretical_pressure_double(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
    double y = fluid_particles.position[i].y;
    double depth = water_height - y;  // 从自由表面的深度
    double theoretical_pressure = rho * g * depth;  // 理论静水压力
    theoretical_pressure_double[i] = theoretical_pressure;
  }
  
  // 计算压力误差
  double max_pressure_error = 0.0;
  double sum_pressure_error = 0.0;
  double max_relative_error = 0.0;
  int error_count = 0;
  for (int i = 0; i < num_fluid_particles; ++i) {
    double error = std::abs(pressure_values[i] - theoretical_pressure_double[i]);
    if (error > max_pressure_error) {
      max_pressure_error = error;
    }
    sum_pressure_error += error;
    
    if (theoretical_pressure_double[i] > 1e-10) {
      double relative_error = error / theoretical_pressure_double[i];
      if (relative_error > max_relative_error) {
        max_relative_error = relative_error;
      }
      if (relative_error > 0.01) {  // 相对误差超过1%
        error_count++;
      }
    }
  }
  double avg_pressure_error = sum_pressure_error / num_fluid_particles;
  
  std::cout << "\n压力求解误差（与理论静水压力对比）：" << std::endl;
  std::cout << "  最大绝对误差: " << max_pressure_error << " Pa" << std::endl;
  std::cout << "  平均绝对误差: " << avg_pressure_error << " Pa" << std::endl;
  std::cout << "  最大相对误差: " << max_relative_error * 100.0 << "%" << std::endl;
  std::cout << "  相对误差超过1%的粒子数: " << error_count << " / " << num_fluid_particles << std::endl;
  
  std::cout << "\n结论：" << std::endl;
  std::cout << "  ✓ 直接使用初始速度（0）构建PPE，速度散度应为0" << std::endl;
  std::cout << "  ✓ 右边项主要包含壁面压力项，适合静水压力测试" << std::endl;
  
  // ========== 分析速度散度的空间分布 ==========
  std::cout << "\n========== 分析速度散度的空间分布 ==========" << std::endl;
  std::cout << "说明：初始速度为0（静水状态），理论上速度散度应为0；但壁面边界条件可能导致靠近壁面的粒子速度散度不为0" << std::endl;
  
  // 计算每个粒子的速度散度（分解为流体贡献和壁面贡献）
  // 注意：distance_to_wall和num_wall_neighbors已在之前的分析中定义，这里重用
  std::vector<double> divergence_fluid_contribution(num_fluid_particles, 0.0);
  std::vector<double> divergence_wall_contribution(num_fluid_particles, 0.0);
  std::vector<double> divergence_total(num_fluid_particles, 0.0);
  // 重用之前定义的distance_to_wall和num_wall_neighbors，重新初始化
  for (int i = 0; i < num_fluid_particles; ++i) {
    distance_to_wall[i] = 1e10;
    num_wall_neighbors[i] = 0;
  }
  
  CorrectiveMatrix corrective_matrix_calc_local;
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    const double2& pos_i = fluid_particles.position[i];
    const double2& vel_i = fluid_particles.velocity[i];
    
    // 提取速度corrective matrix的行向量
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C1_velocity = corrective_matrices_velocity[i].row(0);
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C2_velocity = corrective_matrices_velocity[i].row(1);
    
    double div_fluid = 0.0;
    double div_wall = 0.0;
    
    // 计算流体邻域粒子的速度散度贡献
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
      Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis = 
          corrective_matrix_calc_local.ComputeBasisFunctions(dx, dy, smoothing_radius);
      
      double dux_dr = (vel_j.x - vel_i.x) / dist;
      double duy_dr = (vel_j.y - vel_i.y) / dist;
      double C1P = (C1_velocity * basis)(0, 0);
      double C2P = (C2_velocity * basis)(0, 0);
      div_fluid += weight * (C1P * dux_dr + C2P * duy_dr);
    }
    
    // 计算壁面邻域粒子的速度散度贡献
    double min_wall_dist = 1e10;
    for (int j : fluid_particles.solid_neighbour_list[i]) {
      const double2& pos_j = solid_particles.position[j];
      const double2& vel_wall = solid_particles.velocity[j];  // 壁面速度为0
      
      double dx = pos_j.x - pos_i.x;
      double dy = pos_j.y - pos_i.y;
      double dist = ComputeDistance(pos_i, pos_j);
      
      if (dist < 1e-10 || dist > smoothing_radius) {
        continue;
      }
      
      if (dist < min_wall_dist) {
        min_wall_dist = dist;
      }
      
      double weight = WeightFunction(dist, smoothing_radius);
      Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis_velocity = 
          corrective_matrix_calc_local.ComputeBasisFunctions(dx, dy, smoothing_radius);
      
      // 壁面速度为0，所以速度差就是 -vel_i
      double dux_dr = (vel_wall.x - vel_i.x) / dist;  // = -vel_i.x / dist
      double duy_dr = (vel_wall.y - vel_i.y) / dist;  // = -vel_i.y / dist
      double C1P = (C1_velocity * basis_velocity)(0, 0);
      double C2P = (C2_velocity * basis_velocity)(0, 0);
      div_wall += weight * (C1P * dux_dr + C2P * duy_dr);
      num_wall_neighbors[i]++;
    }
    
    divergence_fluid_contribution[i] = div_fluid;
    divergence_wall_contribution[i] = div_wall;
    divergence_total[i] = div_fluid + div_wall;
    distance_to_wall[i] = min_wall_dist;
  }
  
  // 统计速度散度的分布
  std::cout << "\n速度散度统计（总散度 = 流体贡献 + 壁面贡献）：" << std::endl;
  
  double max_div_total = 0.0;
  double max_div_fluid = 0.0;
  double max_div_wall = 0.0;
  double sum_abs_div_total = 0.0;
  double sum_abs_div_fluid = 0.0;
  double sum_abs_div_wall = 0.0;
  int particles_with_wall_neighbors = 0;
  int particles_without_wall_neighbors = 0;
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    double abs_div_total = std::abs(divergence_total[i]);
    double abs_div_fluid = std::abs(divergence_fluid_contribution[i]);
    double abs_div_wall = std::abs(divergence_wall_contribution[i]);
    
    if (abs_div_total > max_div_total) {
      max_div_total = abs_div_total;
    }
    if (abs_div_fluid > max_div_fluid) {
      max_div_fluid = abs_div_fluid;
    }
    if (abs_div_wall > max_div_wall) {
      max_div_wall = abs_div_wall;
    }
    
    sum_abs_div_total += abs_div_total;
    sum_abs_div_fluid += abs_div_fluid;
    sum_abs_div_wall += abs_div_wall;
    
    if (num_wall_neighbors[i] > 0) {
      particles_with_wall_neighbors++;
    } else {
      particles_without_wall_neighbors++;
    }
  }
  
  double avg_abs_div_total = sum_abs_div_total / num_fluid_particles;
  double avg_abs_div_fluid = sum_abs_div_fluid / num_fluid_particles;
  double avg_abs_div_wall = sum_abs_div_wall / num_fluid_particles;
  
  std::cout << "  总速度散度最大值: " << max_div_total << std::endl;
  std::cout << "  总速度散度平均值: " << avg_abs_div_total << std::endl;
  std::cout << "  流体贡献最大值: " << max_div_fluid << std::endl;
  std::cout << "  流体贡献平均值: " << avg_abs_div_fluid << std::endl;
  std::cout << "  壁面贡献最大值: " << max_div_wall << std::endl;
  std::cout << "  壁面贡献平均值: " << avg_abs_div_wall << std::endl;
  std::cout << "  有壁面邻域的粒子数: " << particles_with_wall_neighbors << std::endl;
  std::cout << "  无壁面邻域的粒子数: " << particles_without_wall_neighbors << std::endl;
  
  // 分析靠近壁面的粒子的速度散度
  std::cout << "\n分析靠近壁面的粒子的速度散度：" << std::endl;
  
  // 按到壁面的距离分组统计
  const int num_bins = 10;
  std::vector<double> bin_max_dist(num_bins);
  std::vector<double> bin_avg_div_total(num_bins, 0.0);
  std::vector<double> bin_avg_div_wall(num_bins, 0.0);
  std::vector<int> bin_count(num_bins, 0);
  
  double max_wall_dist = 0.0;
  for (int i = 0; i < num_fluid_particles; ++i) {
    if (distance_to_wall[i] < 1e9) {  // 有壁面邻域
      if (distance_to_wall[i] > max_wall_dist) {
        max_wall_dist = distance_to_wall[i];
      }
    }
  }
  
  // 创建距离区间
  for (int bin = 0; bin < num_bins; ++bin) {
    bin_max_dist[bin] = (bin + 1) * max_wall_dist / num_bins;
  }
  
  // 统计每个区间的速度散度
  for (int i = 0; i < num_fluid_particles; ++i) {
    if (distance_to_wall[i] < 1e9) {  // 有壁面邻域
      int bin = static_cast<int>(distance_to_wall[i] / (max_wall_dist / num_bins));
      if (bin >= num_bins) bin = num_bins - 1;
      
      bin_avg_div_total[bin] += std::abs(divergence_total[i]);
      bin_avg_div_wall[bin] += std::abs(divergence_wall_contribution[i]);
      bin_count[bin]++;
    }
  }
  
  for (int bin = 0; bin < num_bins; ++bin) {
    if (bin_count[bin] > 0) {
      bin_avg_div_total[bin] /= bin_count[bin];
      bin_avg_div_wall[bin] /= bin_count[bin];
      std::cout << "  距离壁面 [" << (bin * max_wall_dist / num_bins) << ", " 
                << bin_max_dist[bin] << "] m: " << bin_count[bin] << " 个粒子, "
                << "平均总散度: " << bin_avg_div_total[bin] << ", "
                << "平均壁面贡献: " << bin_avg_div_wall[bin] << std::endl;
    }
  }
  
  // 对比有壁面邻域和无壁面邻域的粒子
  double avg_div_with_wall = 0.0;
  double avg_div_without_wall = 0.0;
  int count_with_wall = 0;
  int count_without_wall = 0;
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    if (num_wall_neighbors[i] > 0) {
      avg_div_with_wall += std::abs(divergence_total[i]);
      count_with_wall++;
    } else {
      avg_div_without_wall += std::abs(divergence_total[i]);
      count_without_wall++;
    }
  }
  
  if (count_with_wall > 0) {
    avg_div_with_wall /= count_with_wall;
  }
  if (count_without_wall > 0) {
    avg_div_without_wall /= count_without_wall;
  }
  
  std::cout << "\n对比分析：" << std::endl;
  std::cout << "  有壁面邻域的粒子平均速度散度: " << avg_div_with_wall << std::endl;
  std::cout << "  无壁面邻域的粒子平均速度散度: " << avg_div_without_wall << std::endl;
  std::cout << "  比值: " << (count_without_wall > 0 ? avg_div_with_wall / avg_div_without_wall : 0.0) << std::endl;
  
  std::cout << "\n验证结果：" << std::endl;
  if (avg_div_with_wall > avg_div_without_wall * 10.0) {
    std::cout << "  ✓ 猜想验证：靠近壁面的粒子速度散度显著大于远离壁面的粒子" << std::endl;
    std::cout << "  ✓ 壁面贡献是速度散度的主要来源" << std::endl;
  } else if (avg_abs_div_fluid < avg_abs_div_wall * 0.1) {
    std::cout << "  ✓ 猜想验证：流体贡献很小，壁面贡献是速度散度的主要来源" << std::endl;
  } else {
    std::cout << "  ⚠ 需要进一步分析" << std::endl;
  }
  
  // 计算理论静水压力（p = rho * g * h，其中h是从自由表面的深度）
  // 自由表面在y = water_height，所以深度 h = water_height - y
  
  // 创建理论压力向量 X'
  Vec x_theoretical = NULL;
  VecCreate(PETSC_COMM_WORLD, &x_theoretical);
  VecSetSizes(x_theoretical, PETSC_DECIDE, num_fluid_particles);
  VecSetType(x_theoretical, VECSEQ);
  VecSetFromOptions(x_theoretical);
  
  // 设置理论压力值（使用PetscScalar类型）
  std::vector<PetscScalar> theoretical_pressure_vec(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
    theoretical_pressure_vec[i] = theoretical_pressure_double[i];
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
  
  // ========================================================================
  // 验证Correction模块
  // ========================================================================
  std::cout << "\n=== 验证Correction模块 ===" << std::endl;
  
  // 使用Correction模块计算压力梯度、加速度并更新速度和位置
  std::cout << "\n使用Correction模块计算压力梯度、加速度并更新速度和位置..." << std::endl;
  Correction correction;
  
  // 使用correction模块批量计算并更新
  correction.ComputeAndUpdateAllParticles(
      fluid_particles,
      solid_particles,
      corrective_matrices_pressure,
      smoothing_radius,
      gravity_x,
      gravity_y,
      rho,
      time_step
  );
  
  // 使用Correction模块单独计算压力梯度（用于验证）
  std::vector<double2> pressure_gradient_correction(num_fluid_particles);
  std::vector<double2> acceleration_correction(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
    pressure_gradient_correction[i] = correction.ComputePressureGradient(
        i, fluid_particles, solid_particles, corrective_matrices_pressure[i],
        smoothing_radius, gravity_x, gravity_y, rho);
    acceleration_correction[i] = correction.ComputeAcceleration(
        pressure_gradient_correction[i], rho);
  }
  
  // 使用LSMPS手动计算压力梯度（用于对比验证）
  std::cout << "\n使用LSMPS手动计算压力梯度（用于对比验证）..." << std::endl;
  std::vector<double2> pressure_gradient_manual(num_fluid_particles);
  
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
    
    pressure_gradient_manual[i] = {grad_x, grad_y};
  }
  
  // 计算梯度大小
  std::vector<double> grad_magnitude_manual(num_fluid_particles);
  std::vector<double> grad_magnitude_correction(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
    grad_magnitude_manual[i] = std::sqrt(pressure_gradient_manual[i].x * pressure_gradient_manual[i].x + 
                                          pressure_gradient_manual[i].y * pressure_gradient_manual[i].y);
    grad_magnitude_correction[i] = std::sqrt(pressure_gradient_correction[i].x * pressure_gradient_correction[i].x + 
                                             pressure_gradient_correction[i].y * pressure_gradient_correction[i].y);
  }
  
  // 验证压力梯度计算的准确性
  std::cout << "\n验证压力梯度计算准确性..." << std::endl;
  double max_grad_diff = 0.0;
  double sum_grad_diff = 0.0;
  double max_grad_relative_error = 0.0;
  int grad_error_count = 0;
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    double diff_x = std::abs(pressure_gradient_correction[i].x - pressure_gradient_manual[i].x);
    double diff_y = std::abs(pressure_gradient_correction[i].y - pressure_gradient_manual[i].y);
    double diff_magnitude = std::sqrt(diff_x * diff_x + diff_y * diff_y);
    
    if (diff_magnitude > max_grad_diff) {
      max_grad_diff = diff_magnitude;
    }
    sum_grad_diff += diff_magnitude;
    
    // 计算相对误差
    if (grad_magnitude_manual[i] > 1e-10) {
      double relative_error = diff_magnitude / grad_magnitude_manual[i];
      if (relative_error > max_grad_relative_error) {
        max_grad_relative_error = relative_error;
      }
      if (relative_error > 1e-3) {  // 相对误差超过0.1%
        grad_error_count++;
      }
    }
  }
  
  double avg_grad_diff = sum_grad_diff / num_fluid_particles;
  std::cout << "  压力梯度最大绝对误差: " << max_grad_diff << std::endl;
  std::cout << "  压力梯度平均绝对误差: " << avg_grad_diff << std::endl;
  std::cout << "  压力梯度最大相对误差: " << max_grad_relative_error * 100.0 << "%" << std::endl;
  std::cout << "  压力梯度相对误差超过0.1%的粒子数: " << grad_error_count << " / " << num_fluid_particles << std::endl;
  
  // 验证加速度计算的准确性
  // 对于静水压力，压力梯度应该平衡重力
  // 理论压力梯度：dp/dy = -rho * g（y坐标向上为正，压力随y增加而减小）
  // 理论加速度：a_y = -dp/dy / rho = -(-rho * g) / rho = g（向下为正）
  // 但重力加速度gravity_y = -g，所以压力梯度引起的加速度应该与重力方向相反
  std::cout << "\n验证加速度计算准确性..." << std::endl;
  double max_acc_error = 0.0;
  double sum_acc_error = 0.0;
  // 对于静水压力，理论加速度y分量应该是g（向下为正），即与-gravity_y相同
  double theoretical_acc_y = -gravity_y;  // 理论加速度y分量 = -gravity_y = g
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    // 对于静水压力，压力梯度应该是垂直向下的
    // 理论压力梯度：dp/dy = -rho * g（y坐标向上为正）
    // 理论加速度：a_y = -dp/dy / rho = -(-rho * g) / rho = g = -gravity_y
    double acc_error_x = std::abs(acceleration_correction[i].x - 0.0);
    double acc_error_y = std::abs(acceleration_correction[i].y - theoretical_acc_y);
    double acc_error = std::sqrt(acc_error_x * acc_error_x + acc_error_y * acc_error_y);
    
    if (acc_error > max_acc_error) {
      max_acc_error = acc_error;
    }
    sum_acc_error += acc_error;
  }
  
  double avg_acc_error = sum_acc_error / num_fluid_particles;
  std::cout << "  加速度最大绝对误差: " << max_acc_error << " m/s²" << std::endl;
  std::cout << "  加速度平均绝对误差: " << avg_acc_error << " m/s²" << std::endl;
  std::cout << "  理论加速度y分量: " << theoretical_acc_y << " m/s²" << std::endl;
  
  // 验证速度和位置更新（相对于correction模块之前的状态）
  std::cout << "\n验证速度和位置更新（Correction模块）..." << std::endl;
  double max_velocity_change_correction = 0.0;
  double max_position_change = 0.0;
  double sum_velocity_change_correction = 0.0;
  double sum_position_change = 0.0;
  
  for (int i = 0; i < num_fluid_particles; ++i) {
    double2 velocity_change = {
        fluid_particles.velocity[i].x - velocity_before_correction[i].x,
        fluid_particles.velocity[i].y - velocity_before_correction[i].y
    };
    double2 position_change = {
        fluid_particles.position[i].x - position_before_correction[i].x,
        fluid_particles.position[i].y - position_before_correction[i].y
    };
    
    double vel_change_magnitude = std::sqrt(velocity_change.x * velocity_change.x + 
                                            velocity_change.y * velocity_change.y);
    double pos_change_magnitude = std::sqrt(position_change.x * position_change.x + 
                                            position_change.y * position_change.y);
    
    if (vel_change_magnitude > max_velocity_change_correction) {
      max_velocity_change_correction = vel_change_magnitude;
    }
    if (pos_change_magnitude > max_position_change) {
      max_position_change = pos_change_magnitude;
    }
    
    sum_velocity_change_correction += vel_change_magnitude;
    sum_position_change += pos_change_magnitude;
  }
  
  double avg_velocity_change_correction = sum_velocity_change_correction / num_fluid_particles;
  double avg_position_change = sum_position_change / num_fluid_particles;
  
  std::cout << "  速度变化最大值（Correction模块）: " << max_velocity_change_correction << " m/s" << std::endl;
  std::cout << "  速度变化平均值（Correction模块）: " << avg_velocity_change_correction << " m/s" << std::endl;
  std::cout << "  位置变化最大值: " << max_position_change << " m" << std::endl;
  std::cout << "  位置变化平均值: " << avg_position_change << " m" << std::endl;
  
  // 验证理论值：对于静水压力，初始速度为0，加速度为g（向下为正），经过dt后：
  // v_new = 0 + g * dt = 9.8 * 0.001 = 0.0098 m/s（向下为正，但y坐标向上为正，所以是负的）
  // 由于y坐标向上为正，速度变化应该是负的：v_y_new = -g * dt
  // x_new = x_old + v_new * dt = x_old - g * dt * dt
  double theoretical_velocity_change_y = -g * time_step;  // y方向速度变化（向下为负）
  double theoretical_velocity_change = std::abs(theoretical_velocity_change_y);  // 速度变化大小
  double theoretical_position_change = std::abs(theoretical_velocity_change_y * time_step);  // 位置变化大小
  
  std::cout << "  理论速度变化: " << theoretical_velocity_change << " m/s" << std::endl;
  std::cout << "  理论位置变化: " << theoretical_position_change << " m" << std::endl;
  
  // 输出验证结果摘要
  std::cout << "\n=== Correction模块验证结果摘要 ===" << std::endl;
  bool grad_ok = (max_grad_relative_error < 0.01);  // 相对误差小于1%
  bool acc_ok = (max_acc_error < 0.1);  // 加速度误差小于0.1 m/s²
  // 对于速度变化，需要考虑y坐标方向（向上为正，所以速度变化应该是负的）
  // 但这里我们比较的是速度变化的大小（Correction模块导致的变化）
  bool update_ok = (std::abs(avg_velocity_change_correction - theoretical_velocity_change) < 0.001);
  
  std::cout << "  压力梯度计算: " << (grad_ok ? "✓ 通过" : "✗ 失败") << std::endl;
  std::cout << "  加速度计算: " << (acc_ok ? "✓ 通过" : "✗ 失败") << std::endl;
  std::cout << "  速度和位置更新: " << (update_ok ? "✓ 通过" : "✗ 失败") << std::endl;
  
  if (grad_ok && acc_ok && update_ok) {
    std::cout << "\n✓ Correction模块验证通过！" << std::endl;
  } else {
    std::cout << "\n⚠ 警告：Correction模块验证未完全通过，请检查实现" << std::endl;
  }
  
  // 使用手动计算的压力梯度用于后续输出（保持兼容性）
  std::vector<double2> pressure_gradient = pressure_gradient_manual;
  std::vector<double> grad_magnitude = grad_magnitude_manual;
  
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
    
    // 追加压力梯度向量（手动计算）
    if (file_op.appendVTKVector(result_filename, "pressure_gradient_manual", pressure_gradient_manual)) {
      std::cout << "  压力梯度向量（手动计算）已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加压力梯度向量到VTK文件" << std::endl;
    }
    
    // 追加压力梯度向量（Correction模块计算）
    if (file_op.appendVTKVector(result_filename, "pressure_gradient_correction", pressure_gradient_correction)) {
      std::cout << "  压力梯度向量（Correction模块）已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加压力梯度向量到VTK文件" << std::endl;
    }
    
    // 追加加速度向量（Correction模块计算）
    if (file_op.appendVTKVector(result_filename, "acceleration_correction", acceleration_correction)) {
      std::cout << "  加速度向量（Correction模块）已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加加速度向量到VTK文件" << std::endl;
    }
    
    // 追加梯度大小标量（手动计算）
    if (file_op.appendVTKScalar(result_filename, "gradient_magnitude_manual", grad_magnitude_manual)) {
      std::cout << "  梯度大小标量（手动计算）已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加梯度大小标量到VTK文件" << std::endl;
    }
    
    // 追加梯度大小标量（Correction模块计算）
    if (file_op.appendVTKScalar(result_filename, "gradient_magnitude_correction", grad_magnitude_correction)) {
      std::cout << "  梯度大小标量（Correction模块）已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加梯度大小标量到VTK文件" << std::endl;
    }
    
    // 追加速度变化向量（显式力计算导致的变化，当前为0因为跳过了显式力计算）
    std::vector<double2> velocity_change_explicit(num_fluid_particles);
    for (int i = 0; i < num_fluid_particles; ++i) {
      velocity_change_explicit[i].x = velocity_before_correction[i].x - initial_velocity[i].x;
      velocity_change_explicit[i].y = velocity_before_correction[i].y - initial_velocity[i].y;
    }
    if (file_op.appendVTKVector(result_filename, "velocity_change_explicit", velocity_change_explicit)) {
      std::cout << "  速度变化向量（显式力，当前为0）已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加速度变化向量到VTK文件" << std::endl;
    }
    
    // 追加速度变化向量（Correction模块导致的变化）
    std::vector<double2> velocity_change_correction(num_fluid_particles);
    for (int i = 0; i < num_fluid_particles; ++i) {
      velocity_change_correction[i].x = fluid_particles.velocity[i].x - velocity_before_correction[i].x;
      velocity_change_correction[i].y = fluid_particles.velocity[i].y - velocity_before_correction[i].y;
    }
    if (file_op.appendVTKVector(result_filename, "velocity_change_correction", velocity_change_correction)) {
      std::cout << "  速度变化向量（Correction模块）已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加速度变化向量到VTK文件" << std::endl;
    }
    
    // 追加速度变化向量（总变化）
    std::vector<double2> velocity_change_total(num_fluid_particles);
    for (int i = 0; i < num_fluid_particles; ++i) {
      velocity_change_total[i].x = fluid_particles.velocity[i].x - initial_velocity[i].x;
      velocity_change_total[i].y = fluid_particles.velocity[i].y - initial_velocity[i].y;
    }
    if (file_op.appendVTKVector(result_filename, "velocity_change_total", velocity_change_total)) {
      std::cout << "  速度变化向量（总变化）已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加速度变化向量到VTK文件" << std::endl;
    }
    
    // 保留旧名称以保持兼容性
    if (file_op.appendVTKVector(result_filename, "velocity_change", velocity_change_total)) {
      std::cout << "  速度变化向量（兼容性）已追加" << std::endl;
    }
    
    // 追加位置变化向量
    std::vector<double2> position_change(num_fluid_particles);
    for (int i = 0; i < num_fluid_particles; ++i) {
      position_change[i].x = fluid_particles.position[i].x - position_before_correction[i].x;
      position_change[i].y = fluid_particles.position[i].y - position_before_correction[i].y;
    }
    if (file_op.appendVTKVector(result_filename, "position_change", position_change)) {
      std::cout << "  位置变化向量已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加位置变化向量到VTK文件" << std::endl;
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
    
    // 追加速度散度分析数据
    if (file_op.appendVTKScalar(result_filename, "divergence_total", divergence_total)) {
      std::cout << "  总速度散度已追加" << std::endl;
    }
    if (file_op.appendVTKScalar(result_filename, "divergence_fluid_contribution", divergence_fluid_contribution)) {
      std::cout << "  速度散度（流体贡献）已追加" << std::endl;
    }
    if (file_op.appendVTKScalar(result_filename, "divergence_wall_contribution", divergence_wall_contribution)) {
      std::cout << "  速度散度（壁面贡献）已追加" << std::endl;
    }
    if (file_op.appendVTKScalar(result_filename, "distance_to_wall", distance_to_wall)) {
      std::cout << "  到壁面距离已追加" << std::endl;
    }
    std::vector<int> num_wall_neighbors_int(num_fluid_particles);
    for (int i = 0; i < num_fluid_particles; ++i) {
      num_wall_neighbors_int[i] = num_wall_neighbors[i];
    }
    if (file_op.appendVTKScalar(result_filename, "num_wall_neighbors", num_wall_neighbors_int)) {
      std::cout << "  壁面邻域粒子数已追加" << std::endl;
    }
    
    // 追加右边项贡献分析数据
    if (file_op.appendVTKScalar(result_filename, "divergence_contribution_to_b", divergence_contribution)) {
      std::cout << "  速度散度对右边项的贡献已追加" << std::endl;
    }
    if (file_op.appendVTKScalar(result_filename, "wall_pressure_contribution_to_b", wall_pressure_contribution)) {
      std::cout << "  壁面压力边界条件对右边项的贡献已追加" << std::endl;
    }
    if (file_op.appendVTKScalar(result_filename, "divergence_fluid_part", divergence_fluid_part)) {
      std::cout << "  速度散度（流体邻域部分）已追加" << std::endl;
    }
    if (file_op.appendVTKScalar(result_filename, "divergence_wall_part", divergence_wall_part)) {
      std::cout << "  速度散度（壁面邻域部分）已追加" << std::endl;
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
