#include "../src/explicit_force/ExplicitForce.hpp"
#include "../src/lsmps/CorrectiveMatrix.hpp"
#include "../src/core/Particle.hpp"
#include "core/Types.h"
#include "core/MPSUtils.h"
#include "../src/neighbour_list/NeighborListSearcher.hpp"
#include "../src/surface_detection/SurfaceDetector.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <limits>

using namespace mps2D;

// 生成管道流动测试场景
// 参数：
//   channel_length: 管道长度
//   channel_height: 管道高度
//   particle_spacing: 粒子间距
//   fluid_particles: 输出流体粒子
//   solid_particles: 输出固体粒子（壁面）
void GeneratePipeFlowTest(
    double channel_length,
    double channel_height,
    double particle_spacing,
    FluidParticle& fluid_particles,
    SolidParticle& solid_particles) {
  
  // 计算粒子数量
  int nx_fluid = static_cast<int>(channel_length / particle_spacing) + 1;
  int ny_fluid = static_cast<int>(channel_height / particle_spacing) - 1;  // 排除上下壁面
  int num_fluid = nx_fluid * ny_fluid;
  
  // 生成流体粒子
  fluid_particles.particle_num = num_fluid;
  fluid_particles.position.resize(num_fluid);
  fluid_particles.velocity.resize(num_fluid);
  fluid_particles.pressure.resize(num_fluid);
  fluid_particles.density.resize(num_fluid);
  fluid_particles.surface_type.resize(num_fluid);
  fluid_particles.fluid_neighbour_list.resize(num_fluid);
  fluid_particles.solid_neighbour_list.resize(num_fluid);
  
  int idx = 0;
  for (int j = 1; j <= ny_fluid; ++j) {
    for (int i = 0; i < nx_fluid; ++i) {
      double x = i * particle_spacing;
      double y = j * particle_spacing;
      fluid_particles.position[idx] = {x, y};
      fluid_particles.velocity[idx] = {0.0, 0.0};
      fluid_particles.pressure[idx] = 0.0;
      fluid_particles.density[idx] = 1000.0;
      fluid_particles.surface_type[idx] = SurfaceType::INNER;
      ++idx;
    }
  }
  
  // 生成固体粒子（上下壁面）
  int nx_bottom = static_cast<int>(channel_length / particle_spacing) + 1;
  int nx_top = nx_bottom;
  int num_solid = nx_bottom + nx_top;
  solid_particles.particle_num = num_solid;
  solid_particles.position.resize(num_solid);
  solid_particles.velocity.resize(num_solid);
  solid_particles.normal_vector.resize(num_solid);
  
  idx = 0;
  // 下壁面（法向量向上）
  for (int i = 0; i < nx_bottom; ++i) {
    double x = i * particle_spacing;
    double y = 0;
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {0.0, 1.0};
    ++idx;
  }
  
  // 上壁面（法向量向下）
  for (int i = 0; i < nx_top; ++i) {
    double x = i * particle_spacing;
    double y = channel_height;
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {0.0, -1.0};
    ++idx;
  }
}

// 计算理论速度分布（Poiseuille流动）
// v_x(y) = v_max * (1 - (2y/h - 1)^2)
// v_y = 0
double2 ComputeTheoreticalVelocity(double y, double channel_height, double v_max) {
  double2 velocity;
  if (y < 0.0) y = 0.0;
  if (y > channel_height) y = channel_height;
  double normalized_y = 2.0 * y / channel_height - 1.0;
  velocity.x = v_max * (1.0 - normalized_y * normalized_y);
  velocity.y = 0.0;
  return velocity;
}

// 计算理论速度拉普拉斯算子
// 对于Poiseuille流动：v_x(y) = v_max * (1 - (2y/h - 1)^2)
// ∂^2v_x/∂x^2 = 0
// ∂^2v_x/∂y^2 = -8*v_max/h^2
double2 ComputeTheoreticalVelocityLaplacian(
    double /* y */, double channel_height, double v_max) {
  double laplacian_x = -8.0 * v_max / (channel_height * channel_height);
  double laplacian_y = 0.0;
  return {laplacian_x, laplacian_y};
}

// 输出结果到VTK文件
void WriteResultsToVTK(
    const std::string& filename,
    const FluidParticle& fluid_particles,
    const std::vector<double2>& initial_velocity,
    const std::vector<double2>& final_velocity,
    const std::vector<double2>& initial_position,
    const std::vector<double2>& final_position,
    const std::vector<double2>& velocity_laplacian,
    const std::vector<double2>& viscous_acceleration,
    const std::vector<double2>& gravity_acceleration) {
  
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "错误：无法打开文件 " << filename << std::endl;
    return;
  }
  
  int num_particles = fluid_particles.particle_num;
  file << std::fixed << std::setprecision(15);
  
  // VTK文件头
  file << "# vtk DataFile Version 3.0\n";
  file << "Explicit Force Test - MPS 2D\n";
  file << "ASCII\n";
  file << "DATASET POLYDATA\n";
  
  // 写入点坐标（使用最终位置）
  file << "POINTS " << num_particles << " float\n";
  for (int i = 0; i < num_particles; ++i) {
    const auto& pos = final_position[i];
    file << pos.x << " " << pos.y << " 0.0\n";
  }
  
  // 写入顶点
  file << "VERTICES " << num_particles << " " << (num_particles * 2) << "\n";
  for (int i = 0; i < num_particles; ++i) {
    file << "1 " << i << "\n";
  }
  
  // 写入点数据
  file << "POINT_DATA " << num_particles << "\n";
  
  // 写入初始速度
  file << "VECTORS initial_velocity float\n";
  for (int i = 0; i < num_particles; ++i) {
    file << initial_velocity[i].x << " " 
         << initial_velocity[i].y << " 0.0\n";
  }
  
  // 写入最终速度
  file << "VECTORS final_velocity float\n";
  for (int i = 0; i < num_particles; ++i) {
    file << final_velocity[i].x << " " 
         << final_velocity[i].y << " 0.0\n";
  }
  
  // 写入速度变化
  file << "VECTORS velocity_change float\n";
  for (int i = 0; i < num_particles; ++i) {
    double dv_x = final_velocity[i].x - initial_velocity[i].x;
    double dv_y = final_velocity[i].y - initial_velocity[i].y;
    file << dv_x << " " << dv_y << " 0.0\n";
  }
  
  // 写入初始位置
  file << "VECTORS initial_position float\n";
  for (int i = 0; i < num_particles; ++i) {
    file << initial_position[i].x << " " 
         << initial_position[i].y << " 0.0\n";
  }
  
  // 写入位置变化
  file << "VECTORS position_change float\n";
  for (int i = 0; i < num_particles; ++i) {
    double dx = final_position[i].x - initial_position[i].x;
    double dy = final_position[i].y - initial_position[i].y;
    file << dx << " " << dy << " 0.0\n";
  }
  
  // 写入速度拉普拉斯算子
  file << "VECTORS velocity_laplacian float\n";
  for (int i = 0; i < num_particles; ++i) {
    file << velocity_laplacian[i].x << " " 
         << velocity_laplacian[i].y << " 0.0\n";
  }
  
  // 写入粘性力加速度
  file << "VECTORS viscous_acceleration float\n";
  for (int i = 0; i < num_particles; ++i) {
    file << viscous_acceleration[i].x << " " 
         << viscous_acceleration[i].y << " 0.0\n";
  }
  
  // 写入重力加速度
  file << "VECTORS gravity_acceleration float\n";
  for (int i = 0; i < num_particles; ++i) {
    file << gravity_acceleration[i].x << " " 
         << gravity_acceleration[i].y << " 0.0\n";
  }
  
  // 写入总加速度
  file << "VECTORS total_acceleration float\n";
  for (int i = 0; i < num_particles; ++i) {
    double total_acc_x = viscous_acceleration[i].x + gravity_acceleration[i].x;
    double total_acc_y = viscous_acceleration[i].y + gravity_acceleration[i].y;
    file << total_acc_x << " " << total_acc_y << " 0.0\n";
  }
  
  // 写入Y坐标
  file << "SCALARS y_coordinate float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    file << initial_position[i].y << "\n";
  }
  
  file.close();
  std::cout << "  已输出VTK文件: " << filename << std::endl;
}

int main() {
  std::cout << "=== 显式力计算模块测试 ===" << std::endl;
  
  // 几何参数
  const double channel_length = 2.0;   // 管道长度 (m)
  const double channel_height = 0.5;   // 管道高度 (m)
  
  // 速度参数
  const double v_max = 1.0;  // 最大速度（管道中心）(m/s)
  
  // 物理参数
  const double kinematic_viscosity = 0.001;  // 动力学粘性系数 (m²/s)
  const double gravity_x = 0.0;               // 重力加速度x分量 (m/s²)
  const double gravity_y = -9.8;              // 重力加速度y分量 (m/s²)
  const double time_step = 0.001;             // 时间步长 (s)
  const int num_steps = 1;                    // 时间步数
  
  // 粒子参数
  const double particle_spacing = 0.02;  // 粒子间距 (m)
  const double particle_radius = particle_spacing / 2.0;
  const double smoothing_radius = 2.1 * particle_spacing;
  const double cell_size = 2.0 * smoothing_radius;
  
  std::cout << "\n几何参数:" << std::endl;
  std::cout << "  管道长度: " << channel_length << " m" << std::endl;
  std::cout << "  管道高度: " << channel_height << " m" << std::endl;
  
  std::cout << "\n物理参数:" << std::endl;
  std::cout << "  动力学粘性系数: " << kinematic_viscosity << " m²/s" << std::endl;
  std::cout << "  重力加速度: (" << gravity_x << ", " << gravity_y << ") m/s²" << std::endl;
  std::cout << "  时间步长: " << time_step << " s" << std::endl;
  std::cout << "  时间步数: " << num_steps << std::endl;
  
  std::cout << "\n粒子参数:" << std::endl;
  std::cout << "  粒子间距: " << particle_spacing << " m" << std::endl;
  std::cout << "  平滑半径: " << smoothing_radius << " m" << std::endl;
  
  // 生成测试场景
  FluidParticle fluid_particles("fluid");
  SolidParticle solid_particles("solid");
  
  std::cout << "\n生成测试场景..." << std::endl;
  GeneratePipeFlowTest(
      channel_length, channel_height, particle_spacing,
      fluid_particles, solid_particles);
  
  std::cout << "  流体粒子数: " << fluid_particles.particle_num << std::endl;
  std::cout << "  固体粒子数: " << solid_particles.particle_num << std::endl;
  
  // 构建邻居列表
  std::cout << "\n构建邻居列表..." << std::endl;
  NeighborListSearcher neighbor_searcher;
  neighbor_searcher.BuildNeighborList(
      fluid_particles, solid_particles,
      particle_radius, smoothing_radius, cell_size);
  
  // 判定自由面
  std::cout << "判定自由面..." << std::endl;
  SurfaceDetector surface_detector;
  surface_detector.DetectSurfaceParticles(
      fluid_particles, solid_particles,
      smoothing_radius, particle_spacing);
  
  // 设置初始速度场（Poiseuille流动）
  std::cout << "\n设置初始速度场（Poiseuille流动）..." << std::endl;
  std::vector<double2> initial_velocity(fluid_particles.particle_num);
  std::vector<double2> initial_position(fluid_particles.particle_num);
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    double y = fluid_particles.position[i].y;
    initial_velocity[i] = ComputeTheoreticalVelocity(y, channel_height, v_max);
    initial_position[i] = fluid_particles.position[i];
    fluid_particles.velocity[i] = initial_velocity[i];
  }
  
  // 计算corrective matrix
  std::cout << "\n计算corrective matrix..." << std::endl;
  CorrectiveMatrix corrective_matrix_calculator;
  std::vector<Eigen::Matrix<double, 5, 5>> matrices(fluid_particles.particle_num);
  
  int num_valid_matrix = 0;
  // 使用第一类边界条件（无滑移边界）
  const bool border_condition = false;
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    matrices[i] = corrective_matrix_calculator.ComputeCorrectiveMatrix(
        i, fluid_particles, solid_particles, smoothing_radius, border_condition);
    
    // 检查是否是单位矩阵
    bool is_identity = true;
    const double tolerance = 1e-6;
    for (int row = 0; row < 5; ++row) {
      for (int col = 0; col < 5; ++col) {
        double expected = (row == col) ? 1.0 : 0.0;
        if (std::abs(matrices[i](row, col) - expected) > tolerance) {
          is_identity = false;
          break;
        }
      }
      if (!is_identity) break;
    }
    if (!is_identity) ++num_valid_matrix;
  }
  std::cout << "  有效corrective matrix: " << num_valid_matrix 
            << " / " << fluid_particles.particle_num << std::endl;
  
  // 创建显式力计算器
  ExplicitForce explicit_force;
  
  // 存储中间结果用于分析
  std::vector<double2> velocity_laplacian(fluid_particles.particle_num);
  std::vector<double2> viscous_acceleration(fluid_particles.particle_num);
  double2 gravity_acceleration = explicit_force.ComputeGravityAcceleration(
      gravity_x, gravity_y);
  
  // 计算初始状态的速度拉普拉斯算子和粘性力加速度
  std::cout << "\n计算初始状态的速度拉普拉斯算子和粘性力加速度..." << std::endl;
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    velocity_laplacian[i] = explicit_force.ComputeVelocityLaplacian(
        i,
        fluid_particles.velocity,
        fluid_particles,
        solid_particles,
        matrices[i],
        smoothing_radius);
    
    viscous_acceleration[i] = explicit_force.ComputeViscousAcceleration(
        velocity_laplacian[i],
        kinematic_viscosity);
  }
  
  // 验证速度拉普拉斯算子（与理论值比较）
  std::cout << "\n验证速度拉普拉斯算子（与理论值比较）..." << std::endl;
  double total_error_lap_x = 0.0, total_error_lap_y = 0.0;
  double max_error_lap_x = 0.0, max_error_lap_y = 0.0;
  int num_valid = 0;
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    if (fluid_particles.surface_type[i] == SurfaceType::INNER) {
      double y = initial_position[i].y;
      double2 theoretical_laplacian = ComputeTheoreticalVelocityLaplacian(
          y, channel_height, v_max);
      
      double error_x = std::abs(velocity_laplacian[i].x - theoretical_laplacian.x);
      double error_y = std::abs(velocity_laplacian[i].y - theoretical_laplacian.y);
      
      total_error_lap_x += error_x;
      total_error_lap_y += error_y;
      
      if (error_x > max_error_lap_x) max_error_lap_x = error_x;
      if (error_y > max_error_lap_y) max_error_lap_y = error_y;
      
      ++num_valid;
    }
  }
  
  double avg_error_lap_x = total_error_lap_x / num_valid;
  double avg_error_lap_y = total_error_lap_y / num_valid;
  
  std::cout << "  统计粒子数（内部粒子）: " << num_valid << std::endl;
  std::cout << "  ∇²v_x 平均误差: " << avg_error_lap_x << " 1/(m²·s)" << std::endl;
  std::cout << "  ∇²v_y 平均误差: " << avg_error_lap_y << " 1/(m²·s)" << std::endl;
  std::cout << "  ∇²v_x 最大误差: " << max_error_lap_x << " 1/(m²·s)" << std::endl;
  std::cout << "  ∇²v_y 最大误差: " << max_error_lap_y << " 1/(m²·s)" << std::endl;
  
  // 执行显式时间积分
  std::cout << "\n执行显式时间积分..." << std::endl;
  explicit_force.ComputeAndUpdateAllParticles(
      fluid_particles,
      solid_particles,
      matrices,
      smoothing_radius,
      kinematic_viscosity,
      gravity_x,
      gravity_y,
      time_step);
  
  // 保存最终状态
  std::vector<double2> final_velocity(fluid_particles.particle_num);
  std::vector<double2> final_position(fluid_particles.particle_num);
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    final_velocity[i] = fluid_particles.velocity[i];
    final_position[i] = fluid_particles.position[i];
  }
  
  // 统计速度和位置变化
  std::cout << "\n统计速度和位置变化..." << std::endl;
  double max_velocity_change = 0.0;
  double max_position_change = 0.0;
  double total_velocity_change = 0.0;
  double total_position_change = 0.0;
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    double dv_x = final_velocity[i].x - initial_velocity[i].x;
    double dv_y = final_velocity[i].y - initial_velocity[i].y;
    double velocity_change_mag = std::sqrt(dv_x * dv_x + dv_y * dv_y);
    
    double dx = final_position[i].x - initial_position[i].x;
    double dy = final_position[i].y - initial_position[i].y;
    double position_change_mag = std::sqrt(dx * dx + dy * dy);
    
    total_velocity_change += velocity_change_mag;
    total_position_change += position_change_mag;
    
    if (velocity_change_mag > max_velocity_change) {
      max_velocity_change = velocity_change_mag;
    }
    if (position_change_mag > max_position_change) {
      max_position_change = position_change_mag;
    }
  }
  
  double avg_velocity_change = total_velocity_change / fluid_particles.particle_num;
  double avg_position_change = total_position_change / fluid_particles.particle_num;
  
  std::cout << "  平均速度变化: " << avg_velocity_change << " m/s" << std::endl;
  std::cout << "  最大速度变化: " << max_velocity_change << " m/s" << std::endl;
  std::cout << "  平均位置变化: " << avg_position_change << " m" << std::endl;
  std::cout << "  最大位置变化: " << max_position_change << " m" << std::endl;
  
  // 输出结果到VTK文件
  std::cout << "\n输出结果到VTK文件..." << std::endl;
  WriteResultsToVTK(
      "explicit_force_test.vtk",
      fluid_particles,
      initial_velocity,
      final_velocity,
      initial_position,
      final_position,
      velocity_laplacian,
      viscous_acceleration,
      std::vector<double2>(fluid_particles.particle_num, gravity_acceleration));
  
  std::cout << "\n测试完成!" << std::endl;
  
  return 0;
}
