#include "Correction.hpp"
#include "../lsmps/CorrectiveMatrix.hpp"
#include <cmath>
#include <iostream>

namespace mps2D {

double2 Correction::ComputePressureGradient(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>& corrective_matrix,
    double smoothing_radius,
    double gravity_x,
    double gravity_y,
    double density) {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  double p_i = fluid_particles.pressure[particle_idx];
  
  // 提取压力corrective matrix的前两行（用于梯度计算，第二类边界条件）
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C1 = corrective_matrix.row(0);  // x方向
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C2 = corrective_matrix.row(1);  // y方向
  
  double grad_x = 0.0;
  double grad_y = 0.0;
  
  CorrectiveMatrix corrective_matrix_calc;
  
  // 处理流体邻域粒子
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];
    double p_j = fluid_particles.pressure[j];
    
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
  // 壁面处压力梯度的法向分量：dp/dn = ρ * n · g
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
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
    
    // 壁面处压力梯度的法向分量：dp/dn = ρ * n · g
    // 归一化法向量
    double nn = std::sqrt(normal.x * normal.x + normal.y * normal.y);
    if (nn < 1e-10) {
      continue;  // 跳过无效的法向量
    }
    double n_x = normal.x / nn;
    double n_y = normal.y / nn;
    
    // 计算法向压力梯度：dp/dn = ρ * (n · g) = ρ * (n_x * g_x + n_y * g_y)
    double dp_dn = density * (n_x * gravity_x + n_y * gravity_y);
    
    // 计算梯度贡献
    grad_x += weight * dp_dn * (C1 * basis_wall)(0, 0);
    grad_y += weight * dp_dn * (C2 * basis_wall)(0, 0);
  }
  
  return {grad_x, grad_y};
}

double2 Correction::ComputeAcceleration(
    const double2& pressure_gradient,
    double density) {
  
  // 加速度 = -(压力梯度/密度)
  if (density < 1e-10) {
    std::cerr << "警告：密度过小，可能导致数值不稳定" << std::endl;
    return {0.0, 0.0};
  }
  
  return {-pressure_gradient.x / density, -pressure_gradient.y / density};
}

void Correction::UpdateVelocityAndPosition(
    int particle_idx,
    FluidParticle& fluid_particles,
    const double2& acceleration,
    double time_step) {
  
  // 更新速度：v_new = v_old + a * dt
  double2& velocity = fluid_particles.velocity[particle_idx];
  velocity.x += acceleration.x * time_step;
  velocity.y += acceleration.y * time_step;
  
  // 更新位置：x_new = x_old + v_new * dt
  double2& position = fluid_particles.position[particle_idx];
  position.x += velocity.x * time_step;
  position.y += velocity.y * time_step;
}

void Correction::ComputeAndUpdateAllParticles(
    FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>>& corrective_matrices,
    double smoothing_radius,
    double gravity_x,
    double gravity_y,
    double density,
    double time_step) {
  
  int num_particles = fluid_particles.particle_num;
  
  if (static_cast<int>(corrective_matrices.size()) != num_particles) {
    std::cerr << "错误：corrective_matrices数量与流体粒子数不匹配" << std::endl;
    return;
  }
  
  // 对每个粒子计算压力梯度、加速度并更新速度和位置
  for (int i = 0; i < num_particles; ++i) {
    // 计算压力梯度
    double2 pressure_gradient = ComputePressureGradient(
        i,
        fluid_particles,
        solid_particles,
        corrective_matrices[i],
        smoothing_radius,
        gravity_x,
        gravity_y,
        density);
    
    // 计算加速度
    double2 acceleration = ComputeAcceleration(pressure_gradient, density);
    
    // 更新速度和位置
    UpdateVelocityAndPosition(i, fluid_particles, acceleration, time_step);
  }
}

} // namespace mps2D

