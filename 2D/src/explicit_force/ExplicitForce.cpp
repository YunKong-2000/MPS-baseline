#include "ExplicitForce.hpp"
#include <cmath>

namespace mps2D {

double2 ExplicitForce::ComputeVelocityLaplacian(
    int particle_idx,
    const std::vector<double2>& velocity_field,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const Eigen::Matrix<double, 5, 5>& corrective_matrix,
    double smoothing_radius) {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  const double2& v_i = velocity_field[particle_idx];
  
  // 提取corrective matrix的第3行和第4行（用于拉普拉斯算子计算）
  Eigen::Matrix<double, 1, 5> C3 = corrective_matrix.row(2);  // x²项
  Eigen::Matrix<double, 1, 5> C4 = corrective_matrix.row(3);  // y²项
  
  double laplacian_x = 0.0;
  double laplacian_y = 0.0;
  
  // 标量因子：2.0 / smoothing_radius
  const double scalar_factor = 2.0 / smoothing_radius;
  
  // 处理流体邻域粒子
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    
    // 跳过距离过小或过大的粒子
    if (dist < 1e-10 || dist > smoothing_radius) continue;
    
    const double2& v_j = velocity_field[j];
    double dv_x = v_j.x - v_i.x;
    double dv_y = v_j.y - v_i.y;
    double d_ij_x = dv_x / dist;  // (v_x_j - v_x_i) / r_ij
    double d_ij_y = dv_y / dist;  // (v_y_j - v_y_i) / r_ij
    double weight = WeightFunction(dist, smoothing_radius);
    
    // 计算基函数
    Eigen::Matrix<double, 5, 1> P;
    P << dx / dist,                                    // x/r
         dy / dist,                                    // y/r
         dx * dx / (dist * smoothing_radius),          // x²/(r*r_e)
         dy * dy / (dist * smoothing_radius),          // y²/(r*r_e)
         dx * dy / (dist * smoothing_radius);          // x*y/(r*r_e)
    
    // 计算拉普拉斯算子贡献
    double C3P = (C3 * P)(0, 0);
    double C4P = (C4 * P)(0, 0);
    
    // ∇²v_x = ∂²v_x/∂x² + ∂²v_x/∂y²
    double contrib_x_C3 = scalar_factor * weight * d_ij_x * C3P;  // ∂²v_x/∂x²
    double contrib_x_C4 = scalar_factor * weight * d_ij_x * C4P;  // ∂²v_x/∂y²
    laplacian_x += contrib_x_C3;
    laplacian_x += contrib_x_C4;
    
    // ∇²v_y = ∂²v_y/∂x² + ∂²v_y/∂y²
    double contrib_y_C3 = scalar_factor * weight * d_ij_y * C3P;  // ∂²v_y/∂x²
    double contrib_y_C4 = scalar_factor * weight * d_ij_y * C4P;  // ∂²v_y/∂y²
    laplacian_y += contrib_y_C3;
    laplacian_y += contrib_y_C4;
  }
  
  // 处理固体邻域粒子（壁面粒子）
  // 统一使用第一类边界条件（无滑移边界，壁面速度为零）
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    const double2& pos_j = solid_particles.position[j];
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    
    // 跳过距离过小或过大的粒子
    if (dist < 1e-10 || dist > smoothing_radius) continue;
    // 第一类边界条件：无滑移边界，壁面速度为零
    const double2 v_wall = {0.0, 0.0};
    
    double dv_x = v_wall.x - v_i.x;
    double dv_y = v_wall.y - v_i.y;
    double d_ij_x = dv_x / dist;
    double d_ij_y = dv_y / dist;
    double weight = WeightFunction(dist, smoothing_radius);
    
    // 使用标准基函数（与流体粒子相同）
    Eigen::Matrix<double, 5, 1> P;
    P << dx / dist,                                    // x/r
         dy / dist,                                    // y/r
         dx * dx / (dist * smoothing_radius),          // x²/(r*r_e)
         dy * dy / (dist * smoothing_radius),          // y²/(r*r_e)
         dx * dy / (dist * smoothing_radius);          // x*y/(r*r_e)
    
    // 计算拉普拉斯算子贡献
    double C3P = (C3 * P)(0, 0);
    double C4P = (C4 * P)(0, 0);
    
    double contrib_x_C3 = scalar_factor * weight * d_ij_x * C3P;
    double contrib_x_C4 = scalar_factor * weight * d_ij_x * C4P;
    double contrib_y_C3 = scalar_factor * weight * d_ij_y * C3P;
    double contrib_y_C4 = scalar_factor * weight * d_ij_y * C4P;
    
    laplacian_x += contrib_x_C3;
    laplacian_x += contrib_x_C4;
    laplacian_y += contrib_y_C3;
    laplacian_y += contrib_y_C4;
  }
  
  return {laplacian_x, laplacian_y};
}

double2 ExplicitForce::ComputeViscousAcceleration(
    const double2& velocity_laplacian,
    double kinematic_viscosity) {
  // 粘性力加速度 = ν * ∇²v
  // 其中 ν 是动力学粘性系数
  return {
      kinematic_viscosity * velocity_laplacian.x,
      kinematic_viscosity * velocity_laplacian.y
  };
}

double2 ExplicitForce::ComputeGravityAcceleration(
    double gravity_x,
    double gravity_y) {
  // 重力加速度直接使用输入值
  return {gravity_x, gravity_y};
}

void ExplicitForce::UpdateVelocity(
    int particle_idx,
    FluidParticle& fluid_particles,
    const double2& viscous_acceleration,
    const double2& gravity_acceleration,
    double time_step) {
  
  // 总加速度 = 粘性力加速度 + 重力加速度
  double2 total_acceleration = {
      viscous_acceleration.x + gravity_acceleration.x,
      viscous_acceleration.y + gravity_acceleration.y
  };
  
  // 显式时间积分：只更新速度作为临时速度，不更新位置
  // v_temp = v_old + a * dt
  double2& velocity = fluid_particles.velocity[particle_idx];
  velocity.x += total_acceleration.x * time_step;
  velocity.y += total_acceleration.y * time_step;
}

void ExplicitForce::ComputeAndUpdateVelocity(
    FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const std::vector<Eigen::Matrix<double, 5, 5>>& corrective_matrices,
    double smoothing_radius,
    double kinematic_viscosity,
    double gravity_x,
    double gravity_y,
    double time_step) {
  
  int num_particles = fluid_particles.particle_num;
  
  // 计算重力加速度（对所有粒子相同）
  double2 gravity_acceleration = ComputeGravityAcceleration(gravity_x, gravity_y);
  
  // 第一步：先计算所有粒子的速度拉普拉斯算子和粘性力加速度（使用原始速度）
  std::vector<double2> viscous_accelerations(num_particles);
  for (int i = 0; i < num_particles; ++i) {
    // 计算速度拉普拉斯算子（此时所有粒子的速度都还是原始值）
    double2 velocity_laplacian = ComputeVelocityLaplacian(
        i,
        fluid_particles.velocity,
        fluid_particles,
        solid_particles,
        corrective_matrices[i],
        smoothing_radius);
    
    // 计算粘性力加速度
    viscous_accelerations[i] = ComputeViscousAcceleration(
        velocity_laplacian,
        kinematic_viscosity);
  }
  
  // 第二步：统一更新所有粒子的速度（不更新位置）
  for (int i = 0; i < num_particles; ++i) {
    UpdateVelocity(
        i,
        fluid_particles,
        viscous_accelerations[i],
        gravity_acceleration,
        time_step);
  }
}

void ExplicitForce::ComputeAndUpdateVelocity(
    FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const std::vector<Eigen::Matrix<double, 5, 5>>& corrective_matrices,
    double smoothing_radius,
    double kinematic_viscosity,
    double gravity_x,
    double gravity_y,
    double time_step,
    std::vector<double2>& viscous_acceleration) {
  
  int num_particles = fluid_particles.particle_num;
  
  // 确保输出向量大小正确
  viscous_acceleration.resize(num_particles);
  
  // 计算重力加速度（对所有粒子相同）
  double2 gravity_acceleration = ComputeGravityAcceleration(gravity_x, gravity_y);
  
  // 第一步：先计算所有粒子的速度拉普拉斯算子和粘性力加速度（使用原始速度）
  for (int i = 0; i < num_particles; ++i) {
    // 计算速度拉普拉斯算子（此时所有粒子的速度都还是原始值）
    double2 velocity_laplacian = ComputeVelocityLaplacian(
        i,
        fluid_particles.velocity,
        fluid_particles,
        solid_particles,
        corrective_matrices[i],
        smoothing_radius);
    
    // 计算粘性力加速度
    viscous_acceleration[i] = ComputeViscousAcceleration(
        velocity_laplacian,
        kinematic_viscosity);
  }
  
  // 第二步：统一更新所有粒子的速度（不更新位置）
  for (int i = 0; i < num_particles; ++i) {
    UpdateVelocity(
        i,
        fluid_particles,
        viscous_acceleration[i],
        gravity_acceleration,
        time_step);
  }
}

} // namespace mps2D
