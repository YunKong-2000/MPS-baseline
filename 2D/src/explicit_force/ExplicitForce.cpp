#include "ExplicitForce.hpp"
#include <cmath>

namespace mps2D {

namespace {
constexpr double kSplashRepulsionCoeff = 20.0;
constexpr double kSplashSafetyDistanceRatio = 0.9;
constexpr double kSmallEps = 1e-10;
}  // namespace

double2 ExplicitForce::ComputeVelocityLaplacian(
    int particle_idx,
    const std::vector<double2>& velocity_field,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const Eigen::Matrix<double, 5, 5>& moment_matrix_inverse,
    double smoothing_radius) {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  const double2& v_i = velocity_field[particle_idx];
  
  // 提取 moment matrix 逆矩阵的第3行和第4行（用于拉普拉斯算子计算）
  // 文档中对应 [M_{i,2} + M_{i,3}]
  Eigen::RowVector<double, 5> M2 = moment_matrix_inverse.row(2);  // x² 项
  Eigen::RowVector<double, 5> M3 = moment_matrix_inverse.row(3);  // y² 项
  Eigen::RowVector<double, 5> M2_plus_M3 = M2 + M3;
  
  double laplacian_x = 0.0;
  double laplacian_y = 0.0;
  
  // 标量因子：2 / r_s^2
  // 其中 r_s 为平滑距离（与文档中的 r_s 一致）
  const double scalar_factor = 2.0 / (smoothing_radius * smoothing_radius);

  // 使用与 LSMPS moment 矩阵相同的基函数定义
  // P_ij = [x/r_s, y/r_s, x^2/r_s^2, y^2/r_s^2, x y / r_s^2]
  CorrectiveMatrix corrective_matrix_calc;
  
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
    double weight = WeightFunction(dist, smoothing_radius);

    // 计算基函数 P_ij
    Eigen::Vector<double, 5> basis =
        corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);

    // 根据文档：∇²φ_i = (2/r_s^2) Σ_j w_ij (φ_j - φ_i) [M_{i,2} + M_{i,3}] P_ij
    double M2M3P = (M2_plus_M3 * basis)(0, 0);
    double coeff = scalar_factor * weight * M2M3P;

    laplacian_x += coeff * dv_x;
    laplacian_y += coeff * dv_y;
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

    double weight = WeightFunction(dist, smoothing_radius);

    // 使用与流体粒子相同的基函数 P_ij
    Eigen::Vector<double, 5> basis =
        corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);

    double M2M3P = (M2_plus_M3 * basis)(0, 0);
    double coeff = scalar_factor * weight * M2M3P;

    laplacian_x += coeff * dv_x;
    laplacian_y += coeff * dv_y;
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

double2 ExplicitForce::ComputeSplashRepulsiveAcceleration(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double particle_spacing) const {
  if (particle_idx < 0 ||
      particle_idx >= static_cast<int>(fluid_particles.particle_num) ||
      particle_idx >= static_cast<int>(fluid_particles.surface_type.size()) ||
      particle_idx >= static_cast<int>(fluid_particles.solid_neighbour_list.size()) ||
      fluid_particles.surface_type[particle_idx] != SurfaceType::SPLASH ||
      particle_spacing <= kSmallEps) {
    return {0.0, 0.0};
  }

  const double safety_distance = kSplashSafetyDistanceRatio * particle_spacing;
  if (safety_distance <= kSmallEps) {
    return {0.0, 0.0};
  }

  const double2& pos_i = fluid_particles.position[particle_idx];
  double2 repulsive_acceleration = {0.0, 0.0};

  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    if (j < 0 ||
        j >= static_cast<int>(solid_particles.position.size()) ||
        j >= static_cast<int>(solid_particles.normal_vector.size())) {
      continue;
    }

    const double2& pos_j = solid_particles.position[j];
    const double dist = ComputeDistance(pos_i, pos_j);
    if (dist <= kSmallEps || dist >= safety_distance) {
      continue;
    }

    double2 normal = solid_particles.normal_vector[j];
    const double normal_norm =
        std::sqrt(normal.x * normal.x + normal.y * normal.y);
    if (normal_norm > kSmallEps) {
      normal.x /= normal_norm;
      normal.y /= normal_norm;
    } else {
      const double dx = pos_i.x - pos_j.x;
      const double dy = pos_i.y - pos_j.y;
      const double fallback_norm = std::sqrt(dx * dx + dy * dy);
      if (fallback_norm <= kSmallEps) {
        continue;
      }
      normal = {dx / fallback_norm, dy / fallback_norm};
    }

    const double strength =
        kSplashRepulsionCoeff * (safety_distance - dist) / safety_distance;
    repulsive_acceleration.x += strength * normal.x;
    repulsive_acceleration.y += strength * normal.y;
  }

  return repulsive_acceleration;
}

void ExplicitForce::ComputeAndUpdateVelocity(
    FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const std::vector<Eigen::Matrix<double, 5, 5>>& corrective_matrices,
    double smoothing_radius,
    double particle_spacing,
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
    if (i < static_cast<int>(fluid_particles.surface_type.size()) &&
        fluid_particles.surface_type[i] == SurfaceType::SPLASH) {
      viscous_accelerations[i] = {0.0, 0.0};
      continue;
    }
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
    double2 non_pressure_acceleration = viscous_accelerations[i];
    if (i < static_cast<int>(fluid_particles.surface_type.size()) &&
        fluid_particles.surface_type[i] == SurfaceType::SPLASH) {
      non_pressure_acceleration = ComputeSplashRepulsiveAcceleration(
          i, fluid_particles, solid_particles, particle_spacing);
    }
    UpdateVelocity(
        i,
        fluid_particles,
        non_pressure_acceleration,
        gravity_acceleration,
        time_step);
  }
}

void ExplicitForce::ComputeAndUpdateVelocity(
    FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const std::vector<Eigen::Matrix<double, 5, 5>>& corrective_matrices,
    double smoothing_radius,
    double particle_spacing,
    double kinematic_viscosity,
    double gravity_x,
    double gravity_y,
    double time_step,
    std::vector<double2>& viscous_acceleration) {
  std::vector<double2> splash_repulsive_acceleration;
  ComputeAndUpdateVelocity(
      fluid_particles,
      solid_particles,
      corrective_matrices,
      smoothing_radius,
      particle_spacing,
      kinematic_viscosity,
      gravity_x,
      gravity_y,
      time_step,
      viscous_acceleration,
      splash_repulsive_acceleration);
}

void ExplicitForce::ComputeAndUpdateVelocity(
    FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const std::vector<Eigen::Matrix<double, 5, 5>>& corrective_matrices,
    double smoothing_radius,
    double particle_spacing,
    double kinematic_viscosity,
    double gravity_x,
    double gravity_y,
    double time_step,
    std::vector<double2>& viscous_acceleration,
    std::vector<double2>& splash_repulsive_acceleration) {
  int num_particles = fluid_particles.particle_num;
  
  // 确保输出向量大小正确
  viscous_acceleration.resize(num_particles);
  splash_repulsive_acceleration.resize(num_particles);
  
  // 计算重力加速度（对所有粒子相同）
  double2 gravity_acceleration = ComputeGravityAcceleration(gravity_x, gravity_y);
  
  // 第一步：先计算所有粒子的速度拉普拉斯算子和粘性力加速度（使用原始速度）
  for (int i = 0; i < num_particles; ++i) {
    if (i < static_cast<int>(fluid_particles.surface_type.size()) &&
        fluid_particles.surface_type[i] == SurfaceType::SPLASH) {
      viscous_acceleration[i] = {0.0, 0.0};
      continue;
    }
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
    splash_repulsive_acceleration[i] = {0.0, 0.0};
    double2 non_pressure_acceleration = viscous_acceleration[i];
    if (i < static_cast<int>(fluid_particles.surface_type.size()) &&
        fluid_particles.surface_type[i] == SurfaceType::SPLASH) {
      splash_repulsive_acceleration[i] = ComputeSplashRepulsiveAcceleration(
          i, fluid_particles, solid_particles, particle_spacing);
      non_pressure_acceleration = splash_repulsive_acceleration[i];
    }
    UpdateVelocity(
        i,
        fluid_particles,
        non_pressure_acceleration,
        gravity_acceleration,
        time_step);
  }
}

} // namespace mps2D
