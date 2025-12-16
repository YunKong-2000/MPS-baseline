#pragma once
#include "../core/Particle.hpp"
#include "../lsmps/CorrectiveMatrix.hpp"
#include "core/Types.h"
#include "core/MPSUtils.h"
#include <Eigen/Dense>
#include <vector>

namespace mps2D {

// 显式力计算模块
// 功能：计算粘性力和重力，并根据这些力和当前速度及时间步大小更新位置和速度
class ExplicitForce {
public:
  ExplicitForce() = default;
  ~ExplicitForce() = default;

  // 计算速度拉普拉斯算子（使用LSMPS corrective matrix方法）
  // 参数：
  //   particle_idx: 粒子索引
  //   velocity_field: 速度场（向量场）
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象
  //   corrective_matrix: corrective matrix (5x5)
  //   smoothing_radius: 平滑半径
  // 返回：速度拉普拉斯算子（double2，包含x和y分量）
  // 注意：对于壁面粒子，统一使用第一类边界条件（无滑移边界，壁面速度为零）
  double2 ComputeVelocityLaplacian(
      int particle_idx,
      const std::vector<double2>& velocity_field,
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      const Eigen::Matrix<double, 5, 5>& corrective_matrix,
      double smoothing_radius);

  // 计算粘性力加速度
  // 参数：
  //   velocity_laplacian: 速度拉普拉斯算子
  //   kinematic_viscosity: 动力学粘性系数（ν）
  // 返回：粘性力加速度（double2）
  double2 ComputeViscousAcceleration(
      const double2& velocity_laplacian,
      double kinematic_viscosity);

  // 计算重力加速度
  // 参数：
  //   gravity_x: 重力加速度x分量
  //   gravity_y: 重力加速度y分量
  // 返回：重力加速度（double2）
  double2 ComputeGravityAcceleration(double gravity_x, double gravity_y);

  // 更新单个粒子的速度和位置（显式时间积分）
  // 参数：
  //   particle_idx: 粒子索引
  //   fluid_particles: 流体粒子对象（会被修改）
  //   viscous_acceleration: 粘性力加速度
  //   gravity_acceleration: 重力加速度
  //   time_step: 时间步长
  void UpdateVelocityAndPosition(
      int particle_idx,
      FluidParticle& fluid_particles,
      const double2& viscous_acceleration,
      const double2& gravity_acceleration,
      double time_step);

  // 批量计算并更新所有粒子的速度和位置
  // 参数：
  //   fluid_particles: 流体粒子对象（会被修改）
  //   solid_particles: 固体粒子对象
  //   corrective_matrices: 所有粒子的corrective matrix
  //   smoothing_radius: 平滑半径
  //   kinematic_viscosity: 动力学粘性系数
  //   gravity_x: 重力加速度x分量
  //   gravity_y: 重力加速度y分量
  //   time_step: 时间步长
  void ComputeAndUpdateAllParticles(
      FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      const std::vector<Eigen::Matrix<double, 5, 5>>& corrective_matrices,
      double smoothing_radius,
      double kinematic_viscosity,
      double gravity_x,
      double gravity_y,
      double time_step);
};

} // namespace mps2D
