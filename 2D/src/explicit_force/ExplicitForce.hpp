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

  // 更新单个粒子的速度（只更新速度，不更新位置）
  // 参数：
  //   particle_idx: 粒子索引
  //   fluid_particles: 流体粒子对象（会被修改）
  //   viscous_acceleration: 粘性力加速度
  //   gravity_acceleration: 重力加速度
  //   time_step: 时间步长
  // 注意：此方法只更新速度作为临时速度，位置保持不变
  void UpdateVelocity(
      int particle_idx,
      FluidParticle& fluid_particles,
      const double2& viscous_acceleration,
      const double2& gravity_acceleration,
      double time_step);

  // 批量计算并更新所有粒子的速度（只更新速度，不更新位置）
  // 参数：
  //   fluid_particles: 流体粒子对象（会被修改）
  //   solid_particles: 固体粒子对象
  //   corrective_matrices: 所有粒子的corrective matrix
  //   smoothing_radius: 平滑半径
  //   kinematic_viscosity: 动力学粘性系数
  //   gravity_x: 重力加速度x分量
  //   gravity_y: 重力加速度y分量
  //   time_step: 时间步长
  // 注意：此方法只更新速度作为临时速度，位置保持不变
  void ComputeAndUpdateVelocity(
      FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      const std::vector<Eigen::Matrix<double, 5, 5>>& corrective_matrices,
      double smoothing_radius,
      double particle_spacing,
      double kinematic_viscosity,
      double gravity_x,
      double gravity_y,
      double time_step);

  // 批量计算并更新所有粒子的速度（只更新速度，不更新位置），同时返回粘性力加速度
  // 参数：
  //   fluid_particles: 流体粒子对象（会被修改）
  //   solid_particles: 固体粒子对象
  //   corrective_matrices: 所有粒子的corrective matrix
  //   smoothing_radius: 平滑半径
  //   kinematic_viscosity: 动力学粘性系数
  //   gravity_x: 重力加速度x分量
  //   gravity_y: 重力加速度y分量
  //   time_step: 时间步长
  //   viscous_acceleration: 输出参数，保存每个粒子的粘性力加速度
  // 注意：此方法只更新速度作为临时速度，位置保持不变
  void ComputeAndUpdateVelocity(
      FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      const std::vector<Eigen::Matrix<double, 5, 5>>& corrective_matrices,
      double smoothing_radius,
      double particle_spacing,
      double kinematic_viscosity,
      double gravity_x,
      double gravity_y,
      double time_step,
      std::vector<double2>& viscous_acceleration);

  // 批量计算并更新所有粒子的速度（只更新速度，不更新位置），
  // 同时返回粘性力加速度与飞溅粒子排斥加速度
  // 参数：
  //   splash_repulsive_acceleration: 输出参数，保存每个粒子的排斥加速度；
  //                                 非飞溅粒子写为零向量
  void ComputeAndUpdateVelocity(
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
      std::vector<double2>& splash_repulsive_acceleration);

private:
  // 根据飞溅粒子与壁面粒子的距离，计算弹性排斥加速度修正
  double2 ComputeSplashRepulsiveAcceleration(
      int particle_idx,
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      double particle_spacing) const;

};

} // namespace mps2D
