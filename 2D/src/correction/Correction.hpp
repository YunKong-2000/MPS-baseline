#pragma once
#include "../core/Particle.hpp"
#include "../lsmps/CorrectiveMatrix.hpp"
#include "core/Types.h"
#include "core/MPSUtils.h"
#include <Eigen/Dense>
#include <vector>

namespace mps2D {

// 压力修正模块
// 功能：根据PPE求解出的压力计算压力梯度，并根据压力梯度计算粒子的加速度，
//       从而进一步矫正粒子的速度和位置
// 
// 算法步骤：
// 1. 根据输入的压力计算压力梯度，所有流体粒子使用LSMPS方法来离散压力梯度算子
//    计算压力梯度，壁面边界条件为第二类边界条件（dp/dn = ρ * n · g）
// 2. 计算粒子加速度，也即是 -(压力梯度/密度)
// 3. 根据算得的加速度和时间步长更新粒子速度和位置
class Correction {
public:
  Correction() = default;
  ~Correction() = default;

  // 计算单个粒子的压力梯度（使用LSMPS方法，第二类边界条件）
  // 参数：
  //   particle_idx: 粒子索引
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象
  //   corrective_matrix: corrective matrix (5x5)，使用第二类边界条件计算得到
  //   smoothing_radius: 平滑半径（r_e）
  //   gravity_x: 重力加速度x分量（用于壁面边界条件）
  //   gravity_y: 重力加速度y分量（用于壁面边界条件）
  //   density: 流体密度（用于壁面边界条件）
  // 返回：压力梯度向量 (grad_x, grad_y)
  double2 ComputePressureGradient(
      int particle_idx,
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      const Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>& corrective_matrix,
      double smoothing_radius,
      double gravity_x,
      double gravity_y,
      double density);

  // 计算压力梯度引起的加速度
  // 参数：
  //   pressure_gradient: 压力梯度向量
  //   density: 流体密度
  // 返回：加速度向量 (a_x, a_y) = -(压力梯度/密度)
  double2 ComputeAcceleration(
      const double2& pressure_gradient,
      double density);

  // 更新单个粒子的速度和位置（根据压力梯度引起的加速度）
  // 参数：
  //   particle_idx: 粒子索引
  //   fluid_particles: 流体粒子对象（会被修改）
  //   acceleration: 压力梯度引起的加速度
  //   time_step: 时间步长
  void UpdateVelocityAndPosition(
      int particle_idx,
      FluidParticle& fluid_particles,
      const double2& acceleration,
      double time_step);

  // 批量计算并更新所有粒子的速度和位置
  // 参数：
  //   fluid_particles: 流体粒子对象（会被修改）
  //   solid_particles: 固体粒子对象
  //   corrective_matrices: 所有粒子的corrective matrix（使用第二类边界条件计算得到）
  //   smoothing_radius: 平滑半径
  //   gravity_x: 重力加速度x分量
  //   gravity_y: 重力加速度y分量
  //   density: 流体密度
  //   time_step: 时间步长
  void ComputeAndUpdateAllParticles(
      FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      const std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>>& corrective_matrices,
      double smoothing_radius,
      double gravity_x,
      double gravity_y,
      double density,
      double time_step);

  // 批量计算并更新所有粒子的速度和位置，同时保存压力梯度（在更新位置之前）
  // 参数：
  //   fluid_particles: 流体粒子对象（会被修改）
  //   solid_particles: 固体粒子对象
  //   corrective_matrices: 所有粒子的corrective matrix（使用第二类边界条件计算得到）
  //   smoothing_radius: 平滑半径
  //   gravity_x: 重力加速度x分量
  //   gravity_y: 重力加速度y分量
  //   density: 流体密度
  //   time_step: 时间步长
  //   pressure_gradients: 输出参数，保存所有粒子的压力梯度（在位置更新之前计算的）
  void ComputeAndUpdateAllParticles(
      FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      const std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>>& corrective_matrices,
      double smoothing_radius,
      double gravity_x,
      double gravity_y,
      double density,
      double time_step,
      std::vector<double2>& pressure_gradients);
};

} // namespace mps2D

