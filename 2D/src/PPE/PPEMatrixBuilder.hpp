#pragma once
#include "../core/Particle.hpp"
#include "../lsmps/CorrectiveMatrix.hpp"
#include "core/Types.h"
#include "core/MPSUtils.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>

namespace mps2D {

// PPE（压力泊松方程）系数矩阵构建器
// 基于LSMPS算子离散方法构建压力泊松方程的系数矩阵A和右边项b
class PPEMatrixBuilder {
public:
  PPEMatrixBuilder() = default;
  ~PPEMatrixBuilder() = default;

  // 构建PPE系数矩阵A和右边项b
  // 参数：
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象
  //   corrective_matrices: 每个粒子的corrective matrix（5x5）
  //   smoothing_radius: 平滑半径（r_e）
  //   density: 流体密度（ρ）
  //   time_step: 时间步长（Δt）
  //   gravity_x: 重力加速度x分量（用于壁面压力边界条件）
  //   gravity_y: 重力加速度y分量（用于壁面压力边界条件）
  //   A: 输出的系数矩阵（稀疏矩阵）
  //   b: 输出的右边项向量
  // 返回：是否成功构建
  bool BuildPPEMatrix(
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      const std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>>& corrective_matrices,
      double smoothing_radius,
      double density,
      double time_step,
      double gravity_x,
      double gravity_y,
      Eigen::SparseMatrix<double, Eigen::ColMajor>& A,
      Eigen::VectorXd& b);

  // 统计非零元素个数（用于预分配稀疏矩阵内存）
  // 参数：
  //   fluid_particles: 流体粒子对象
  // 返回：非零元素个数
  int CountNonZeros(const FluidParticle& fluid_particles) const;

  // 调试信息结构体
  struct DebugInfo {
    Eigen::VectorXd diagonal_coefficients;      // 对角线系数
    Eigen::VectorXd off_diagonal_row_sums;       // 每行非对角线系数之和
    Eigen::VectorXd divergence_terms;            // 速度散度项（未除以时间步长）
    Eigen::VectorXd wall_pressure_terms;        // 壁面压力边界条件项（未乘以系数因子）
    Eigen::VectorXd right_hand_side;            // 完整的右边项
  };

  // 构建PPE系数矩阵A和右边项b（带调试信息）
  // 参数：
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象
  //   corrective_matrices: 每个粒子的corrective matrix（5x5）
  //   smoothing_radius: 平滑半径（r_e）
  //   density: 流体密度（ρ）
  //   time_step: 时间步长（Δt）
  //   gravity_x: 重力加速度x分量（用于壁面压力边界条件）
  //   gravity_y: 重力加速度y分量（用于壁面压力边界条件）
  //   A: 输出的系数矩阵（稀疏矩阵）
  //   b: 输出的右边项向量
  //   debug_info: 输出的调试信息
  // 返回：是否成功构建
  bool BuildPPEMatrixWithDebug(
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      const std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>>& corrective_matrices,
      double smoothing_radius,
      double density,
      double time_step,
      double gravity_x,
      double gravity_y,
      Eigen::SparseMatrix<double, Eigen::ColMajor>& A,
      Eigen::VectorXd& b,
      DebugInfo& debug_info);

private:
  // 为单个粒子构建系数矩阵行和右边项
  // 参数：
  //   particle_idx: 粒子索引
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象
  //   corrective_matrix: 该粒子的corrective matrix
  //   smoothing_radius: 平滑半径
  //   density: 流体密度
  //   time_step: 时间步长
  //   gravity_x: 重力加速度x分量
  //   gravity_y: 重力加速度y分量
  //   triplets: 输出的三元组列表（用于构建稀疏矩阵）
  //   b: 输出的右边项向量
  // 说明：此函数在一个循环中同时计算对角线系数、非对角线系数和右边项，避免重复计算
  void BuildParticleRow(
      int particle_idx,
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      const Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>& corrective_matrix,
      double smoothing_radius,
      double density,
      double time_step,
      double gravity_x,
      double gravity_y,
      std::vector<Eigen::Triplet<double>>& triplets,
      Eigen::VectorXd& b) const;

  // 为单个粒子构建系数矩阵行和右边项（带调试信息）
  // 参数：与BuildParticleRow相同，额外添加调试信息输出
  //   debug_info: 输出的调试信息
  void BuildParticleRowWithDebug(
      int particle_idx,
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      const Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>& corrective_matrix,
      double smoothing_radius,
      double density,
      double time_step,
      double gravity_x,
      double gravity_y,
      std::vector<Eigen::Triplet<double>>& triplets,
      Eigen::VectorXd& b,
      DebugInfo& debug_info) const;

};

} // namespace mps2D
