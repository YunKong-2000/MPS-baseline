#pragma once
#include "../core/Particle.hpp"
#include "core/Types.h"
#include "core/MPSUtils.h"
#include <Eigen/Dense>
#include <Eigen/Core>

namespace mps2D {

// LSMPS moment matrix / corrective matrix 计算模块
// 用于恢复MPS方法中离散化的一致性，提高计算精度
// 2D moment/corrective 矩阵是5x5的矩阵
class CorrectiveMatrix {
public:
  // 矩阵大小（2D情况：5x5）
  static constexpr int MATRIX_SIZE = 5;
  
  // 基函数数量
  static constexpr int BASIS_SIZE = 5;

  CorrectiveMatrix() = default;
  ~CorrectiveMatrix() = default;

  // 兼容旧接口：为单个粒子计算corrective matrix
  // 内部直接调用ComputeMomentMatrix（两者数值等价）
  Eigen::Matrix<double, MATRIX_SIZE, MATRIX_SIZE>
  ComputeCorrectiveMatrix(
      int particle_idx,
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      double smoothing_radius,
      bool border_condition = false);

  // 兼容旧接口：仅考虑流体粒子的corrective matrix
  Eigen::Matrix<double, MATRIX_SIZE, MATRIX_SIZE>
  ComputeCorrectiveMatrixFluidOnly(
      int particle_idx,
      const FluidParticle& fluid_particles,
      double smoothing_radius);

  // 为单个粒子计算moment matrix
  // 参数：
  //   particle_idx: 粒子索引
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象
  //   smoothing_radius: 平滑半径（r_e）
  //   border_condition: 边界条件类型，false表示第一类边界条件，true表示第二类边界条件
  // 返回：moment matrix（5x5），如果计算失败返回单位矩阵
  Eigen::Matrix<double, MATRIX_SIZE, MATRIX_SIZE>
  ComputeMomentMatrix(
      int particle_idx,
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      double smoothing_radius,
      bool border_condition = false);

  // 为单个粒子计算moment matrix（仅考虑流体粒子）
  // 该方法不考虑壁面粒子的影响，仅使用流体粒子之间的相互作用
  // 参数：
  //   particle_idx: 粒子索引
  //   fluid_particles: 流体粒子对象
  //   smoothing_radius: 平滑半径（r_e）
  // 返回：moment matrix（5x5），如果计算失败返回单位矩阵
  Eigen::Matrix<double, MATRIX_SIZE, MATRIX_SIZE>
  ComputeMomentMatrixFluidOnly(
      int particle_idx,
      const FluidParticle& fluid_particles,
      double smoothing_radius);

  // 计算基函数值（用于流体邻域粒子）
  // 参数：
  //   dx, dy: 相对于中心粒子的相对位置
  //   smoothing_radius: 平滑半径（用于归一化）
  // 返回：基函数值向量 [x/r_e, y/r_e, (x/r_e)^2, (y/r_e)^2, (x/r_e)*(y/r_e)]
  Eigen::Vector<double, BASIS_SIZE> ComputeBasisFunctions(
      double dx, double dy, double smoothing_radius);

  // 计算基函数值（用于固体邻域粒子/壁面粒子）
  // 参数：
  //   dx, dy: 相对于中心粒子的相对位置
  //   normal_x, normal_y: 壁面法向量的x和y分量
  //   smoothing_radius: 平滑半径（用于归一化）
  // 返回：基函数值向量 [n_x, n_y, 2*n_x*x/r_e, 2*n_y*y/r_e, (n_x*x + n_y*y)/r_e]
  Eigen::Vector<double, BASIS_SIZE> ComputeBasisFunctionsForWall(
      double dx, double dy,
      double normal_x, double normal_y,
      double smoothing_radius);

  // 诊断用：构建系数矩阵C（用于计算条件数等诊断信息）
  // 参数：
  //   particle_idx: 粒子索引
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象
  //   smoothing_radius: 平滑半径
  //   border_condition: 边界条件类型，false表示第一类边界条件，true表示第二类边界条件
  // 返回：系数矩阵C（5x5）
  Eigen::MatrixXd BuildCoefficientMatrixForDiagnostics(
      int particle_idx,
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      double smoothing_radius,
      bool border_condition = false);

private:
  // 构建系数矩阵C
  // 参数：
  //   particle_idx: 粒子索引
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象
  //   smoothing_radius: 平滑半径
  //   border_condition: 边界条件类型，false表示第一类边界条件，true表示第二类边界条件
  // 返回：系数矩阵C（5x5）
  Eigen::MatrixXd BuildCoefficientMatrix(
      int particle_idx,
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      double smoothing_radius,
      bool border_condition = false);

  // 检查矩阵是否可逆
  // 参数：
  //   matrix: 待检查的矩阵
  //   tolerance: 行列式阈值
  // 返回：如果矩阵可逆返回true，否则返回false
  bool IsMatrixInvertible(
      const Eigen::Matrix<double, MATRIX_SIZE, MATRIX_SIZE>& matrix,
      double tolerance = 1e-10) const;
};

} // namespace mps2D

