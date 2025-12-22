#pragma once
#include "../core/Particle.hpp"
#include "../lsmps/CorrectiveMatrix.hpp"
#include "core/Types.h"
#include "core/MPSUtils.h"
#include <Eigen/Dense>
#include <vector>

// PETSc前向声明（避免在头文件中包含PETSc头文件）
// 注意：这些类型在PETSc中是typedef，这里使用前向声明
struct _p_Mat;
struct _p_Vec;
typedef struct _p_Mat* Mat;
typedef struct _p_Vec* Vec;

namespace mps2D {

// PPE（压力泊松方程）系数矩阵构建器
// 基于LSMPS算子离散方法构建压力泊松方程的系数矩阵A和右边项b（PETSc格式）
class PPEMatrixBuilder {
public:
  PPEMatrixBuilder() = default;
  ~PPEMatrixBuilder() = default;

  // 直接构建PETSc格式的PPE系数矩阵A和右边项b
  // 参数：
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象
  //   corrective_matrices: 每个粒子的corrective matrix（5x5）
  //   smoothing_radius: 平滑半径（r_e）
  //   density: 流体密度（ρ）
  //   time_step: 时间步长（Δt）
  //   particle_spacing: 粒子初始间距（l_0，用于自由面粒子计算）
  //   gravity_x: 重力加速度x分量（用于壁面压力边界条件）
  //   gravity_y: 重力加速度y分量（用于壁面压力边界条件）
  //   A_petsc: 输出的PETSc系数矩阵（必须在调用前初始化为NULL或已创建的Mat对象）
  //   b_petsc: 输出的PETSc右边项向量（必须在调用前初始化为NULL或已创建的Vec对象）
  // 返回：是否成功构建
  // 注意：调用者负责销毁返回的Mat和Vec对象（使用MatDestroy和VecDestroy）
  bool BuildPPEMatrixPetsc(
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      const std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>>& corrective_matrices,
      double smoothing_radius,
      double density,
      double time_step,
      double particle_spacing,
      double gravity_x,
      double gravity_y,
      Mat& A_petsc,
      Vec& b_petsc);

  // 调试函数：将系数矩阵的对角线元素和右边项输出到VTK文件
  // 参数：
  //   fluid_particles: 流体粒子对象（用于输出粒子位置）
  //   A_petsc: PETSc系数矩阵
  //   b_petsc: PETSc右边项向量
  //   filename: 输出VTK文件名
  // 返回：是否成功输出
  // 注意：VTK文件包含粒子位置、速度、对角线元素和右边项等信息
  bool WriteDebugInfoToVTK(
      const FluidParticle& fluid_particles,
      Mat A_petsc,
      Vec b_petsc,
      const std::string& filename) const;

private:
  // 初始化PETSc矩阵和向量
  void InitializePetscMatrixAndVector(
      int num_particles,
      const std::vector<int>& nnz_per_row,
      Mat& A_petsc,
      Vec& b_petsc) const;
  
  // 构建自由面粒子的矩阵行和右边项
  void BuildSurfaceParticleRow(
      int particle_idx,
      const FluidParticle& fluid_particles,
      double smoothing_radius,
      double density,
      double time_step,
      double reference_density,
      double lambda,
      double w_l0,
      Mat& A_petsc,
      Vec& b_petsc) const;
  
  // 构建内部粒子的矩阵行和右边项
  void BuildInnerParticleRow(
      int particle_idx,
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      const Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>& corrective_matrix,
      CorrectiveMatrix& corrective_matrix_calc,
      double smoothing_radius,
      double density,
      double time_step,
      double gravity_x,
      double gravity_y,
      double coeff_factor,
      Mat& A_petsc,
      Vec& b_petsc) const;
};

} // namespace mps2D
