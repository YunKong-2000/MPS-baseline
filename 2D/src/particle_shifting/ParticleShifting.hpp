#pragma once

#include "../core/Particle.hpp"
#include "../lsmps/CorrectiveMatrix.hpp"
#include "core/MPSUtils.h"
#include <Eigen/Dense>
#include <vector>

namespace mps2D {

class ParticleShifting {
 public:
  ParticleShifting() = default;
  ~ParticleShifting() = default;

  void Apply(
      FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      const std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE,
                                      CorrectiveMatrix::MATRIX_SIZE>>&
          corrective_matrices_velocity,
      double smoothing_radius,
      double particle_spacing,
      std::vector<double2>* out_ps_displacement = nullptr,
      std::vector<double>* out_rve = nullptr,
      std::vector<double>* out_temp_number_density = nullptr,
      // PS 前（即更新速度/位置之前）的自由面法向量（非自由面为零向量）
      std::vector<double2>* out_free_surface_normals = nullptr,
      // PS 前的速度梯度张量分量（二维矩阵的 4 个标量分量）
      std::vector<double>* out_velocity_grad_xx = nullptr,
      std::vector<double>* out_velocity_grad_xy = nullptr,
      std::vector<double>* out_velocity_grad_yx = nullptr,
      std::vector<double>* out_velocity_grad_yy = nullptr) const;

 private:
  static constexpr int kDimension = 2;
  static constexpr double kLambdaShift = 0.1;
  static constexpr double kMaxShiftFactor = 0.1;
  static constexpr double kSafeDistanceFactor = 0.25;
  static constexpr double kEps = 1e-12;

  double ComputeReferenceNumberDensity(
      double particle_spacing,
      double re_ps) const;

  void ComputeMinWallDistanceAndRve(
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      double smoothing_radius,
      double particle_spacing,
      std::vector<double>& rve) const;

  void ComputeFreeSurfaceNormals(
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      const std::vector<double>& rve,
      std::vector<double2>& normals) const;

  double2 ComputeBasicShift(
      int particle_idx,
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      double rve,
      double l0,
      double n0) const;

  // 根据文档 Liu_shifting 的公式计算滑移矢量 δr^SL
  // 切向方向由速度相对项在切向投影后的单位向量给出。
  double2 ComputeSlipShift(
      int particle_idx,
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      double rve,
      double l0) const;

  double2 ProjectForSurfaceConstraint(
      const double2& shift,
      const double2& normal) const;

  double2 ApplyMagnitudeLimiter(
      const double2& shift,
      double particle_spacing) const;

  double2 ApplyWallCollisionLimiter(
      int particle_idx,
      const double2& shift,
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      double particle_spacing) const;

  Eigen::Matrix2d ComputeVelocityGradient(
      int particle_idx,
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      const Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE,
                          CorrectiveMatrix::MATRIX_SIZE>& corrective_matrix,
      double smoothing_radius) const;
};

}  // namespace mps2D

