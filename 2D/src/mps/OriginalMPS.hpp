#pragma once
#include "../core/Particle.hpp"
#include "core/Types.h"
#include "core/MPSUtils.h"
#include <vector>

namespace mps2D {

// 原始MPS方法算子计算模块
// 用于自由面粒子和近自由面粒子，这些粒子的邻域粒子数较少，无法使用LSMPS方法
// 
// 功能：
// 1. 计算标量场的梯度
// 2. 计算向量场的散度
// 3. 计算标量场的拉普拉斯算子
class OriginalMPS {
public:
  OriginalMPS() = default;
  ~OriginalMPS() = default;

  // 计算标量场的梯度（原始MPS方法）
  // 参数：
  //   particle_idx: 粒子索引
  //   scalar_field: 标量场值（如压力）
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象（可选，用于边界条件）
  //   smoothing_radius: 平滑半径（r_e）
  //   reference_density: 参考粒子数密度（n_0）
  // 返回：梯度向量 (grad_x, grad_y)
  double2 ComputeGradient(
      int particle_idx,
      const std::vector<double>& scalar_field,
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      double smoothing_radius,
      double reference_density) const;

  // 计算向量场的散度（原始MPS方法）
  // 参数：
  //   particle_idx: 粒子索引
  //   vector_field: 向量场值（如速度场）
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象（可选，用于边界条件）
  //   smoothing_radius: 平滑半径（r_e）
  //   reference_density: 参考粒子数密度（n_0）
  // 返回：散度值（标量）
  double ComputeDivergence(
      int particle_idx,
      const std::vector<double2>& vector_field,
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      double smoothing_radius,
      double reference_density) const;

  // 计算标量场的拉普拉斯算子（原始MPS方法）
  // 参数：
  //   particle_idx: 粒子索引
  //   scalar_field: 标量场值（如压力）
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象（可选，用于边界条件）
  //   smoothing_radius: 平滑半径（r_e）
  //   reference_density: 参考粒子数密度（n_0）
  //   particle_spacing: 粒子间距（用于计算离散lambda）
  // 返回：拉普拉斯算子值（标量）
  double ComputeLaplacian(
      int particle_idx,
      const std::vector<double>& scalar_field,
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      double smoothing_radius,
      double reference_density,
      double particle_spacing) const;

  // 计算参考粒子数密度（n_0）
  // 参数：
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象
  //   smoothing_radius: 平滑半径
  // 返回：参考粒子数密度（所有粒子的平均粒子数密度）
  double ComputeReferenceDensity(
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      double smoothing_radius) const;

  // 计算拉普拉斯归一化参数lambda（离散分布方法，对所有粒子都相同）
  // 直接生成均匀分布的粒子配置来计算，不使用现有粒子
  // 参数：
  //   particle_spacing: 粒子间距（与算例相同）
  //   smoothing_radius: 平滑半径（r_e，与算例相同）
  // 返回：lambda值，使用离散求和公式计算（对所有粒子都相同）
  double ComputeLambdaDiscreteUniform(
      double particle_spacing,
      double smoothing_radius) const;

  // 比较两种lambda计算方法的差值
  // 参数：
  //   particle_spacing: 粒子间距（与算例相同）
  //   smoothing_radius: 平滑半径（r_e，与算例相同）
  // 返回：差值统计信息
  struct LambdaComparisonResult {
    double analytical_lambda;      // 解析方法计算的lambda
    double discrete_lambda;        // 离散方法计算的lambda（使用均匀分布）
    double difference;              // 差值
    double relative_difference;    // 相对差值（百分比）
    int num_neighbors;             // 均匀分布中的邻域粒子数
  };
  
  LambdaComparisonResult CompareLambdaMethods(
      double particle_spacing,
      double smoothing_radius) const;

private:
  // 计算拉普拉斯归一化参数lambda（解析公式，对所有粒子都相同）
  // 参数：
  //   smoothing_radius: 平滑半径（r_e）
  // 返回：lambda值，计算公式为 λ = (1/5) * r_e^2
  double ComputeLambda(double smoothing_radius) const;

  // 计算拉普拉斯归一化参数lambda（离散分布方法，针对特定粒子）
  // 参数：
  //   particle_idx: 粒子索引
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象
  //   smoothing_radius: 平滑半径
  // 返回：lambda值，使用离散求和公式计算
  double ComputeLambdaDiscrete(
      int particle_idx,
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      double smoothing_radius) const;

  // 生成均匀分布的粒子配置并计算离散lambda
  // 参数：
  //   particle_spacing: 粒子间距
  //   smoothing_radius: 平滑半径
  //   center_pos: 中心粒子位置（默认原点）
  // 返回：离散lambda值
  double ComputeLambdaFromUniformDistribution(
      double particle_spacing,
      double smoothing_radius,
      const double2& center_pos = {0.0, 0.0}) const;

  // 计算粒子数密度（用于计算参考密度）
  // 参数：
  //   particle_idx: 粒子索引
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象
  //   smoothing_radius: 平滑半径
  // 返回：粒子数密度
  double ComputeParticleNumberDensity(
      int particle_idx,
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      double smoothing_radius) const;

  // 维度（2D）
  static constexpr int DIMENSION = 2;
};

} // namespace mps2D
