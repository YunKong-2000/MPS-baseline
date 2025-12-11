#pragma once
#include "../core/Particle.hpp"
#include "core/Types.h"
#include "core/MPSUtils.h"
#include <vector>
#include <array>
#include <cmath>

namespace mps {

// 自由面粒子判定器
// 功能：根据粒子数密度和几何方法判定流体粒子的自由面类型
class SurfaceDetector {
public:
  SurfaceDetector() = default;
  ~SurfaceDetector() = default;

  // 判定流体粒子的自由面类型
  // 参数：
  //   fluid_particles: 流体粒子对象（需要已构建邻域列表）
  //   solid_particles: 固体粒子对象
  //   smoothing_radius: 平滑半径（r_e）
  //   particle_radius: 粒子半径
  //   use_virtual_light: 是否使用虚拟光源法（true）或锥形区域法（false），默认true
  void DetectSurfaceParticles(FluidParticle& fluid_particles,
                              const SolidParticle& solid_particles,
                              const double smoothing_radius,
                              const double particle_radius,
                              bool use_virtual_light = true);

  // 虚拟光源法：计算球面幕布上的阴影面积比例（公开方法，用于输出到VTK）
  // 参数：
  //   particle_idx: 粒子索引
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象
  //   smoothing_radius: 平滑半径（球面幕布半径）
  //   particle_radius: 粒子半径（用于计算遮挡角度）
  // 返回：阴影面积占总球面面积的比例（0.0-1.0）
  double ComputeShadowAreaRatio(int particle_idx,
                                const FluidParticle& fluid_particles,
                                const SolidParticle& solid_particles,
                                const double smoothing_radius,
                                const double particle_radius) const;

private:
  // 计算粒子数密度
  // 参数：
  //   particle_idx: 粒子索引
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象
  //   smoothing_radius: 平滑半径
  double ComputeParticleNumberDensity(int particle_idx,
                                      const FluidParticle& fluid_particles,
                                      const SolidParticle& solid_particles,
                                      const double smoothing_radius) const;

  // 计算偏置向量（指向最稀疏方向的向量）
  // 参数：
  //   particle_idx: 粒子索引
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象
  //   smoothing_radius: 平滑半径
  // 返回：偏置向量（归一化）
  double3 ComputeBiasVector(int particle_idx,
                            const FluidParticle& fluid_particles,
                            const SolidParticle& solid_particles,
                            const double smoothing_radius) const;

  // 判定是否为飞溅粒子（基于粒子数密度）
  // 参数：
  //   number_density: 粒子数密度
  //   reference_density: 参考粒子数密度（通常为内部粒子的平均密度）
  bool IsSplashParticle(double number_density, double reference_density) const;

  // 判定是否为自由面粒子（基于几何方法）
  // 参数：
  //   particle_idx: 粒子索引
  //   bias_vector: 偏置向量（归一化）
  //   bias_magnitude: 偏置向量模长
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象
  //   smoothing_radius: 平滑半径
  //   bias_magnitude_threshold: 偏置向量模长阈值
  //   use_virtual_light: 是否使用虚拟光源法（true）或锥形区域法（false）
  //   particle_radius: 粒子半径（虚拟光源法需要）
  bool IsSurfaceParticle(int particle_idx,
                         const double3& bias_vector,
                         double bias_magnitude,
                         const FluidParticle& fluid_particles,
                         const SolidParticle& solid_particles,
                         const double smoothing_radius,
                         double bias_magnitude_threshold,
                         bool use_virtual_light = false,
                         double particle_radius = 0.0) const;

  // 在锥形区域内搜索邻域粒子数量
  // 参数：
  //   particle_idx: 粒子索引
  //   bias_vector: 偏置向量（归一化，指向锥形区域方向）
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象
  //   smoothing_radius: 平滑半径
  // 返回：锥形区域内的邻域粒子数量（包括流体和固体粒子）
  int CountParticlesInCone(int particle_idx,
                           const double3& bias_vector,
                           const FluidParticle& fluid_particles,
                           const SolidParticle& solid_particles,
                           const double smoothing_radius) const;

  // 判定是否为近自由面粒子
  // 参数：
  //   particle_idx: 粒子索引
  //   fluid_particles: 流体粒子对象
  //   smoothing_radius: 平滑半径
  bool IsNearSurfaceParticle(int particle_idx,
                             const FluidParticle& fluid_particles,
                             const double smoothing_radius) const;

  // 计算参考粒子数密度（内部粒子的平均密度）
  // 参数：
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象
  //   smoothing_radius: 平滑半径
  double ComputeReferenceDensity(const FluidParticle& fluid_particles,
                                 const SolidParticle& solid_particles,
                                 const double smoothing_radius) const;

  // 常量定义（避免魔法数字）
  static constexpr double SPLASH_DENSITY_RATIO = 0.3;  // 飞溅粒子密度阈值比例
  static constexpr double CONE_ANGLE_DEGREES = 45.0;   // 锥形区域角度（度）
  static constexpr double CONE_COS_THRESHOLD = 0.7071067811865476;  // cos(45°) 锥形区域余弦阈值
  static constexpr double BIAS_MAGNITUDE_THRESHOLD_RATIO = 0.4;  // 偏置向量模长阈值比例（相对于粒子间距）
  static constexpr int CONE_PARTICLE_THRESHOLD = 2;  // 锥形区域内粒子数量阈值
  static constexpr double NEAR_SURFACE_DISTANCE_RATIO = 1.5;  // 近自由面距离比例（相对于smoothing_radius）
  
  // 虚拟光源法相关常量
  static constexpr int SPHERE_GRID_THETA_RESOLUTION = 36;  // 球面网格theta方向分辨率（0到π）
  static constexpr int SPHERE_GRID_PHI_RESOLUTION = 72;   // 球面网格phi方向分辨率（0到2π）
  static constexpr double SHADOW_AREA_THRESHOLD = 0.88;    // 阴影面积比例阈值（大于等于此值判定为内部粒子）
  static constexpr int MIN_NEIGHBOR_COUNT_FOR_INNER = 8;  // 最小邻域粒子数（用于区分角落/棱边粒子和自由面粒子）
  static constexpr double CORNER_EDGE_SHADOW_THRESHOLD = 0.60;  // 角落/棱边粒子阴影面积上限
  static constexpr double FREE_SURFACE_SHADOW_THRESHOLD = 0.65;  // 自由面粒子阴影面积上限
  static constexpr double MAX_RADIUS_RATIO = 0.15;         // 最大半径比（用于限制遮挡角度）
  static constexpr double MAX_OCCLUSION_ANGLE_RATIO = 0.15;  // 最大遮挡角度比例（相对于π）
};

} // namespace mps
