#pragma once
#include "../core/Particle.hpp"
#include "core/Types.h"
#include "core/MPSUtils.h"
#include <vector>
#include <array>
#include <cmath>

namespace mps2D {

// 自由面粒子判定器（2D版本）
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
  //   particle_spacing: 粒子间距
  void DetectSurfaceParticles(FluidParticle& fluid_particles,
                              const SolidParticle& solid_particles,
                              const double smoothing_radius,
                              const double particle_spacing);

  // 虚拟光源法：计算圆形幕布上的阴影面积比例（公开方法，用于输出到VTK）
  // 参数：
  //   particle_idx: 粒子索引
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象
  //   smoothing_radius: 平滑半径（圆形幕布半径）
  //   particle_spacing: 粒子间距（用于计算遮挡角度）
  // 返回：阴影面积占总圆形面积的比例（0.0-1.0）
  double ComputeShadowAreaRatio(int particle_idx,
                                const FluidParticle& fluid_particles,
                                const SolidParticle& solid_particles,
                                const double smoothing_radius,
                                const double particle_spacing) const;

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

  // 计算参考粒子数密度（内部粒子的平均密度）
  // 参数：
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象
  //   smoothing_radius: 平滑半径
  double ComputeReferenceDensity(const FluidParticle& fluid_particles,
                                 const SolidParticle& solid_particles,
                                 const double smoothing_radius) const;

  // 粗筛：判定飞溅粒子和内部粒子
  // 参数：
  //   particle_idx: 粒子索引
  //   number_density: 粒子数密度
  //   reference_density: 参考粒子数密度
  //   neighbor_count: 邻域粒子总数
  // 返回：true表示已判定（飞溅或内部），false表示需要细筛
  bool CoarseScreening(int particle_idx,
                      double number_density,
                      double reference_density,
                      int neighbor_count,
                      FluidParticle& fluid_particles) const;

  // 细筛：使用虚拟光源法判定自由面粒子
  // 参数：
  //   particle_idx: 粒子索引
  //   fluid_particles: 流体粒子对象
  //   solid_particles: 固体粒子对象
  //   smoothing_radius: 平滑半径
  //   particle_spacing: 粒子间距
  // 返回：true表示自由面粒子，false表示内部粒子
  bool FineScreening(int particle_idx,
                    const FluidParticle& fluid_particles,
                    const SolidParticle& solid_particles,
                    const double smoothing_radius,
                    const double particle_spacing) const;

  // 判定是否为近自由面粒子
  // 参数：
  //   particle_idx: 粒子索引
  //   fluid_particles: 流体粒子对象
  //   particle_spacing: 粒子间距
  bool IsNearSurfaceParticle(int particle_idx,
                             const FluidParticle& fluid_particles,
                             const double particle_spacing) const;

  // 常量定义（避免魔法数字）
  static constexpr double SPLASH_DENSITY_RATIO = 0.3;  // 飞溅粒子密度阈值比例
  static constexpr double INNER_DENSITY_RATIO = 0.97;  // 内部粒子密度阈值比例（接近参考密度）
  static constexpr int MIN_NEIGHBOR_COUNT_FOR_INNER = 8;  // 最小邻域粒子数（用于粗筛）
  static constexpr int MIN_NEIGHBOR_COUNT_FOR_SPLASH = 4;  // 邻域粒子数 < 5 判为飞溅（粗筛）
  static constexpr double NEAR_SURFACE_DISTANCE_RATIO = 1.5;  // 近自由面距离比例（相对于粒子间距）
  
  // 虚拟光源法相关常量（2D版本使用圆形幕布）
  static constexpr int CIRCLE_GRID_RESOLUTION = 180;   // 圆形网格角度分辨率（0到2π）
  static constexpr double SHADOW_AREA_THRESHOLD = 0.92;   // 阴影面积比例阈值（shadow_ratio < 该值判为自由面）
  static constexpr double PARTICLE_RADIUS_RATIO = 0.5;  // 粒子半径相对于粒子间距的比例
};

} // namespace mps2D

