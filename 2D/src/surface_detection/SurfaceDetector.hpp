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
  //   fluid_particles: 流体粒子对象
  //   smoothing_radius: 平滑半径
  // 返回：true表示已判定为飞溅粒子，false表示需要细筛
  bool CoarseScreening(int particle_idx,
                      FluidParticle& fluid_particles,
                      const double smoothing_radius) const;

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
  static constexpr int MIN_FLUID_NEIGHBOR_COUNT_FOR_SPLASH = 5;  // 流体邻域粒子数阈值（小于该值判定为飞溅）
  static constexpr double NEAR_SURFACE_DISTANCE_RATIO = 1.5;  // 近自由面距离比例（相对于粒子间距）
  
  // 虚拟光源法相关常量（2D版本使用圆形幕布）
  static constexpr int CIRCLE_GRID_RESOLUTION = 180;   // 圆形网格角度分辨率（0到2π）
  static constexpr double SHADOW_AREA_THRESHOLD = 0.9;    // 阴影面积比例阈值（大于等于此值判定为内部粒子）
  static constexpr double PARTICLE_RADIUS_RATIO = 0.5;  // 粒子半径相对于粒子间距的比例
};

} // namespace mps2D

