#include "SurfaceDetector.hpp"
#include <algorithm>
#include <limits>
#include <cmath>

namespace mps2D {

void SurfaceDetector::DetectSurfaceParticles(
    FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const double smoothing_radius,
    const double particle_spacing) {
  
  const int num_particles = fluid_particles.particle_num;
  if (num_particles == 0) {
    return;
  }

  // 步骤1：计算所有粒子的粒子数密度
  std::vector<double> number_densities(num_particles);
  for (int i = 0; i < num_particles; ++i) {
    number_densities[i] = ComputeParticleNumberDensity(
        i, fluid_particles, solid_particles, smoothing_radius);
  }

  // 步骤2：计算参考粒子数密度（用于判定飞溅粒子）
  double reference_density = ComputeReferenceDensity(
      fluid_particles, solid_particles, smoothing_radius);

  // 初始化所有粒子为待判定状态
  for (int i = 0; i < num_particles; ++i) {
    fluid_particles.surface_type[i] = SurfaceType::INNER;  // 默认值，后续会更新
  }

  // 步骤3：粗筛
  // 粒子数密度特别小或邻域内流体粒子数很少的直接判定为飞溅粒子
  // 其余粒子进入细筛（虚拟光源法）判定 SURFACE/INNER
  // 其余为待细筛粒子
  std::vector<bool> need_fine_screening(num_particles, false);
  for (int i = 0; i < num_particles; ++i) {
    // SPLASH 只参考“邻域内流体粒子数”
    const int fluid_neighbors =
        static_cast<int>(fluid_particles.fluid_neighbour_list[i].size());
    
    if (!CoarseScreening(i, number_densities[i], reference_density, 
                          fluid_neighbors, fluid_particles)) {
      // 需要细筛
      need_fine_screening[i] = true;
    }
  }

  // 步骤4：细筛 - 使用虚拟光源法进行进一步判定
  for (int i = 0; i < num_particles; ++i) {
    if (need_fine_screening[i]) {
      if (FineScreening(i, fluid_particles, solid_particles, 
                       smoothing_radius, particle_spacing)) {
        fluid_particles.surface_type[i] = SurfaceType::SURFACE;
      } else {
        fluid_particles.surface_type[i] = SurfaceType::INNER;
      }
    }
  }

  // 步骤5：近自由面判定
  // 内部粒子中，1.5倍粒子间距内有自由面粒子的为近自由面粒子
  for (int i = 0; i < num_particles; ++i) {
    if (fluid_particles.surface_type[i] == SurfaceType::INNER) {
      if (IsNearSurfaceParticle(i, fluid_particles, particle_spacing)) {
        fluid_particles.surface_type[i] = SurfaceType::NEAR_SURFACE;
      }
    }
  }
}

double SurfaceDetector::ComputeParticleNumberDensity(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const double smoothing_radius) const {
  
  double density = 0.0;
  const double2& pos_i = fluid_particles.position[particle_idx];

  // 计算流体粒子贡献
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    double dist = ComputeDistance(pos_i, fluid_particles.position[j]);
    density += WeightFunction(dist, smoothing_radius);
  }

  // 计算固体粒子贡献
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    double dist = ComputeDistance(pos_i, solid_particles.position[j]);
    density += WeightFunction(dist, smoothing_radius);
  }

  return density;
}

double SurfaceDetector::ComputeReferenceDensity(
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const double smoothing_radius) const {
  
  // 计算所有粒子的平均粒子数密度作为参考
  double total_density = 0.0;
  int count = 0;

  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    double density = ComputeParticleNumberDensity(
        i, fluid_particles, solid_particles, smoothing_radius);
    total_density += density;
    ++count;
  }

  if (count > 0) {
    return total_density / count;
  } else {
    // 如果没有粒子，返回一个默认值
    return 1.0;
  }
}

bool SurfaceDetector::CoarseScreening(int particle_idx,
                                     double number_density,
                                     double reference_density,
                                     int neighbor_count,
                                     FluidParticle& fluid_particles) const {
  
  // 判定飞溅粒子：粒子数密度特别小或邻域粒子数很少
  if (number_density < (reference_density * SPLASH_DENSITY_RATIO) ||
      neighbor_count <= MIN_NEIGHBOR_COUNT_FOR_SPLASH) {
    fluid_particles.surface_type[particle_idx] = SurfaceType::SPLASH;
    return true;  // 已判定
  }
  
  // 不在粗筛阶段直接判定为内部粒子，避免“邻域足够大就直接内点”；
  // 其余粒子交给细筛（虚拟光源法）做最终判断。
  return false;
}

bool SurfaceDetector::FineScreening(int particle_idx,
                                   const FluidParticle& fluid_particles,
                                   const SolidParticle& solid_particles,
                                   const double smoothing_radius,
                                   const double particle_spacing) const {
  
  // 使用虚拟光源法：计算圆形幕布上的阴影面积比例
  double shadow_ratio = ComputeShadowAreaRatio(
      particle_idx, fluid_particles, solid_particles, smoothing_radius, particle_spacing);
  
  // 阴影面积比例小于阈值，判定为自由面粒子
  // 内部粒子周围有很多邻域粒子，阴影面积大（接近1.0）
  // 自由面粒子某个方向缺少邻域粒子，阴影面积小（接近0.0）
  return shadow_ratio < SHADOW_AREA_THRESHOLD;
}

bool SurfaceDetector::IsNearSurfaceParticle(int particle_idx,
                                           const FluidParticle& fluid_particles,
                                           const double particle_spacing) const {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  const double near_surface_distance = particle_spacing * NEAR_SURFACE_DISTANCE_RATIO;

  // 检查邻域内是否有自由面粒子
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    if (fluid_particles.surface_type[j] == SurfaceType::SURFACE) {
      double dist = ComputeDistance(pos_i, fluid_particles.position[j]);
      if (dist <= near_surface_distance) {
        return true;
      }
    }
  }

  return false;
}

// 辅助函数：标记圆形幕布上的阴影区域（2D版本）
namespace {
void mark_shadow_circle(std::vector<bool>& shadow_circle,
                       const double2& ri,
                       const double2& rj,
                       const double particle_radius,
                       const double smoothing_radius,
                       const int num_angles,
                       const double dangle) {
  // 计算粒子间距离
  const double2 r_ij = {rj.x - ri.x, rj.y - ri.y};
  const double rr = ComputeVectorMagnitude(r_ij);
  
  // 检查距离是否有效
  if (rr < 1e-10) {
    return;
  }
  
  // 计算角度（atan2返回[-π, π]，需要转换到[0, 2π]）
  double angle = std::atan2(r_ij.y, r_ij.x);
  // atan2返回值范围：[-π, π]
  // 需要转换到[0, 2π]以便用于角度网格索引
  if (angle < 0.0) {
    angle += 2.0 * M_PI;
  }
  
  // 计算角度对应的网格索引
  int angle_idx_center = static_cast<int>(angle / dangle);
  if (angle_idx_center < 0) angle_idx_center = 0;
  if (angle_idx_center >= num_angles) angle_idx_center = num_angles - 1;
  
  // 计算遮挡角度（2D版本：粒子在圆形幕布上形成的阴影角度）
  const double radius_ratio = particle_radius / rr;
  double occlusion_angle;
  if (radius_ratio >= 1.0) {
    // 如果粒子间距小于粒子半径，将目标粒子视为平行光源
    occlusion_angle = std::asin(particle_radius / smoothing_radius);
  } else {
    occlusion_angle = std::asin(radius_ratio);
  }
  
  // 计算遮挡范围对应的网格索引范围
  const int angle_range = static_cast<int>(std::ceil(occlusion_angle / dangle));
  
  // 标记被遮挡的网格
  for (int da = -angle_range; da <= angle_range; ++da) {
    int angle_idx = angle_idx_center + da;
    // 处理角度索引的周期性边界
    if (angle_idx < 0) angle_idx += num_angles;
    if (angle_idx >= num_angles) angle_idx -= num_angles;
    
    if (angle_idx >= 0 && angle_idx < num_angles) {
      shadow_circle[angle_idx] = true;
    }
  }
}
} // 匿名命名空间

double SurfaceDetector::ComputeShadowAreaRatio(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const double smoothing_radius,
    const double particle_spacing) const {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  
  // 步骤1：将圆形幕布划分为角度网格
  // 2D版本使用圆形幕布，只需要角度theta (0到2π)
  const int num_angles = CIRCLE_GRID_RESOLUTION;
  const double dangle = 2.0 * M_PI / num_angles;
  
  // 计算粒子半径（相对于粒子间距）
  const double particle_radius = particle_spacing * PARTICLE_RADIUS_RATIO;
  
  // 标记哪些角度网格被遮挡（使用一维数组）
  std::vector<bool> shadow_circle(num_angles, false);
  
  // 步骤2：遍历邻域粒子，标记被遮挡的角度网格
  // 处理流体粒子
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];
    const double dist = ComputeDistance(pos_i, pos_j);
    if (dist > 1e-10 && dist <= smoothing_radius) {
      mark_shadow_circle(shadow_circle, pos_i, pos_j, particle_radius,
                        smoothing_radius, num_angles, dangle);
    }
  }
  
  // 处理固体粒子
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    const double2& pos_j = solid_particles.position[j];
    const double dist = ComputeDistance(pos_i, pos_j);
    if (dist > 1e-10 && dist <= smoothing_radius) {
      mark_shadow_circle(shadow_circle, pos_i, pos_j, particle_radius,
                        smoothing_radius, num_angles, dangle);
    }
  }
  
  // 步骤3：统计被遮挡网格的总面积
  // 在2D圆形幕布上，每个角度网格的面积相等（dangle * R）
  // 但比例计算时R会约掉，所以只需要统计被遮挡的角度数量
  int shadowed_count = 0;
  for (int i = 0; i < num_angles; ++i) {
    if (shadow_circle[i]) {
      shadowed_count++;
    }
  }
  
  // 返回阴影面积比例（被遮挡角度数 / 总角度数）
  if (num_angles > 0) {
    return static_cast<double>(shadowed_count) / num_angles;
  } else {
    return 0.0;
  }
}

} // namespace mps2D

