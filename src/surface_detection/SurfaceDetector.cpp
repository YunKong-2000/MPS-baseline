#include "SurfaceDetector.hpp"
#include <algorithm>
#include <limits>
#include <cmath>

namespace mps {

void SurfaceDetector::DetectSurfaceParticles(
    FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const double smoothing_radius,
    const double particle_radius,
    bool use_virtual_light) {
  
  // 计算偏置向量模长阈值（相对于粒子间距）
  // 粒子间距约为 particle_radius * 2
  double particle_spacing = particle_radius * 2.0;
  double bias_magnitude_threshold = particle_spacing * BIAS_MAGNITUDE_THRESHOLD_RATIO;
  
  const int num_particles = fluid_particles.particle_num;
  if (num_particles == 0) {
    return;
  }

  // 步骤1：计算所有粒子的粒子数密度和偏置向量
  std::vector<double> number_densities(num_particles);
  std::vector<double3> bias_vectors(num_particles);
  std::vector<double> bias_magnitudes(num_particles);

  for (int i = 0; i < num_particles; ++i) {
    number_densities[i] = ComputeParticleNumberDensity(
        i, fluid_particles, solid_particles, smoothing_radius);
    
    // 只调用一次计算偏置向量（未归一化）
    double3 bias_vec_raw = ComputeBiasVector(
        i, fluid_particles, solid_particles, smoothing_radius);
    
    // 从偏置向量计算模长和归一化向量
    bias_magnitudes[i] = mps::ComputeVectorMagnitude(bias_vec_raw);
    if (bias_magnitudes[i] > 1e-10) {
      bias_vectors[i] = mps::NormalizeVector(bias_vec_raw);
    } else {
      bias_vectors[i] = {0.0, 0.0, 0.0};
    }
  }

  // 步骤2：计算参考粒子数密度（用于判定飞溅粒子）
  double reference_density = ComputeReferenceDensity(
      fluid_particles, solid_particles, smoothing_radius);

  // 步骤3：初步判定飞溅粒子
  for (int i = 0; i < num_particles; ++i) {
    if (IsSplashParticle(number_densities[i], reference_density)) {
      fluid_particles.surface_type[i] = SurfaceType::SPLASH;
    } else {
      // 暂时标记为内部粒子，后续步骤会更新
      fluid_particles.surface_type[i] = SurfaceType::INNER;
    }
  }

  // 步骤4：根据几何方法判定内部粒子和自由面粒子
  for (int i = 0; i < num_particles; ++i) {
    // 跳过已判定的飞溅粒子
    if (fluid_particles.surface_type[i] == SurfaceType::SPLASH) {
      continue;
    }

    if (use_virtual_light) {
      // 使用虚拟光源法：直接基于阴影面积判定
      // 检查邻域粒子数量
      int total_neighbors = fluid_particles.fluid_neighbour_list[i].size() + 
                           fluid_particles.solid_neighbour_list[i].size();
      
      // 如果没有邻域粒子，直接判定为自由面粒子（飞溅粒子）
      if (total_neighbors == 0) {
        // 这种情况应该已经在步骤2中被判定为飞溅粒子了
        // 但为了安全起见，这里也检查一下
        if (fluid_particles.surface_type[i] != SurfaceType::SPLASH) {
          fluid_particles.surface_type[i] = SurfaceType::SURFACE;
        }
        continue;
      }
      
      // 计算阴影面积比例
      double shadow_ratio = ComputeShadowAreaRatio(
          i, fluid_particles, solid_particles, smoothing_radius, particle_radius);
      
      // 阴影面积比例小于阈值，判定为自由面粒子
      // 内部粒子周围有很多邻域粒子，阴影面积大（接近1.0）
      // 自由面粒子某个方向缺少邻域粒子，阴影面积小（接近0.0）
      if (shadow_ratio < SHADOW_AREA_THRESHOLD) {
        fluid_particles.surface_type[i] = SurfaceType::SURFACE;
      } else {
        fluid_particles.surface_type[i] = SurfaceType::INNER;
      }
    } else {
      // 使用锥形区域法：需要偏置向量模长检查
      // 根据偏置向量模长初步判断
      if (bias_magnitudes[i] < bias_magnitude_threshold) {
        // 偏置向量较小，判定为内部粒子
        fluid_particles.surface_type[i] = SurfaceType::INNER;
      } else {
        // 偏置向量较大，检查是否为自由面粒子
        if (IsSurfaceParticle(i, bias_vectors[i], bias_magnitudes[i],
                             fluid_particles, solid_particles, smoothing_radius,
                             bias_magnitude_threshold, use_virtual_light, particle_radius)) {
          fluid_particles.surface_type[i] = SurfaceType::SURFACE;
        } else {
          fluid_particles.surface_type[i] = SurfaceType::INNER;
        }
      }
    }
  }

  // 步骤5：判定近自由面粒子
  for (int i = 0; i < num_particles; ++i) {
    // 只对已判定为内部粒子的粒子进行检查
    if (fluid_particles.surface_type[i] == SurfaceType::INNER) {
      if (IsNearSurfaceParticle(i, fluid_particles, smoothing_radius)) {
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
  const double3& pos_i = fluid_particles.position[particle_idx];

  // 计算流体粒子贡献
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    double dist = mps::ComputeDistance(pos_i, fluid_particles.position[j]);
    density += mps::WeightFunction(dist, smoothing_radius);
  }

  // 计算固体粒子贡献
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    double dist = mps::ComputeDistance(pos_i, solid_particles.position[j]);
    density += mps::WeightFunction(dist, smoothing_radius);
  }

  return density;
}

double3 SurfaceDetector::ComputeBiasVector(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const double smoothing_radius) const {
  
  double3 weighted_center = {0.0, 0.0, 0.0};
  double total_weight = 0.0;
  const double3& pos_i = fluid_particles.position[particle_idx];

  // 计算流体粒子贡献
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    double dist = mps::ComputeDistance(pos_i, fluid_particles.position[j]);
    double weight = mps::WeightFunction(dist, smoothing_radius);
    
    // 计算加权中心位置
    weighted_center[0] += weight * fluid_particles.position[j][0];
    weighted_center[1] += weight * fluid_particles.position[j][1];
    weighted_center[2] += weight * fluid_particles.position[j][2];
    total_weight += weight;
  }

  // 计算固体粒子贡献
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    double dist = mps::ComputeDistance(pos_i, solid_particles.position[j]);
    double weight = mps::WeightFunction(dist, smoothing_radius);
    
    weighted_center[0] += weight * solid_particles.position[j][0];
    weighted_center[1] += weight * solid_particles.position[j][1];
    weighted_center[2] += weight * solid_particles.position[j][2];
    total_weight += weight;
  }

  // 计算偏置向量：从加权中心指向当前粒子（指向最稀疏的方向）
  if (total_weight > 1e-10) {
    weighted_center[0] /= total_weight;
    weighted_center[1] /= total_weight;
    weighted_center[2] /= total_weight;
    
    double3 bias_vec = {
      pos_i[0] - weighted_center[0],
      pos_i[1] - weighted_center[1],
      pos_i[2] - weighted_center[2]
    };
    
    return bias_vec;
  }
  
  // 如果没有邻域粒子，返回零向量
  return {0.0, 0.0, 0.0};
}

bool SurfaceDetector::IsSplashParticle(
    double number_density,
    double reference_density) const {
  
  // 如果粒子数密度小于参考密度的阈值比例，判定为飞溅粒子
  return number_density < (reference_density * SPLASH_DENSITY_RATIO);
}

bool SurfaceDetector::IsSurfaceParticle(
    int particle_idx,
    const double3& bias_vector,
    double bias_magnitude,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const double smoothing_radius,
    double bias_magnitude_threshold,
    bool use_virtual_light,
    double particle_radius) const {
  
  if (use_virtual_light) {
    // 使用虚拟光源法（不需要偏置向量模长检查，直接基于阴影面积判定）
    double shadow_ratio = ComputeShadowAreaRatio(
        particle_idx, fluid_particles, solid_particles, smoothing_radius, particle_radius);
    // 阴影面积比例小于阈值，判定为自由面粒子
    // 注意：内部粒子阴影面积大（接近1.0），自由面粒子阴影面积小（接近0.0）
    return shadow_ratio < SHADOW_AREA_THRESHOLD;
  } else {
    // 使用锥形区域法（需要偏置向量模长检查）
    // 偏置向量必须足够大
    if (bias_magnitude < bias_magnitude_threshold) {
      return false;
    }
    
    int cone_particle_count = CountParticlesInCone(
        particle_idx, bias_vector, fluid_particles, solid_particles, smoothing_radius);
    // 如果锥形区域内粒子数量很少，判定为自由面粒子
    return cone_particle_count <= CONE_PARTICLE_THRESHOLD;
  }
}

int SurfaceDetector::CountParticlesInCone(
    int particle_idx,
    const double3& bias_vector,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const double smoothing_radius) const {
  
  // smoothing_radius参数保留用于未来可能的扩展（如距离加权）
  (void)smoothing_radius;
  
  int count = 0;
  const double3& pos_i = fluid_particles.position[particle_idx];

  // 检查流体粒子邻域
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    double3 r_ij = {
      fluid_particles.position[j][0] - pos_i[0],
      fluid_particles.position[j][1] - pos_i[1],
      fluid_particles.position[j][2] - pos_i[2]
    };
    double dist = mps::ComputeVectorMagnitude(r_ij);
    
    if (dist > 1e-10) {
      double3 r_ij_normalized = mps::NormalizeVector(r_ij);
      double cos_angle = mps::DotProduct(bias_vector, r_ij_normalized);
      
      // 如果粒子在锥形区域内（cos值大于阈值）
      if (cos_angle >= CONE_COS_THRESHOLD) {
        ++count;
      }
    }
  }

  // 检查固体粒子邻域
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    double3 r_ij = {
      solid_particles.position[j][0] - pos_i[0],
      solid_particles.position[j][1] - pos_i[1],
      solid_particles.position[j][2] - pos_i[2]
    };
    double dist = mps::ComputeVectorMagnitude(r_ij);
    
    if (dist > 1e-10) {
      double3 r_ij_normalized = mps::NormalizeVector(r_ij);
      double cos_angle = mps::DotProduct(bias_vector, r_ij_normalized);
      
      if (cos_angle >= CONE_COS_THRESHOLD) {
        ++count;
      }
    }
  }

  return count;
}

bool SurfaceDetector::IsNearSurfaceParticle(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const double smoothing_radius) const {
  
  const double3& pos_i = fluid_particles.position[particle_idx];
  const double near_surface_distance = smoothing_radius * NEAR_SURFACE_DISTANCE_RATIO;

  // 检查邻域内是否有自由面粒子
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    if (fluid_particles.surface_type[j] == SurfaceType::SURFACE) {
      double dist = mps::ComputeDistance(pos_i, fluid_particles.position[j]);
      if (dist <= near_surface_distance) {
        return true;
      }
    }
  }

  return false;
}

double SurfaceDetector::ComputeReferenceDensity(
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const double smoothing_radius) const {
  
  // 计算所有粒子的平均粒子数密度作为参考
  // 这里使用所有粒子的平均密度，实际应用中可能需要更精确的方法
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

void mark_shadow_grid(std::vector<std::vector<bool>>& shadow_grid,
                      const double3& ri,
                      const double3& rj,
                      const double particle_radius,
                      const double smoothing_radius,
                      const int num_theta,
                      const int num_phi,
                      const double dtheta,
                      const double dphi) {
  // 计算粒子间距离及三方向分量长度
  const double3 r_ij = {
    rj[0] - ri[0],
    rj[1] - ri[1],
    rj[2] - ri[2]
  };
  const double rr = mps::ComputeVectorMagnitude(r_ij);
  
  // 检查距离是否有效
  if (rr < 1e-10) {
    return;
  }
  
  const double rx = r_ij[0];
  const double ry = r_ij[1];
  const double rz = r_ij[2];
  
  // 计算方位角theta及phi
  // theta: 从z轴正方向的夹角，范围[0, π]
  const double z_component = std::max(-1.0, std::min(1.0, rz / rr));
  const double theta = std::acos(z_component);
  
  // phi: 绕z轴的旋转角，范围[0, 2π]
  // 注意：原代码中 phi = atan2(ry, rx) + M_PI，这里保持原逻辑
  double phi = std::atan2(ry, rx) + M_PI;
  // 确保phi在[0, 2π]范围内
  if (phi >= 2.0 * M_PI) {
    phi -= 2.0 * M_PI;
  }
  if (phi < 0.0) {
    phi += 2.0 * M_PI;
  }
  
  // 计算theta及phi对应的网格索引
  int theta_idx_center = static_cast<int>(theta / dtheta);
  if (theta_idx_center < 0) theta_idx_center = 0;
  if (theta_idx_center >= num_theta) theta_idx_center = num_theta - 1;
  
  int phi_idx_center = static_cast<int>(phi / dphi);
  if (phi_idx_center < 0) phi_idx_center = 0;
  if (phi_idx_center >= num_phi) phi_idx_center = num_phi - 1;
  
  // 计算遮挡角度
  // 检查particle_radius / rr是否在有效范围内
  const double radius_ratio = particle_radius / rr;
  double occlusion_angle;
  if (radius_ratio >= 1.0) {
    // 如果粒子间距小于粒子半径，将目标粒子视为平行光源，遮挡角会比较小
    occlusion_angle = std::asin(particle_radius / smoothing_radius);
  } else {
    occlusion_angle = std::asin(radius_ratio);
  }
  
  // 计算遮挡范围对应的网格索引范围
  const int theta_range = static_cast<int>(std::ceil(occlusion_angle / dtheta));
  const int phi_range = static_cast<int>(std::ceil(occlusion_angle / dphi));
  
  // 标记被遮挡的网格
  for (int dt = -theta_range; dt <= theta_range; ++dt) {
    int theta_idx = theta_idx_center + dt;
    // theta索引越界，忽略越界部分
    if (theta_idx < 0 || theta_idx >= num_theta) continue;
    
    for (int dp = -phi_range; dp <= phi_range; ++dp) {
      int phi_idx = phi_idx_center + dp;
      // phi索引越界，负数索引转化为倒数索引，越界大索引转化为小索引
      if (phi_idx < 0) phi_idx += num_phi;
      if (phi_idx >= num_phi) phi_idx -= num_phi;
      shadow_grid[theta_idx][phi_idx] = true;
    }
  }
}

double SurfaceDetector::ComputeShadowAreaRatio(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const double smoothing_radius,
    const double particle_radius) const {
  
  const double3& pos_i = fluid_particles.position[particle_idx];
  
  // 步骤1：将球面幕布划分为网格
  // 使用等面积网格划分方式，确保对称性
  // 使用球面坐标系统：theta (0到π), phi (0到2π)
  const int num_theta = SPHERE_GRID_THETA_RESOLUTION;
  const int num_phi = SPHERE_GRID_PHI_RESOLUTION;
  const double dtheta = M_PI / num_theta;
  const double dphi = 2.0 * M_PI / num_phi;
  
  // 标记哪些网格被遮挡（使用二维数组）
  std::vector<std::vector<bool>> shadow_grid(num_theta, std::vector<bool>(num_phi, false));
  
  // 统计被遮挡的网格数量（用于验证）
  int shadowed_grid_count = 0;
  
  // 步骤2：遍历邻域粒子，标记被遮挡的网格
  // 处理流体粒子
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    double3 r_ij = {
      fluid_particles.position[j][0] - pos_i[0],
      fluid_particles.position[j][1] - pos_i[1],
      fluid_particles.position[j][2] - pos_i[2]
    };    
    const double dist = mps::ComputeVectorMagnitude(r_ij);
    if (dist > 1e-10 && dist <= smoothing_radius) {
      mark_shadow_grid(shadow_grid, pos_i, fluid_particles.position[j], particle_radius,
                      smoothing_radius, num_theta, num_phi, dtheta, dphi);
    }
  }
  // 处理固体粒子
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    double3 r_ij = {
      solid_particles.position[j][0] - pos_i[0],
      solid_particles.position[j][1] - pos_i[1],
      solid_particles.position[j][2] - pos_i[2]
    };
    double dist = mps::ComputeVectorMagnitude(r_ij);
    
    if (dist > 1e-10 && dist <= smoothing_radius) {
      mark_shadow_grid(shadow_grid, pos_i, solid_particles.position[j], particle_radius,
                      smoothing_radius, num_theta, num_phi, dtheta, dphi);
    }
  }
  
  // 步骤3：统计被遮挡网格的总面积
  double total_shadow_area = 0.0;
  double total_sphere_area = 0.0;
  
  for (int i = 0; i < num_theta; ++i) {
    double theta = (i + 0.5) * dtheta;
    // 球面网格面积 = sin(theta) * dtheta * dphi * R^2
    // 这里R = smoothing_radius，但比例计算时R^2会约掉
    double grid_area_weight = std::sin(theta) * dtheta * dphi;
    for (int j = 0; j < num_phi; ++j) {
      total_sphere_area += grid_area_weight;
      if (shadow_grid[i][j]) {
        total_shadow_area += grid_area_weight;
      }
    }
  }
  
  // 返回阴影面积比例
  if (total_sphere_area > 1e-10) {
    return total_shadow_area / total_sphere_area;
  } else {
    return 0.0;
  }
}

} // namespace mps
