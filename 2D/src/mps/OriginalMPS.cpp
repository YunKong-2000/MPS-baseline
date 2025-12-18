#include "OriginalMPS.hpp"
#include <cmath>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace mps2D {

double2 OriginalMPS::ComputeGradient(
    int particle_idx,
    const std::vector<double>& scalar_field,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double smoothing_radius,
    double reference_density) const {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  double phi_i = scalar_field[particle_idx];
  
  double grad_x = 0.0;
  double grad_y = 0.0;
  
  // 处理流体邻域粒子
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    
    // 避免除零和超出平滑半径的情况
    if (dist < 1e-10 || dist > smoothing_radius) {
      continue;
    }
    
    double phi_j = scalar_field[j];
    double weight = WeightFunction(dist, smoothing_radius);
    
    // 原始MPS梯度公式：∇φ = (d/n0) * Σ[(φj - φi) / |rij|^2 * rij * w(rij)]
    // 其中d=2（2D情况）
    double factor = DIMENSION * weight * (phi_j - phi_i) / (reference_density * dist * dist);
    grad_x += factor * dx;
    grad_y += factor * dy;
  }
  
  // 处理固体邻域粒子（边界条件）
  // 对于固体粒子，通常使用边界值（如壁面处速度为零）
  // 这里假设固体粒子处的标量值为0（可以根据实际边界条件修改）
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    const double2& pos_j = solid_particles.position[j];
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    
    if (dist < 1e-10 || dist > smoothing_radius) {
      continue;
    }
    
    // 假设固体边界处的标量值为0（Dirichlet边界条件）
    // 可以根据实际需求修改边界值
    double phi_boundary = 0.0;
    double weight = WeightFunction(dist, smoothing_radius);
    
    double factor = DIMENSION * weight * (phi_boundary - phi_i) / (reference_density * dist * dist);
    grad_x += factor * dx;
    grad_y += factor * dy;
  }
  
  return {grad_x, grad_y};
}

double OriginalMPS::ComputeDivergence(
    int particle_idx,
    const std::vector<double2>& vector_field,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double smoothing_radius,
    double reference_density) const {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  const double2& u_i = vector_field[particle_idx];
  
  double divergence = 0.0;
  
  // 处理流体邻域粒子
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    
    if (dist < 1e-10 || dist > smoothing_radius) {
      continue;
    }
    
    const double2& u_j = vector_field[j];
    double2 du = {u_j.x - u_i.x, u_j.y - u_i.y};
    double2 r_ij = {dx, dy};
    
    double weight = WeightFunction(dist, smoothing_radius);
    
    // 原始MPS散度公式：∇·u = (d/n0) * Σ[(uj - ui) · rij / |rij|^2 * w(rij)]
    // 其中d=2（2D情况）
    double dot_product = du.x * r_ij.x + du.y * r_ij.y;
    divergence += DIMENSION * weight * dot_product / (reference_density * dist * dist);
  }
  
  // 处理固体邻域粒子（边界条件）
  // 假设壁面处速度为零（无滑移边界条件）
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    const double2& pos_j = solid_particles.position[j];
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    
    if (dist < 1e-10 || dist > smoothing_radius) {
      continue;
    }
    
    // 壁面处速度为零
    double2 u_boundary = {0.0, 0.0};
    double2 du = {u_boundary.x - u_i.x, u_boundary.y - u_i.y};
    double2 r_ij = {dx, dy};
    
    double weight = WeightFunction(dist, smoothing_radius);
    
    double dot_product = du.x * r_ij.x + du.y * r_ij.y;
    divergence += DIMENSION * weight * dot_product / (reference_density * dist * dist);
  }
  
  return divergence;
}

double OriginalMPS::ComputeLaplacian(
    int particle_idx,
    const std::vector<double>& scalar_field,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double smoothing_radius,
    double reference_density,
    double particle_spacing) const {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  double phi_i = scalar_field[particle_idx];
  
  // 计算lambda参数（对所有粒子都相同，使用离散方法计算）
  double lambda = ComputeLambdaDiscreteUniform(particle_spacing, smoothing_radius);
  
  double laplacian = 0.0;
  
  // 处理流体邻域粒子
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];
    double dist = ComputeDistance(pos_i, pos_j);
    
    if (dist < 1e-10 || dist > smoothing_radius) {
      continue;
    }
    
    double phi_j = scalar_field[j];
    double weight = WeightFunction(dist, smoothing_radius);
    
    // 原始MPS拉普拉斯公式：Δφ = (2d/(n0*λ)) * Σ[(φj - φi) * w(rij)]
    // 其中d=2（2D情况）
    laplacian += (2.0 * DIMENSION * weight * (phi_j - phi_i)) / (reference_density * lambda);
  }
  
  // 处理固体邻域粒子（边界条件）
  // 假设固体边界处的标量值为0
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    const double2& pos_j = solid_particles.position[j];
    double dist = ComputeDistance(pos_i, pos_j);
    
    if (dist < 1e-10 || dist > smoothing_radius) {
      continue;
    }
    
    double phi_boundary = 0.0;
    double weight = WeightFunction(dist, smoothing_radius);
    
    laplacian += (2.0 * DIMENSION * weight * (phi_boundary - phi_i)) / (reference_density * lambda);
  }
  
  return laplacian;
}

double OriginalMPS::ComputeReferenceDensity(
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double smoothing_radius) const {
  
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

double OriginalMPS::ComputeLambda(double smoothing_radius) const {
  // lambda参数直接从对作用域的定积分计算获得
  // λ = (1/5) * r_e^2
  // 这是从积分公式推导出的解析解，对所有粒子都相同
  return (1.0 / 5.0) * smoothing_radius * smoothing_radius;
}

double OriginalMPS::ComputeLambdaDiscrete(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double smoothing_radius) const {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  
  double sum_weight_dist2 = 0.0;
  double sum_weight = 0.0;
  
  // 处理流体邻域粒子
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];
    double dist = ComputeDistance(pos_i, pos_j);
    
    if (dist < 1e-10 || dist > smoothing_radius) {
      continue;
    }
    
    double weight = WeightFunction(dist, smoothing_radius);
    sum_weight_dist2 += weight * dist * dist;
    sum_weight += weight;
  }
  
  // 处理固体邻域粒子
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    const double2& pos_j = solid_particles.position[j];
    double dist = ComputeDistance(pos_i, pos_j);
    
    if (dist < 1e-10 || dist > smoothing_radius) {
      continue;
    }
    
    double weight = WeightFunction(dist, smoothing_radius);
    sum_weight_dist2 += weight * dist * dist;
    sum_weight += weight;
  }
  
  // lambda = Σ(w_ij * |rij|^2) / Σ(w_ij)
  if (sum_weight > 1e-10) {
    return sum_weight_dist2 / sum_weight;
  } else {
    // 如果没有邻域粒子，返回解析解作为默认值
    return ComputeLambda(smoothing_radius);
  }
}

double OriginalMPS::ComputeLambdaDiscreteUniform(
    double particle_spacing,
    double smoothing_radius) const {
  
  // 直接生成均匀分布的粒子配置来计算离散lambda
  return ComputeLambdaFromUniformDistribution(
      particle_spacing, smoothing_radius, {0.0, 0.0});
}

double OriginalMPS::ComputeLambdaFromUniformDistribution(
    double particle_spacing,
    double smoothing_radius,
    const double2& center_pos) const {
  
  // 生成均匀分布的粒子配置
  // 在平滑半径范围内生成规则网格的粒子
  double sum_weight_dist2 = 0.0;
  double sum_weight = 0.0;
  int num_neighbors = 0;
  
  // 计算需要覆盖的范围（以平滑半径为半径的圆形区域）
  // 为了确保圆形区域内的所有粒子都被覆盖，需要考虑对角线方向
  // 对角线距离 = grid_size * particle_spacing * sqrt(2)
  // 要确保对角线距离 >= smoothing_radius，所以：
  // grid_size >= smoothing_radius / (particle_spacing * sqrt(2))
  // 增加额外的边界层（2-3个粒子间距）确保完全覆盖
  const double sqrt2 = 1.4142135623730951;
  int grid_size = static_cast<int>(std::ceil(smoothing_radius / (particle_spacing * sqrt2))) + 3;
  
  // 生成均匀分布的粒子网格
  for (int i = -grid_size; i <= grid_size; ++i) {
    for (int j = -grid_size; j <= grid_size; ++j) {
      // 计算粒子位置（相对于中心）
      double2 pos = {
        center_pos.x + i * particle_spacing,
        center_pos.y + j * particle_spacing
      };
      
      // 计算距离
      double dx = pos.x - center_pos.x;
      double dy = pos.y - center_pos.y;
      double dist = std::sqrt(dx * dx + dy * dy);
      
      // 跳过中心粒子本身
      if (dist < 1e-10) {
        continue;
      }
      
      // 只考虑在平滑半径内的粒子
      if (dist > smoothing_radius) {
        continue;
      }
      
      // 计算权重
      double weight = WeightFunction(dist, smoothing_radius);
      
      // 累加
      sum_weight_dist2 += weight * dist * dist;
      sum_weight += weight;
      ++num_neighbors;
    }
  }
  
  // lambda = Σ(w_ij * |rij|^2) / Σ(w_ij)
  if (sum_weight > 1e-10) {
    return sum_weight_dist2 / sum_weight;
  } else {
    // 如果没有邻域粒子，返回解析解作为默认值
    return ComputeLambda(smoothing_radius);
  }
}

OriginalMPS::LambdaComparisonResult OriginalMPS::CompareLambdaMethods(
    double particle_spacing,
    double smoothing_radius) const {
  
  LambdaComparisonResult result;
  
  // 计算解析lambda
  result.analytical_lambda = ComputeLambda(smoothing_radius);
  
  // 使用均匀分布计算离散lambda，同时统计邻域粒子数
  double2 center_pos = {0.0, 0.0};
  double sum_weight_dist2 = 0.0;
  double sum_weight = 0.0;
  int num_neighbors = 0;
  
  // 计算需要覆盖的范围（以平滑半径为半径的圆形区域）
  // 为了确保圆形区域内的所有粒子都被覆盖，需要考虑对角线方向
  // 对角线距离 = grid_size * particle_spacing * sqrt(2)
  // 要确保对角线距离 >= smoothing_radius，所以：
  // grid_size >= smoothing_radius / (particle_spacing * sqrt(2))
  // 增加额外的边界层（2-3个粒子间距）确保完全覆盖
  const double sqrt2 = 1.4142135623730951;
  int grid_size = static_cast<int>(std::ceil(smoothing_radius / (particle_spacing * sqrt2))) + 3;
  
  // 生成均匀分布的粒子网格
  for (int i = -grid_size; i <= grid_size; ++i) {
    for (int j = -grid_size; j <= grid_size; ++j) {
      // 计算粒子位置（相对于中心）
      double2 pos = {
        center_pos.x + i * particle_spacing,
        center_pos.y + j * particle_spacing
      };
      
      // 计算距离
      double dx = pos.x - center_pos.x;
      double dy = pos.y - center_pos.y;
      double dist = std::sqrt(dx * dx + dy * dy);
      
      // 跳过中心粒子本身
      if (dist < 1e-10) {
        continue;
      }
      
      // 只考虑在平滑半径内的粒子
      if (dist > smoothing_radius) {
        continue;
      }
      
      // 计算权重
      double weight = WeightFunction(dist, smoothing_radius);
      
      // 累加
      sum_weight_dist2 += weight * dist * dist;
      sum_weight += weight;
      ++num_neighbors;
    }
  }
  
  // 计算离散lambda
  if (sum_weight > 1e-10) {
    result.discrete_lambda = sum_weight_dist2 / sum_weight;
  } else {
    result.discrete_lambda = result.analytical_lambda;
  }
  
  result.num_neighbors = num_neighbors;
  
  // 计算差值
  result.difference = std::abs(result.discrete_lambda - result.analytical_lambda);
  
  // 计算相对差值（百分比）
  if (result.analytical_lambda > 1e-10) {
    result.relative_difference = (result.difference / result.analytical_lambda) * 100.0;
  } else {
    result.relative_difference = 0.0;
  }
  
  return result;
}

double OriginalMPS::ComputeParticleNumberDensity(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double smoothing_radius) const {
  
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

} // namespace mps2D
