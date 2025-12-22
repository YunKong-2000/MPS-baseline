#include "../src/mps/OriginalMPS.hpp"
#include "../src/core/Particle.hpp"
#include "core/Types.h"
#include "core/MPSUtils.h"
#include "../src/neighbour_list/NeighborListSearcher.hpp"
#include "../src/surface_detection/SurfaceDetector.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <limits>

using namespace mps2D;

// 生成静水压强测试场景
// 参数：
//   container_width: 容器宽度
//   container_height: 容器高度
//   water_height: 水位高度
//   particle_spacing: 粒子间距
//   fluid_particles: 输出流体粒子
//   solid_particles: 输出固体粒子
void GenerateHydrostaticTest(
    double container_width,
    double container_height,
    double water_height,
    double particle_spacing,
    FluidParticle& fluid_particles,
    SolidParticle& solid_particles) {
  
  // 计算粒子数量
  int nx_fluid = static_cast<int>(container_width / particle_spacing) + 1;
  int ny_fluid = static_cast<int>(water_height / particle_spacing) + 1;
  int num_fluid = nx_fluid * ny_fluid;
  
  // 生成流体粒子（水）
  fluid_particles.particle_num = num_fluid;
  fluid_particles.position.resize(num_fluid);
  fluid_particles.velocity.resize(num_fluid);
  fluid_particles.pressure.resize(num_fluid);
  fluid_particles.density.resize(num_fluid);
  fluid_particles.surface_type.resize(num_fluid);
  fluid_particles.fluid_neighbour_list.resize(num_fluid);
  fluid_particles.solid_neighbour_list.resize(num_fluid);
  
  int idx = 0;
  for (int j = 0; j < ny_fluid; ++j) {
    for (int i = 0; i < nx_fluid; ++i) {
      double x = i * particle_spacing;
      double y = j * particle_spacing;
      fluid_particles.position[idx] = {x, y};
      fluid_particles.velocity[idx] = {0.0, 0.0};
      fluid_particles.pressure[idx] = 0.0;  // 稍后计算
      fluid_particles.density[idx] = 1000.0;  // 水的密度
      fluid_particles.surface_type[idx] = SurfaceType::INNER;
      ++idx;
    }
  }
  
  // 生成固体粒子（容器壁面）
  // 底部壁面
  int nx_bottom = static_cast<int>(container_width / particle_spacing) + 1;
  // 左侧壁面
  int ny_left = static_cast<int>(container_height / particle_spacing) + 1;
  // 右侧壁面
  int ny_right = ny_left;
  
  int num_solid = nx_bottom + 2 * ny_left;
  solid_particles.particle_num = num_solid;
  solid_particles.position.resize(num_solid);
  solid_particles.velocity.resize(num_solid);
  solid_particles.normal_vector.resize(num_solid);
  
  idx = 0;
  
  // 底部壁面（法向量向上）
  for (int i = 0; i < nx_bottom; ++i) {
    double x = i * particle_spacing;
    double y = -particle_spacing;  // 距离底部流体粒子一个粒子间距
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {0.0, 1.0};  // 向上
    ++idx;
  }
  
  // 左侧壁面（法向量向右）
  for (int j = 0; j < ny_left; ++j) {
    double x = -particle_spacing;  // 距离左侧流体粒子一个粒子间距
    double y = j * particle_spacing;
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {1.0, 0.0};  // 向右
    ++idx;
  }
  
  // 右侧壁面（法向量向左）
  for (int j = 0; j < ny_right; ++j) {
    double x = container_width;  // 距离右侧流体粒子一个粒子间距
    double y = j * particle_spacing;
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {-1.0, 0.0};  // 向左
    ++idx;
  }
}

// 计算理论压力（静水压强）
// 参数：
//   rho: 流体密度
//   g: 重力加速度
//   h: 水深（从底部到粒子位置）
// 返回：压力值
double ComputeTheoreticalPressure(double rho, double g, double h) {
  return rho * g * h;
}

// 权重函数（原始MPS）- 使用命名空间限定的全局函数以避免冲突
// 注意：由于MPSUtils.h中已经定义了WeightFunction，这里直接使用mps2D::WeightFunction

// 计算自由面粒子的拉普拉斯算子（根据ppeMatrix.md文档）
// 参数：
//   particle_idx: 粒子索引
//   scalar_field: 标量场（压力）
//   fluid_particles: 流体粒子对象
//   solid_particles: 固体粒子对象
//   smoothing_radius: 平滑半径
//   reference_density: 参考粒子数密度 n_0
//   lambda: lambda参数
// 返回：拉普拉斯算子值
double ComputeLaplacianFreeSurface(
    int particle_idx,
    const std::vector<double>& scalar_field,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double smoothing_radius,
    double reference_density,
    double lambda) {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  double phi_i = scalar_field[particle_idx];
  
  // 计算粒子数密度 n_i^*
  double n_i_star = 0.0;
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    double dist = std::sqrt(
        std::pow(fluid_particles.position[j].x - pos_i.x, 2) +
        std::pow(fluid_particles.position[j].y - pos_i.y, 2));
    if (dist < smoothing_radius) {
      n_i_star += mps2D::WeightFunction(dist, smoothing_radius);
    }
  }
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    double dist = std::sqrt(
        std::pow(solid_particles.position[j].x - pos_i.x, 2) +
        std::pow(solid_particles.position[j].y - pos_i.y, 2));
    if (dist < smoothing_radius) {
      n_i_star += mps2D::WeightFunction(dist, smoothing_radius);
    }
  }
  
  // 计算 w_{l_0}（参考距离处的权重，通常取 r_e 处的权重值，即0）
  double w_l0 = mps2D::WeightFunction(smoothing_radius, smoothing_radius);  // 应该为0
  
  // 计算 \tilde{n}_i = n_0 + 0.1 * sum_{j≠i}^{N_r} (w_{ij} - w_{l_0})
  double sum_weight_diff = 0.0;
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    double dist = std::sqrt(
        std::pow(fluid_particles.position[j].x - pos_i.x, 2) +
        std::pow(fluid_particles.position[j].y - pos_i.y, 2));
    if (dist < smoothing_radius) {
      double w_ij = mps2D::WeightFunction(dist, smoothing_radius);
      sum_weight_diff += (w_ij - w_l0);
    }
  }
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    double dist = std::sqrt(
        std::pow(solid_particles.position[j].x - pos_i.x, 2) +
        std::pow(solid_particles.position[j].y - pos_i.y, 2));
    if (dist < smoothing_radius) {
      double w_ij = mps2D::WeightFunction(dist, smoothing_radius);
      sum_weight_diff += (w_ij - w_l0);
    }
  }
  double n_i_tilde = reference_density + 0.1 * sum_weight_diff;
  
  // 计算 n_i' = max(n_i^*, \tilde{n}_i)
  double n_i_prime = std::max(n_i_star, n_i_tilde);
  
  // 计算拉普拉斯算子
  // ⟨Δp⟩_i = (4/(n_0*λ)) * sum_{j∈fluid} (p_j - p_i) * w_{ij} 
  //          - (4/(n_0*λ)) * (n_i' - n_i^*) * p_i
  double laplacian = 0.0;
  double coeff = 4.0 / (reference_density * lambda);
  
  // 第一项：sum_{j∈fluid} (p_j - p_i) * w_{ij}
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    double dist = std::sqrt(
        std::pow(fluid_particles.position[j].x - pos_i.x, 2) +
        std::pow(fluid_particles.position[j].y - pos_i.y, 2));
    if (dist < smoothing_radius) {
      double w_ij = mps2D::WeightFunction(dist, smoothing_radius);
      double phi_j = scalar_field[j];
      laplacian += coeff * (phi_j - phi_i) * w_ij;
    }
  }
  
  // 第二项：- (4/(n_0*λ)) * (n_i' - n_i^*) * p_i
  laplacian -= coeff * (n_i_prime - n_i_star) * phi_i;
  
  return laplacian;
}

// 计算自由面粒子的速度散度（根据ppeMatrix.md文档）
// 参数：
//   particle_idx: 粒子索引
//   vector_field: 向量场（速度）
//   fluid_particles: 流体粒子对象
//   solid_particles: 固体粒子对象
//   smoothing_radius: 平滑半径
//   reference_density: 参考粒子数密度 n_0
// 返回：散度值
double ComputeDivergenceFreeSurface(
    int particle_idx,
    const std::vector<double2>& vector_field,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double smoothing_radius,
    double reference_density) {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  const double2& u_i = vector_field[particle_idx];
  
  // ⟨∇·u⟩_i = (2/n_0) * sum_{j≠i} [(u_j - u_i)/r_{ij} · r_{ij}/r_{ij} * w_{ij}]
  double divergence = 0.0;
  double coeff = 2.0 / reference_density;
  
  // 处理流体邻域粒子
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = std::sqrt(dx * dx + dy * dy);
    
    if (dist < 1e-10 || dist > smoothing_radius) {
      continue;
    }
    
    const double2& u_j = vector_field[j];
    double2 du = {u_j.x - u_i.x, u_j.y - u_i.y};
    double2 r_ij = {dx, dy};
    double r_ij_norm = dist;
    
    double w_ij = mps2D::WeightFunction(dist, smoothing_radius);
    
    // (u_j - u_i) / r_{ij} · r_{ij} / r_{ij}
    double dot_product = (du.x * r_ij.x + du.y * r_ij.y) / (r_ij_norm * r_ij_norm);
    divergence += coeff * dot_product * w_ij;
  }
  
  // 处理固体邻域粒子（壁面处速度为零）
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    const double2& pos_j = solid_particles.position[j];
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = std::sqrt(dx * dx + dy * dy);
    
    if (dist < 1e-10 || dist > smoothing_radius) {
      continue;
    }
    
    double2 u_boundary = {0.0, 0.0};
    double2 du = {u_boundary.x - u_i.x, u_boundary.y - u_i.y};
    double2 r_ij = {dx, dy};
    double r_ij_norm = dist;
    
    double w_ij = mps2D::WeightFunction(dist, smoothing_radius);
    
    double dot_product = (du.x * r_ij.x + du.y * r_ij.y) / (r_ij_norm * r_ij_norm);
    divergence += coeff * dot_product * w_ij;
  }
  
  return divergence;
}

// 输出压力数据到VTK文件（用于调试）
void WritePressureToVTK(
    const std::string& filename,
    const FluidParticle& fluid_particles,
    const std::vector<double>& theoretical_pressure,
    const std::vector<double2>& computed_gradient,
    const std::vector<double2>& theoretical_gradient,
    const std::vector<double>& computed_laplacian) {
  
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "错误：无法打开文件 " << filename << std::endl;
    return;
  }
  
  int num_particles = fluid_particles.particle_num;
  
  file << std::fixed << std::setprecision(15);
  
  // VTK文件头
  file << "# vtk DataFile Version 3.0\n";
  file << "Hydrostatic Pressure Test - Original MPS 2D\n";
  file << "ASCII\n";
  file << "DATASET POLYDATA\n";
  
  // 写入点坐标（2D数据添加z=0）
  file << "POINTS " << num_particles << " float\n";
  for (int i = 0; i < num_particles; ++i) {
    const auto& pos = fluid_particles.position[i];
    file << pos.x << " " << pos.y << " 0.0\n";
  }
  
  // 写入顶点（每个点作为一个顶点）
  file << "VERTICES " << num_particles << " " << (num_particles * 2) << "\n";
  for (int i = 0; i < num_particles; ++i) {
    file << "1 " << i << "\n";
  }
  
  // 写入点数据
  file << "POINT_DATA " << num_particles << "\n";
  
  // 写入理论压力
  file << "SCALARS theoretical_pressure float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    file << theoretical_pressure[i] << "\n";
  }
  
  // 写入计算出的压力梯度
  file << "VECTORS computed_pressure_gradient float\n";
  for (int i = 0; i < num_particles; ++i) {
    file << computed_gradient[i].x << " " 
         << computed_gradient[i].y << " 0.0\n";
  }
  
  // 写入理论梯度（压力梯度）
  file << "VECTORS theoretical_pressure_gradient float\n";
  for (int i = 0; i < num_particles; ++i) {
    file << theoretical_gradient[i].x << " " 
         << theoretical_gradient[i].y << " 0.0\n";
  }
  
  // 写入梯度误差
  file << "VECTORS gradient_error float\n";
  for (int i = 0; i < num_particles; ++i) {
    double error_x = computed_gradient[i].x - theoretical_gradient[i].x;
    double error_y = computed_gradient[i].y - theoretical_gradient[i].y;
    file << error_x << " " << error_y << " 0.0\n";
  }
  
  // 写入梯度误差大小
  file << "SCALARS gradient_error_magnitude float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    double error_x = computed_gradient[i].x - theoretical_gradient[i].x;
    double error_y = computed_gradient[i].y - theoretical_gradient[i].y;
    double error_mag = std::sqrt(error_x * error_x + error_y * error_y);
    file << error_mag << "\n";
  }
  
  // 写入压力拉普拉斯算子
  file << "SCALARS pressure_laplacian float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    file << computed_laplacian[i] << "\n";
  }
  
  // 写入Y坐标（用于检查压力分布）
  file << "SCALARS y_coordinate float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    file << fluid_particles.position[i].y << "\n";
  }
  
  // 写入表面类型（用于可视化自由面粒子）
  file << "SCALARS surface_type int\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    int type_value = static_cast<int>(fluid_particles.surface_type[i]);
    file << type_value << "\n";
  }
  
  file.close();
  std::cout << "  已输出VTK文件: " << filename << std::endl;
}

int main() {
  std::cout << "=== 静水压强测例（原始MPS方法） ===" << std::endl;
  
  // 物理参数
  const double rho = 1000.0;  // 水的密度 (kg/m³)
  const double g = 9.8;       // 重力加速度 (m/s²)
  
  // 几何参数
  const double container_width = 1.0;   // 容器宽度 (m)
  const double container_height = 1.0;  // 容器高度 (m)
  const double water_height = 0.5;      // 水位高度 (m)
  
  // 粒子参数
  const double particle_spacing = 0.007;  // 粒子间距 (m)
  const double particle_radius = particle_spacing / 2.0;
  const double smoothing_radius = 2.1 * particle_spacing;
  const double cell_size = 2.0 * smoothing_radius;
  
  std::cout << "\n物理参数:" << std::endl;
  std::cout << "  流体密度: " << rho << " kg/m³" << std::endl;
  std::cout << "  重力加速度: " << g << " m/s²" << std::endl;
  
  std::cout << "\n几何参数:" << std::endl;
  std::cout << "  容器宽度: " << container_width << " m" << std::endl;
  std::cout << "  容器高度: " << container_height << " m" << std::endl;
  std::cout << "  水位高度: " << water_height << " m" << std::endl;
  
  std::cout << "\n粒子参数:" << std::endl;
  std::cout << "  粒子间距: " << particle_spacing << " m" << std::endl;
  std::cout << "  平滑半径: " << smoothing_radius << " m" << std::endl;
  
  // 生成测试场景
  FluidParticle fluid_particles("fluid");
  SolidParticle solid_particles("solid");
  
  std::cout << "\n生成测试场景..." << std::endl;
  GenerateHydrostaticTest(
      container_width, container_height, water_height,
      particle_spacing, fluid_particles, solid_particles);
  
  std::cout << "  流体粒子数: " << fluid_particles.particle_num << std::endl;
  std::cout << "  固体粒子数: " << solid_particles.particle_num << std::endl;
  std::cout << "  总粒子数: " << (fluid_particles.particle_num + solid_particles.particle_num) << std::endl;
  
  // 构建邻居列表
  std::cout << "\n构建邻居列表..." << std::endl;
  NeighborListSearcher neighbor_searcher;
  neighbor_searcher.BuildNeighborList(
      fluid_particles, solid_particles,
      particle_radius, smoothing_radius, cell_size);
  
  // 计算理论压力
  std::cout << "\n计算理论压力..." << std::endl;
  std::vector<double> theoretical_pressure(fluid_particles.particle_num);
  double max_pressure = 0.0;
  double min_pressure = std::numeric_limits<double>::max();
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    // 水深 = 自由面高度 - 粒子y坐标（粒子到自由面的距离）
    double h = water_height - fluid_particles.position[i].y;
    // 确保水深不为负（防止浮点误差）
    if (h < 0.0) h = 0.0;
    theoretical_pressure[i] = ComputeTheoreticalPressure(rho, g, h);
    fluid_particles.pressure[i] = theoretical_pressure[i];
    if (theoretical_pressure[i] > max_pressure) {
      max_pressure = theoretical_pressure[i];
    }
    if (theoretical_pressure[i] < min_pressure) {
      min_pressure = theoretical_pressure[i];
    }
  }
  std::cout << "  最大压力（底部）: " << max_pressure << " Pa" << std::endl;
  std::cout << "  最小压力（表面）: " << min_pressure << " Pa" << std::endl;
  
  // 检测自由面粒子
  std::cout << "\n检测自由面粒子..." << std::endl;
  SurfaceDetector surface_detector;
  surface_detector.DetectSurfaceParticles(
      fluid_particles, solid_particles, smoothing_radius, particle_spacing);
  
  // 统计自由面粒子数量
  int num_surface = 0;
  int num_inner = 0;
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    if (fluid_particles.surface_type[i] == SurfaceType::SURFACE) {
      num_surface++;
    } else if (fluid_particles.surface_type[i] == SurfaceType::INNER) {
      num_inner++;
    }
  }
  std::cout << "  自由面粒子数: " << num_surface << std::endl;
  std::cout << "  内部粒子数: " << num_inner << std::endl;
  
  // 计算参考粒子数密度
  std::cout << "\n计算参考粒子数密度..." << std::endl;
  OriginalMPS mps_calculator;
  double reference_density = mps_calculator.ComputeReferenceDensity(
      fluid_particles, solid_particles, smoothing_radius);
  std::cout << "  参考粒子数密度: " << reference_density << std::endl;
  
  // 比较两种lambda计算方法
  std::cout << "\n比较lambda计算方法（解析公式 vs 离散分布）..." << std::endl;
  auto lambda_comparison = mps_calculator.CompareLambdaMethods(
      particle_spacing, smoothing_radius);
  std::cout << "  解析lambda (λ = (1/5) * r_e^2，对所有粒子相同): " 
            << std::fixed << std::setprecision(10) 
            << lambda_comparison.analytical_lambda << std::endl;
  std::cout << "  离散lambda (使用均匀分布，粒子间距: " << particle_spacing 
            << "，估算邻域粒子数: " << lambda_comparison.num_neighbors 
            << "，对所有粒子相同): " 
            << lambda_comparison.discrete_lambda << std::endl;
  std::cout << "  绝对差值: " << lambda_comparison.difference << std::endl;
  std::cout << "  相对差值: " << std::setprecision(4) 
            << lambda_comparison.relative_difference << "%" << std::endl;
  
  // 计算压力梯度
  std::cout << "\n计算压力梯度（原始MPS方法）..." << std::endl;
  std::vector<double2> computed_gradient(fluid_particles.particle_num);
  std::vector<double2> theoretical_gradient(fluid_particles.particle_num);
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    // 使用原始MPS方法计算梯度
    computed_gradient[i] = mps_calculator.ComputeGradient(
        i, theoretical_pressure, fluid_particles, solid_particles,
        smoothing_radius, reference_density);
    
    // 理论梯度（静水压强）
    // 压力 p = rho * g * h，其中 h = water_height - y
    // dp/dy = rho * g * d(water_height - y)/dy = rho * g * (-1) = -rho * g
    // 所以y方向梯度应该是负数（压力随y增加而减小）
    theoretical_gradient[i] = {0.0, -rho * g};
  }
  
  // 计算压力拉普拉斯算子
  std::cout << "\n计算压力拉普拉斯算子..." << std::endl;
  std::cout << "  （内部粒子：原始MPS方法；自由面粒子：特殊公式）" << std::endl;
  std::vector<double> computed_laplacian(fluid_particles.particle_num);
  
  // 计算lambda参数（用于自由面粒子）
  double lambda = mps_calculator.CompareLambdaMethods(particle_spacing, smoothing_radius).discrete_lambda;
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    if (fluid_particles.surface_type[i] == SurfaceType::SURFACE) {
      // 自由面粒子：使用特殊公式
      computed_laplacian[i] = ComputeLaplacianFreeSurface(
          i, theoretical_pressure, fluid_particles, solid_particles,
          smoothing_radius, reference_density, lambda);
    } else {
      // 内部粒子：使用原始MPS方法计算拉普拉斯算子（使用离散lambda）
      computed_laplacian[i] = mps_calculator.ComputeLaplacian(
          i, theoretical_pressure, fluid_particles, solid_particles,
          smoothing_radius, reference_density, particle_spacing);
    }
  }
  
  // 统计自由面粒子和内部粒子的拉普拉斯算子
  std::cout << "  自由面粒子拉普拉斯算子统计:" << std::endl;
  double surface_laplacian_sum = 0.0;
  double inner_laplacian_sum = 0.0;
  int surface_count = 0;
  int inner_count = 0;
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    if (fluid_particles.surface_type[i] == SurfaceType::SURFACE) {
      surface_laplacian_sum += computed_laplacian[i];
      surface_count++;
    } else if (fluid_particles.surface_type[i] == SurfaceType::INNER) {
      inner_laplacian_sum += computed_laplacian[i];
      inner_count++;
    }
  }
  if (surface_count > 0) {
    std::cout << "    平均拉普拉斯算子: " << surface_laplacian_sum / surface_count << std::endl;
  }
  if (inner_count > 0) {
    std::cout << "  内部粒子拉普拉斯算子统计:" << std::endl;
    std::cout << "    平均拉普拉斯算子: " << inner_laplacian_sum / inner_count << std::endl;
  }
  
  // 统计梯度误差（相对误差）
  std::cout << "\n梯度误差统计（相对误差）:" << std::endl;
  double total_relative_error_x = 0.0;
  double total_relative_error_y = 0.0;
  double max_relative_error_x = 0.0;
  double max_relative_error_y = 0.0;
  int num_valid_gradient = fluid_particles.particle_num;
  int num_valid_x = 0;
  int num_valid_y = 0;
  
  const double theoretical_grad_x = 0.0;  // 理论X方向梯度
  const double theoretical_grad_y = -rho * g;  // 理论Y方向梯度
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    // X方向相对误差（理论值为0，使用绝对误差）
    double abs_error_x = std::abs(computed_gradient[i].x - theoretical_gradient[i].x);
    double relative_error_x = (std::abs(theoretical_grad_x) > 1e-10) ? 
        (abs_error_x / std::abs(theoretical_grad_x) * 100.0) : abs_error_x;
    total_relative_error_x += relative_error_x;
    if (relative_error_x > max_relative_error_x) max_relative_error_x = relative_error_x;
    ++num_valid_x;
    
    // Y方向相对误差
    double abs_error_y = std::abs(computed_gradient[i].y - theoretical_gradient[i].y);
    double relative_error_y = (std::abs(theoretical_grad_y) > 1e-10) ? 
        (abs_error_y / std::abs(theoretical_grad_y) * 100.0) : abs_error_y;
    total_relative_error_y += relative_error_y;
    if (relative_error_y > max_relative_error_y) max_relative_error_y = relative_error_y;
    ++num_valid_y;
  }
  
  double avg_relative_error_x = total_relative_error_x / num_valid_x;
  double avg_relative_error_y = total_relative_error_y / num_valid_y;
  
  std::cout << "  统计粒子数: " << num_valid_gradient << std::endl;
  std::cout << "  X方向平均相对误差: " << std::fixed << std::setprecision(4) 
            << avg_relative_error_x << "%" << std::endl;
  std::cout << "  Y方向平均相对误差: " << avg_relative_error_y << "%" << std::endl;
  std::cout << "  X方向最大相对误差: " << max_relative_error_x << "%" << std::endl;
  std::cout << "  Y方向最大相对误差: " << max_relative_error_y << "%" << std::endl;
  std::cout << "  理论X方向梯度: " << theoretical_grad_x << " Pa/m" << std::endl;
  std::cout << "  理论Y方向梯度: " << theoretical_grad_y << " Pa/m (压力随y增加而减小)" << std::endl;
  
  // 显示前几个粒子的详细信息（用于调试）
  std::cout << "\n前5个粒子的详细信息（用于调试）:" << std::endl;
  std::cout << std::setw(6) << "粒子" 
            << std::setw(12) << "位置Y" 
            << std::setw(12) << "压力"
            << std::setw(12) << "理论GradY"
            << std::setw(12) << "计算GradY"
            << std::setw(12) << "相对误差Y(%)"
            << std::setw(8) << "邻域数" << std::endl;
  
  int count = 0;
  for (int i = 0; i < fluid_particles.particle_num && count < 5; ++i) {
    double error_y = computed_gradient[i].y - theoretical_gradient[i].y;
    int num_neighbors = fluid_particles.fluid_neighbour_list[i].size() + 
                       fluid_particles.solid_neighbour_list[i].size();
    
      double relative_error_y = (std::abs(theoretical_gradient[i].y) > 1e-10) ?
          (std::abs(error_y) / std::abs(theoretical_gradient[i].y) * 100.0) : std::abs(error_y);
      
      std::cout << std::setw(6) << i
              << std::setw(12) << std::fixed << std::setprecision(4) 
              << fluid_particles.position[i].y
              << std::setw(12) << theoretical_pressure[i]
              << std::setw(12) << theoretical_gradient[i].y
              << std::setw(12) << computed_gradient[i].y
              << std::setw(12) << relative_error_y << "%"
              << std::setw(8) << num_neighbors << std::endl;
    ++count;
  }
  
  // 输出压力数据到VTK文件（用于调试）
  std::cout << "\n输出VTK文件..." << std::endl;
  WritePressureToVTK(
      "original_mps_hydrostatic_pressure.vtk",
      fluid_particles,
      theoretical_pressure,
      computed_gradient,
      theoretical_gradient,
      computed_laplacian);
  
  std::cout << "\n测试完成!" << std::endl;
  
  return 0;
}
