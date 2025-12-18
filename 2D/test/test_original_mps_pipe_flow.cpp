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

// 速度梯度结构（用于输出）
struct VelocityGradient {
  double grad_xx;  // ∂v_x/∂x
  double grad_xy;  // ∂v_x/∂y
  double grad_yx;  // ∂v_y/∂x
  double grad_yy;  // ∂v_y/∂y
};

// 速度拉普拉斯算子结构
struct VelocityLaplacian {
  double laplacian_x;  // ∂^2v_x/∂x^2 + ∂^2v_x/∂y^2
  double laplacian_y;  // ∂^2v_y/∂x^2 + ∂^2v_y/∂y^2
};

// 从原始MPS梯度计算速度梯度张量
// 原始MPS方法计算的是标量场的梯度，对于向量场需要对每个分量分别计算
VelocityGradient ComputeVelocityGradient(
    int particle_idx,
    const std::vector<double2>& vector_field,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double smoothing_radius,
    double reference_density,
    OriginalMPS& mps_calculator) {
  
  VelocityGradient grad;
  
  // 提取速度的x和y分量
  std::vector<double> vx_field(fluid_particles.particle_num);
  std::vector<double> vy_field(fluid_particles.particle_num);
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    vx_field[i] = vector_field[i].x;
    vy_field[i] = vector_field[i].y;
  }
  
  // 计算v_x的梯度
  double2 grad_vx = mps_calculator.ComputeGradient(
      particle_idx, vx_field, fluid_particles, solid_particles,
      smoothing_radius, reference_density);
  
  // 计算v_y的梯度
  double2 grad_vy = mps_calculator.ComputeGradient(
      particle_idx, vy_field, fluid_particles, solid_particles,
      smoothing_radius, reference_density);
  
  grad.grad_xx = grad_vx.x;  // ∂v_x/∂x
  grad.grad_xy = grad_vx.y;  // ∂v_x/∂y
  grad.grad_yx = grad_vy.x;  // ∂v_y/∂x
  grad.grad_yy = grad_vy.y;  // ∂v_y/∂y
  
  return grad;
}

// 从原始MPS拉普拉斯算子计算速度拉普拉斯算子
VelocityLaplacian ComputeVelocityLaplacian(
    int particle_idx,
    const std::vector<double2>& vector_field,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double smoothing_radius,
    double reference_density,
    double particle_spacing,
    OriginalMPS& mps_calculator) {
  
  VelocityLaplacian laplacian;
  
  // 提取速度的x和y分量
  std::vector<double> vx_field(fluid_particles.particle_num);
  std::vector<double> vy_field(fluid_particles.particle_num);
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    vx_field[i] = vector_field[i].x;
    vy_field[i] = vector_field[i].y;
  }
  
  // 计算v_x的拉普拉斯算子（使用离散lambda）
  laplacian.laplacian_x = mps_calculator.ComputeLaplacian(
      particle_idx, vx_field, fluid_particles, solid_particles,
      smoothing_radius, reference_density, particle_spacing);
  
  // 计算v_y的拉普拉斯算子（使用离散lambda）
  laplacian.laplacian_y = mps_calculator.ComputeLaplacian(
      particle_idx, vy_field, fluid_particles, solid_particles,
      smoothing_radius, reference_density, particle_spacing);
  
  return laplacian;
}

// 生成管道流动测试场景
void GeneratePipeFlowTest(
    double channel_length,
    double channel_height,
    double particle_spacing,
    FluidParticle& fluid_particles,
    SolidParticle& solid_particles) {
  
  // 计算粒子数量
  int nx_fluid = static_cast<int>(channel_length / particle_spacing) + 1;
  int ny_fluid = static_cast<int>(channel_height / particle_spacing) - 1;  // 排除上下壁面
  int num_fluid = nx_fluid * ny_fluid;
  
  // 生成流体粒子
  fluid_particles.particle_num = num_fluid;
  fluid_particles.position.resize(num_fluid);
  fluid_particles.velocity.resize(num_fluid);
  fluid_particles.pressure.resize(num_fluid);
  fluid_particles.density.resize(num_fluid);
  fluid_particles.surface_type.resize(num_fluid);
  fluid_particles.fluid_neighbour_list.resize(num_fluid);
  fluid_particles.solid_neighbour_list.resize(num_fluid);
  
  int idx = 0;
  for (int j = 1; j <= ny_fluid; ++j) {
    for (int i = 0; i < nx_fluid; ++i) {
      double x = i * particle_spacing;
      double y = j * particle_spacing;
      fluid_particles.position[idx] = {x, y};
      fluid_particles.velocity[idx] = {0.0, 0.0};  // 稍后设置速度
      fluid_particles.pressure[idx] = 0.0;
      fluid_particles.density[idx] = 1000.0;
      fluid_particles.surface_type[idx] = SurfaceType::INNER;
      ++idx;
    }
  }
  
  // 生成固体粒子（上下壁面）
  int nx_bottom = static_cast<int>(channel_length / particle_spacing) + 1;
  int nx_top = nx_bottom;
  
  int num_solid = nx_bottom + nx_top;
  solid_particles.particle_num = num_solid;
  solid_particles.position.resize(num_solid);
  solid_particles.velocity.resize(num_solid);
  solid_particles.normal_vector.resize(num_solid);
  
  idx = 0;
  
  // 下壁面（法向量向上）
  for (int i = 0; i < nx_bottom; ++i) {
    double x = i * particle_spacing;
    double y = 0;
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {0.0, 1.0};  // 向上
    ++idx;
  }
  
  // 上壁面（法向量向下）
  for (int i = 0; i < nx_top; ++i) {
    double x = i * particle_spacing;
    double y = channel_height;
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {0.0, -1.0};  // 向下
    ++idx;
  }
}

// 计算理论速度分布（Poiseuille流动）
double2 ComputeTheoreticalVelocity(double y, double channel_height, double v_max) {
  double2 velocity;
  if (y < 0.0) y = 0.0;
  if (y > channel_height) y = channel_height;
  double normalized_y = 2.0 * y / channel_height - 1.0;  // 归一化到[-1, 1]
  velocity.x = v_max * (1.0 - normalized_y * normalized_y);
  velocity.y = 0.0;
  return velocity;
}

// 计算理论速度梯度
VelocityGradient ComputeTheoreticalVelocityGradient(
    double y, double channel_height, double v_max) {
  VelocityGradient grad;
  grad.grad_xx = 0.0;
  if (y < 0.0) y = 0.0;
  if (y > channel_height) y = channel_height;
  double normalized_y = 2.0 * y / channel_height - 1.0;
  grad.grad_xy = -4.0 * v_max * normalized_y / channel_height;
  grad.grad_yx = 0.0;
  grad.grad_yy = 0.0;
  return grad;
}

// 计算理论散度（不可压缩流动，散度应为0）
double ComputeTheoreticalDivergence() {
  return 0.0;
}

// 计算理论速度拉普拉斯算子
VelocityLaplacian ComputeTheoreticalVelocityLaplacian(
    double y, double channel_height, double v_max) {
  VelocityLaplacian laplacian;
  laplacian.laplacian_x = -8.0 * v_max / (channel_height * channel_height);
  laplacian.laplacian_y = 0.0;
  return laplacian;
}

// 输出速度数据到VTK文件
void WriteVelocityToVTK(
    const std::string& filename,
    const FluidParticle& fluid_particles,
    const std::vector<double2>& theoretical_velocity,
    const std::vector<double2>& computed_velocity,
    const std::vector<VelocityGradient>& computed_gradient,
    const std::vector<VelocityGradient>& theoretical_gradient,
    const std::vector<double>& computed_divergence,
    const std::vector<double>& theoretical_divergence,
    const std::vector<VelocityLaplacian>& computed_laplacian,
    const std::vector<VelocityLaplacian>& theoretical_laplacian) {
  
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "错误：无法打开文件 " << filename << std::endl;
    return;
  }
  
  int num_particles = fluid_particles.particle_num;
  
  file << std::fixed << std::setprecision(15);
  
  // VTK文件头
  file << "# vtk DataFile Version 3.0\n";
  file << "Pipe Flow Test - Original MPS 2D\n";
  file << "ASCII\n";
  file << "DATASET POLYDATA\n";
  
  // 写入点坐标
  file << "POINTS " << num_particles << " float\n";
  for (int i = 0; i < num_particles; ++i) {
    const auto& pos = fluid_particles.position[i];
    file << pos.x << " " << pos.y << " 0.0\n";
  }
  
  // 写入顶点
  file << "VERTICES " << num_particles << " " << (num_particles * 2) << "\n";
  for (int i = 0; i < num_particles; ++i) {
    file << "1 " << i << "\n";
  }
  
  // 写入点数据
  file << "POINT_DATA " << num_particles << "\n";
  
  // 写入理论速度
  file << "VECTORS theoretical_velocity float\n";
  for (int i = 0; i < num_particles; ++i) {
    file << theoretical_velocity[i].x << " " 
         << theoretical_velocity[i].y << " 0.0\n";
  }
  
  // 写入计算出的速度
  file << "VECTORS computed_velocity float\n";
  for (int i = 0; i < num_particles; ++i) {
    file << computed_velocity[i].x << " " 
         << computed_velocity[i].y << " 0.0\n";
  }
  
  // 写入速度梯度（计算值）
  file << "TENSORS computed_velocity_gradient float\n";
  for (int i = 0; i < num_particles; ++i) {
    file << computed_gradient[i].grad_xx << " " 
         << computed_gradient[i].grad_xy << " 0.0\n";
    file << computed_gradient[i].grad_yx << " " 
         << computed_gradient[i].grad_yy << " 0.0\n";
    file << "0.0 0.0 0.0\n";
  }
  
  // 写入速度梯度（理论值）
  file << "TENSORS theoretical_velocity_gradient float\n";
  for (int i = 0; i < num_particles; ++i) {
    file << theoretical_gradient[i].grad_xx << " " 
         << theoretical_gradient[i].grad_xy << " 0.0\n";
    file << theoretical_gradient[i].grad_yx << " " 
         << theoretical_gradient[i].grad_yy << " 0.0\n";
    file << "0.0 0.0 0.0\n";
  }
  
  // 写入速度梯度误差
  file << "TENSORS velocity_gradient_error float\n";
  for (int i = 0; i < num_particles; ++i) {
    double error_xx = computed_gradient[i].grad_xx - theoretical_gradient[i].grad_xx;
    double error_xy = computed_gradient[i].grad_xy - theoretical_gradient[i].grad_xy;
    double error_yx = computed_gradient[i].grad_yx - theoretical_gradient[i].grad_yx;
    double error_yy = computed_gradient[i].grad_yy - theoretical_gradient[i].grad_yy;
    file << error_xx << " " << error_xy << " 0.0\n";
    file << error_yx << " " << error_yy << " 0.0\n";
    file << "0.0 0.0 0.0\n";
  }
  
  // 写入散度（计算值）
  file << "SCALARS computed_divergence float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    file << computed_divergence[i] << "\n";
  }
  
  // 写入散度（理论值）
  file << "SCALARS theoretical_divergence float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    file << theoretical_divergence[i] << "\n";
  }
  
  // 写入散度误差
  file << "SCALARS divergence_error float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    double error = computed_divergence[i] - theoretical_divergence[i];
    file << error << "\n";
  }

  // 写入拉普拉斯算子（计算值）- 向量形式
  file << "VECTORS computed_velocity_laplacian float\n";
  for (int i = 0; i < num_particles; ++i) {
    file << computed_laplacian[i].laplacian_x << " " 
         << computed_laplacian[i].laplacian_y << " 0.0\n";
  }
  
  // 写入拉普拉斯算子（理论值）- 向量形式
  file << "VECTORS theoretical_velocity_laplacian float\n";
  for (int i = 0; i < num_particles; ++i) {
    file << theoretical_laplacian[i].laplacian_x << " " 
         << theoretical_laplacian[i].laplacian_y << " 0.0\n";
  }
  
  // 写入拉普拉斯算子误差
  file << "VECTORS laplacian_error float\n";
  for (int i = 0; i < num_particles; ++i) {
    double error_x = computed_laplacian[i].laplacian_x - theoretical_laplacian[i].laplacian_x;
    double error_y = computed_laplacian[i].laplacian_y - theoretical_laplacian[i].laplacian_y;
    file << error_x << " " << error_y << " 0.0\n";
  }
  
  // 写入拉普拉斯算子误差大小
  file << "SCALARS laplacian_error_magnitude float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    double error_x = computed_laplacian[i].laplacian_x - theoretical_laplacian[i].laplacian_x;
    double error_y = computed_laplacian[i].laplacian_y - theoretical_laplacian[i].laplacian_y;
    double error_mag = std::sqrt(error_x * error_x + error_y * error_y);
    file << error_mag << "\n";
  }
  
  // 写入Y坐标
  file << "SCALARS y_coordinate float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    file << fluid_particles.position[i].y << "\n";
  }
  
  file.close();
  std::cout << "  已输出VTK文件: " << filename << std::endl;
}

int main() {
  std::cout << "=== 管道流动测例（原始MPS方法） ===" << std::endl;
  
  // 几何参数
  const double channel_length = 2.0;   // 管道长度 (m)
  const double channel_height = 0.5;   // 管道高度 (m)
  
  // 速度参数
  const double v_max = 1.0;  // 最大速度（管道中心）(m/s)
  
  // 粒子参数
  const double particle_spacing = 0.02;  // 粒子间距 (m)
  const double particle_radius = particle_spacing / 2.0;
  const double smoothing_radius = 2.1 * particle_spacing;
  const double cell_size = 2.0 * smoothing_radius;
  
  std::cout << "\n几何参数:" << std::endl;
  std::cout << "  管道长度: " << channel_length << " m" << std::endl;
  std::cout << "  管道高度: " << channel_height << " m" << std::endl;
  
  std::cout << "\n速度参数:" << std::endl;
  std::cout << "  最大速度: " << v_max << " m/s" << std::endl;
  
  std::cout << "\n粒子参数:" << std::endl;
  std::cout << "  粒子间距: " << particle_spacing << " m" << std::endl;
  std::cout << "  平滑半径: " << smoothing_radius << " m" << std::endl;
  
  // 生成测试场景
  FluidParticle fluid_particles("fluid");
  SolidParticle solid_particles("solid");
  
  std::cout << "\n生成测试场景..." << std::endl;
  GeneratePipeFlowTest(
      channel_length, channel_height, particle_spacing,
      fluid_particles, solid_particles);
  
  std::cout << "  流体粒子数: " << fluid_particles.particle_num << std::endl;
  std::cout << "  固体粒子数: " << solid_particles.particle_num << std::endl;
  
  // 构建邻居列表
  std::cout << "\n构建邻居列表..." << std::endl;
  NeighborListSearcher neighbor_searcher;
  neighbor_searcher.BuildNeighborList(
      fluid_particles, solid_particles,
      particle_radius, smoothing_radius, cell_size);
  
  // 设置理论速度场
  std::cout << "\n设置理论速度场..." << std::endl;
  std::vector<double2> theoretical_velocity(fluid_particles.particle_num);
  std::vector<double2> computed_velocity(fluid_particles.particle_num);
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    double y = fluid_particles.position[i].y;
    theoretical_velocity[i] = ComputeTheoreticalVelocity(y, channel_height, v_max);
    computed_velocity[i] = theoretical_velocity[i];  // 使用理论值作为输入
    fluid_particles.velocity[i] = theoretical_velocity[i];
  }
  
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
  
  // 计算速度梯度
  std::cout << "\n计算速度梯度（原始MPS方法）..." << std::endl;
  std::vector<VelocityGradient> computed_gradient(fluid_particles.particle_num);
  std::vector<VelocityGradient> theoretical_gradient(fluid_particles.particle_num);
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    computed_gradient[i] = ComputeVelocityGradient(
        i, computed_velocity, fluid_particles, solid_particles,
        smoothing_radius, reference_density, mps_calculator);
    
    double y = fluid_particles.position[i].y;
    theoretical_gradient[i] = ComputeTheoreticalVelocityGradient(
        y, channel_height, v_max);
  }
  
  // 计算散度
  std::cout << "\n计算速度散度（原始MPS方法）..." << std::endl;
  std::vector<double> computed_divergence(fluid_particles.particle_num);
  std::vector<double> theoretical_divergence(fluid_particles.particle_num);
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    computed_divergence[i] = mps_calculator.ComputeDivergence(
        i, computed_velocity, fluid_particles, solid_particles,
        smoothing_radius, reference_density);
    theoretical_divergence[i] = ComputeTheoreticalDivergence();
  }

  // 计算速度拉普拉斯算子
  std::cout << "\n计算速度拉普拉斯算子（原始MPS方法，使用离散lambda）..." << std::endl;
  std::vector<VelocityLaplacian> computed_laplacian(fluid_particles.particle_num);
  std::vector<VelocityLaplacian> theoretical_laplacian(fluid_particles.particle_num);
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    computed_laplacian[i] = ComputeVelocityLaplacian(
        i, computed_velocity, fluid_particles, solid_particles,
        smoothing_radius, reference_density, particle_spacing, mps_calculator);
    
    double y = fluid_particles.position[i].y;
    theoretical_laplacian[i] = ComputeTheoreticalVelocityLaplacian(
        y, channel_height, v_max);
  }
  
  // 统计误差（相对误差）
  std::cout << "\n速度梯度误差统计（相对误差）:" << std::endl;
  double total_relative_error_xx = 0.0, total_relative_error_xy = 0.0;
  double total_relative_error_yx = 0.0, total_relative_error_yy = 0.0;
  double max_relative_error_xx = 0.0, max_relative_error_xy = 0.0;
  double max_relative_error_yx = 0.0, max_relative_error_yy = 0.0;
  int num_valid = fluid_particles.particle_num;
  int num_valid_xx = 0, num_valid_xy = 0, num_valid_yx = 0, num_valid_yy = 0;
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    // ∂v_x/∂x 相对误差
    double abs_error_xx = std::abs(computed_gradient[i].grad_xx - theoretical_gradient[i].grad_xx);
    double relative_error_xx = (std::abs(theoretical_gradient[i].grad_xx) > 1e-10) ?
        (abs_error_xx / std::abs(theoretical_gradient[i].grad_xx) * 100.0) : abs_error_xx;
    total_relative_error_xx += relative_error_xx;
    if (relative_error_xx > max_relative_error_xx) max_relative_error_xx = relative_error_xx;
    ++num_valid_xx;
    
    // ∂v_x/∂y 相对误差
    double abs_error_xy = std::abs(computed_gradient[i].grad_xy - theoretical_gradient[i].grad_xy);
    double relative_error_xy = (std::abs(theoretical_gradient[i].grad_xy) > 1e-10) ?
        (abs_error_xy / std::abs(theoretical_gradient[i].grad_xy) * 100.0) : abs_error_xy;
    total_relative_error_xy += relative_error_xy;
    if (relative_error_xy > max_relative_error_xy) max_relative_error_xy = relative_error_xy;
    ++num_valid_xy;
    
    // ∂v_y/∂x 相对误差
    double abs_error_yx = std::abs(computed_gradient[i].grad_yx - theoretical_gradient[i].grad_yx);
    double relative_error_yx = (std::abs(theoretical_gradient[i].grad_yx) > 1e-10) ?
        (abs_error_yx / std::abs(theoretical_gradient[i].grad_yx) * 100.0) : abs_error_yx;
    total_relative_error_yx += relative_error_yx;
    if (relative_error_yx > max_relative_error_yx) max_relative_error_yx = relative_error_yx;
    ++num_valid_yx;
    
    // ∂v_y/∂y 相对误差
    double abs_error_yy = std::abs(computed_gradient[i].grad_yy - theoretical_gradient[i].grad_yy);
    double relative_error_yy = (std::abs(theoretical_gradient[i].grad_yy) > 1e-10) ?
        (abs_error_yy / std::abs(theoretical_gradient[i].grad_yy) * 100.0) : abs_error_yy;
    total_relative_error_yy += relative_error_yy;
    if (relative_error_yy > max_relative_error_yy) max_relative_error_yy = relative_error_yy;
    ++num_valid_yy;
  }
  
  double avg_relative_error_xx = total_relative_error_xx / num_valid_xx;
  double avg_relative_error_xy = total_relative_error_xy / num_valid_xy;
  double avg_relative_error_yx = total_relative_error_yx / num_valid_yx;
  double avg_relative_error_yy = total_relative_error_yy / num_valid_yy;
  
  std::cout << "  统计粒子数: " << num_valid << std::endl;
  std::cout << "  ∂v_x/∂x 平均相对误差: " << std::fixed << std::setprecision(4) 
            << avg_relative_error_xx << "%" << std::endl;
  std::cout << "  ∂v_x/∂y 平均相对误差: " << avg_relative_error_xy << "%" << std::endl;
  std::cout << "  ∂v_y/∂x 平均相对误差: " << avg_relative_error_yx << "%" << std::endl;
  std::cout << "  ∂v_y/∂y 平均相对误差: " << avg_relative_error_yy << "%" << std::endl;
  std::cout << "  ∂v_x/∂x 最大相对误差: " << max_relative_error_xx << "%" << std::endl;
  std::cout << "  ∂v_x/∂y 最大相对误差: " << max_relative_error_xy << "%" << std::endl;
  std::cout << "  ∂v_y/∂x 最大相对误差: " << max_relative_error_yx << "%" << std::endl;
  std::cout << "  ∂v_y/∂y 最大相对误差: " << max_relative_error_yy << "%" << std::endl;
  
  // 散度误差统计（理论值为0，使用绝对误差）
  std::cout << "\n速度散度误差统计（理论值为0，显示绝对误差）:" << std::endl;
  double total_error_div = 0.0;
  double max_error_div = 0.0;
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    double error = std::abs(computed_divergence[i] - theoretical_divergence[i]);
    total_error_div += error;
    if (error > max_error_div) max_error_div = error;
  }
  
  double avg_error_div = total_error_div / num_valid;
  std::cout << "  平均散度误差: " << std::fixed << std::setprecision(6) 
            << avg_error_div << " 1/s" << std::endl;
  std::cout << "  最大散度误差: " << max_error_div << " 1/s" << std::endl;
  std::cout << "  理论散度: 0.0 1/s（不可压缩流动，无法计算相对误差）" << std::endl;
  
  // 拉普拉斯算子误差统计（相对误差）
  std::cout << "\n速度拉普拉斯算子误差统计（相对误差）:" << std::endl;
  double total_relative_error_lap_x = 0.0, total_relative_error_lap_y = 0.0;
  double max_relative_error_lap_x = 0.0, max_relative_error_lap_y = 0.0;
  int num_valid_lap_x = 0, num_valid_lap_y = 0;
  
  // 计算理论拉普拉斯算子值（用于参考）
  double theoretical_lap_x = -8.0 * v_max / (channel_height * channel_height);
  const double theoretical_lap_y = 0.0;
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    // ∇²v_x 相对误差
    double abs_error_x = std::abs(computed_laplacian[i].laplacian_x - theoretical_laplacian[i].laplacian_x);
    double relative_error_x = (std::abs(theoretical_lap_x) > 1e-10) ?
        (abs_error_x / std::abs(theoretical_lap_x) * 100.0) : abs_error_x;
    total_relative_error_lap_x += relative_error_x;
    if (relative_error_x > max_relative_error_lap_x) max_relative_error_lap_x = relative_error_x;
    ++num_valid_lap_x;
    
    // ∇²v_y 相对误差（理论值为0，使用绝对误差）
    double abs_error_y = std::abs(computed_laplacian[i].laplacian_y - theoretical_laplacian[i].laplacian_y);
    double relative_error_y = (std::abs(theoretical_lap_y) > 1e-10) ?
        (abs_error_y / std::abs(theoretical_lap_y) * 100.0) : abs_error_y;
    total_relative_error_lap_y += relative_error_y;
    if (relative_error_y > max_relative_error_lap_y) max_relative_error_lap_y = relative_error_y;
    ++num_valid_lap_y;
  }
  
  double avg_relative_error_lap_x = total_relative_error_lap_x / num_valid_lap_x;
  double avg_relative_error_lap_y = total_relative_error_lap_y / num_valid_lap_y;
  
  std::cout << "  统计粒子数: " << num_valid << std::endl;
  std::cout << "  ∇²v_x 平均相对误差: " << std::fixed << std::setprecision(4) 
            << avg_relative_error_lap_x << "%" << std::endl;
  std::cout << "  ∇²v_y 平均误差（理论值为0，显示绝对误差）: " 
            << std::setprecision(6) << avg_relative_error_lap_y << " 1/(m²·s)" << std::endl;
  std::cout << "  ∇²v_x 最大相对误差: " << std::setprecision(4) 
            << max_relative_error_lap_x << "%" << std::endl;
  std::cout << "  ∇²v_y 最大误差（理论值为0，显示绝对误差）: " 
            << std::setprecision(6) << max_relative_error_lap_y << " 1/(m²·s)" << std::endl;
  std::cout << "  理论 ∇²v_x: " << std::setprecision(10) << theoretical_lap_x << " 1/(m²·s)" << std::endl;
  std::cout << "  理论 ∇²v_y: " << theoretical_lap_y << " 1/(m²·s)" << std::endl;
  
  // 输出到VTK文件
  std::cout << "\n输出VTK文件..." << std::endl;
  WriteVelocityToVTK(
      "original_mps_pipe_flow.vtk",
      fluid_particles,
      theoretical_velocity,
      computed_velocity,
      computed_gradient,
      theoretical_gradient,
      computed_divergence,
      theoretical_divergence,
      computed_laplacian,
      theoretical_laplacian);
  
  std::cout << "\n测试完成!" << std::endl;
  
  return 0;
}
