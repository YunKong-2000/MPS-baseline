#include "../src/lsmps/CorrectiveMatrix.hpp"
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

using namespace mps2D;
void printMatrix(const Eigen::Matrix<double, 5, 5>& matrix, std::string name) {
  std::cout << name << ":" << std::endl;
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      std::cout << std::setw(10) << matrix(i, j) << " ";
    }
    std::cout << std::endl;
  }
}

// 计算向量场的梯度（使用LSMPS corrective matrix方法）
// 参数：
//   particle_idx: 粒子索引
//   vector_field: 向量场值（速度）
//   fluid_particles: 流体粒子对象
//   solid_particles: 固体粒子对象
//   corrective_matrix: corrective matrix (5x5)
//   smoothing_radius: 平滑半径
// 返回：梯度张量的分量 (grad_xx, grad_xy, grad_yx, grad_yy)
struct VelocityGradient {
  double grad_xx;  // ∂v_x/∂x
  double grad_xy;  // ∂v_x/∂y
  double grad_yx;  // ∂v_y/∂x
  double grad_yy;  // ∂v_y/∂y
};

VelocityGradient ComputeVelocityGradient(
    int particle_idx,
    const std::vector<double2>& vector_field,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const Eigen::Matrix<double, 5, 5>& corrective_matrix,
    double smoothing_radius) {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  const double2& v_i = vector_field[particle_idx];
  
  // 提取corrective matrix的前两行（C1和C2）
  Eigen::Matrix<double, 1, 5> C1 = corrective_matrix.row(0);  // x方向
  Eigen::Matrix<double, 1, 5> C2 = corrective_matrix.row(1);  // y方向
  
  VelocityGradient grad;
  grad.grad_xx = 0.0;  // ∂v_x/∂x
  grad.grad_xy = 0.0;  // ∂v_x/∂y
  grad.grad_yx = 0.0;  // ∂v_y/∂x
  grad.grad_yy = 0.0;  // ∂v_y/∂y
  
  // 使用与 LSMPS moment 矩阵相同的基函数定义
  // P_ij = [x/r_e, y/r_e, x^2/r_e^2, y^2/r_e^2, x y / r_e^2]
  CorrectiveMatrix corrective_matrix_calc;
  
  // 处理流体邻域粒子
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    
    if (dist < 1e-10 || dist > smoothing_radius) continue;
    
    const double2& v_j = vector_field[j];
    double dv_x = v_j.x - v_i.x;  // 直接使用速度差，不除以dist
    double dv_y = v_j.y - v_i.y;
    double weight = WeightFunction(dist, smoothing_radius);
    
    // 计算基函数 P_ij（使用标准LSMPS基函数）
    Eigen::Vector<double, 5> basis =
        corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
    
    // 根据文档：∇φ_i = (1/r_e) Σ_j w_ij (φ_j - φ_i) [M_{i,0}; M_{i,1}] P_ij
    double M0P = (C1 * basis)(0, 0);  // M_{i,0} * P_ij
    double M1P = (C2 * basis)(0, 0);  // M_{i,1} * P_ij
    
    grad.grad_xx += weight * dv_x * M0P;  // ∂v_x/∂x
    grad.grad_xy += weight * dv_x * M1P;  // ∂v_x/∂y
    grad.grad_yx += weight * dv_y * M0P;  // ∂v_y/∂x
    grad.grad_yy += weight * dv_y * M1P;  // ∂v_y/∂y
  }
  
  // 处理固体邻域粒子（壁面粒子）
  // 对于管道流动，壁面处速度为零（无滑移边界条件，第一类边界条件）
  // 第一类边界条件使用标准基函数，而不是壁面基函数
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    const double2& pos_j = solid_particles.position[j];
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    
    if (dist < 1e-10 || dist > smoothing_radius) continue;
    
    // 壁面处速度为零（第一类边界条件：Dirichlet边界条件）
    const double2 v_wall = {0.0, 0.0};
    double dv_x = v_wall.x - v_i.x;  // 直接使用速度差，不除以dist
    double dv_y = v_wall.y - v_i.y;
    double weight = WeightFunction(dist, smoothing_radius);
    
    // 第一类边界条件使用标准基函数（与流体粒子相同）
    Eigen::Vector<double, 5> basis =
        corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
    
    // 计算梯度贡献
    double M0P = (C1 * basis)(0, 0);
    double M1P = (C2 * basis)(0, 0);
    
    grad.grad_xx += weight * dv_x * M0P;
    grad.grad_xy += weight * dv_x * M1P;
    grad.grad_yx += weight * dv_y * M0P;
    grad.grad_yy += weight * dv_y * M1P;
  }
  
  // 统一乘以 1 / r_e（根据文档公式）
  double inv_re = 1.0 / smoothing_radius;
  grad.grad_xx *= inv_re;
  grad.grad_xy *= inv_re;
  grad.grad_yx *= inv_re;
  grad.grad_yy *= inv_re;
  
  return grad;
}

// 计算向量场的拉普拉斯算子（使用LSMPS corrective matrix方法）
// 参数：
//   particle_idx: 粒子索引
//   vector_field: 向量场值（速度）
//   fluid_particles: 流体粒子对象
//   solid_particles: 固体粒子对象
//   corrective_matrix: corrective matrix (5x5)
//   smoothing_radius: 平滑半径
// 返回：拉普拉斯算子

struct VelocityLaplacian {
  double laplacian_x;  // ∂^2v_x/∂x^2 + ∂^2v_x/∂y^2
  double laplacian_y;  // ∂^2v_y/∂x^2 + ∂^2v_y/∂y^2
};

VelocityLaplacian ComputeVelocityLaplacian(
  int particle_idx,
  const std::vector<double2>& vector_field,
  const FluidParticle& fluid_particles,
  const SolidParticle& solid_particles,
  const Eigen::Matrix<double, 5, 5>& corrective_matrix,
  double smoothing_radius) {

const double2& pos_i = fluid_particles.position[particle_idx];
const double2& v_i = vector_field[particle_idx];

// 提取 moment matrix 逆矩阵的第3行和第4行（用于拉普拉斯算子计算）
// 文档中对应 [M_{i,2} + M_{i,3}]
Eigen::RowVector<double, 5> M2 = corrective_matrix.row(2);  // x² 项
Eigen::RowVector<double, 5> M3 = corrective_matrix.row(3);  // y² 项
Eigen::RowVector<double, 5> M2_plus_M3 = M2 + M3;

VelocityLaplacian laplacian;
laplacian.laplacian_x = 0.0;
laplacian.laplacian_y = 0.0;

// 标量因子：2 / r_e^2（根据文档公式）
const double scalar_factor = 2.0 / (smoothing_radius * smoothing_radius);

// 使用与 LSMPS moment 矩阵相同的基函数定义
CorrectiveMatrix corrective_matrix_calc;

// 处理流体邻域粒子
for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
  const double2& pos_j = fluid_particles.position[j];
  double dx = pos_j.x - pos_i.x;
  double dy = pos_j.y - pos_i.y;
  double dist = ComputeDistance(pos_i, pos_j);
  
  if (dist < 1e-10 || dist > smoothing_radius) continue;
  
  const double2& v_j = vector_field[j];
  double dv_x = v_j.x - v_i.x;  // 直接使用速度差，不除以dist
  double dv_y = v_j.y - v_i.y;
  double weight = WeightFunction(dist, smoothing_radius);
  
  // 计算基函数 P_ij（使用标准LSMPS基函数）
  Eigen::Vector<double, 5> basis =
      corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
  
  // 根据文档：∇²φ_i = (2/r_e^2) Σ_j w_ij (φ_j - φ_i) [M_{i,2} + M_{i,3}] P_ij
  double M2M3P = (M2_plus_M3 * basis)(0, 0);
  double coeff = scalar_factor * weight * M2M3P;
  
  laplacian.laplacian_x += coeff * dv_x;  // ∂^2v_x/∂x^2 + ∂^2v_x/∂y^2
  laplacian.laplacian_y += coeff * dv_y;  // ∂^2v_y/∂x^2 + ∂^2v_y/∂y^2
}

// 处理固体邻域粒子（壁面粒子）
// 对于管道流动，壁面处速度为零（无滑移边界条件，第一类边界条件）
// 第一类边界条件使用标准基函数，而不是壁面基函数
for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
  const double2& pos_j = solid_particles.position[j];
  double dx = pos_j.x - pos_i.x;
  double dy = pos_j.y - pos_i.y;
  double dist = ComputeDistance(pos_i, pos_j);
  
  if (dist < 1e-10 || dist > smoothing_radius) continue;
  
  // 壁面处速度为零（第一类边界条件：Dirichlet边界条件）
  const double2 v_wall = {0.0, 0.0};
  double dv_x = v_wall.x - v_i.x;  // 直接使用速度差，不除以dist
  double dv_y = v_wall.y - v_i.y;
  double weight = WeightFunction(dist, smoothing_radius);
  
  // 第一类边界条件使用标准基函数（与流体粒子相同）
  Eigen::Vector<double, 5> basis =
      corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
  
  // 计算拉普拉斯算子贡献
  double M2M3P = (M2_plus_M3 * basis)(0, 0);
  double coeff = scalar_factor * weight * M2M3P;
  
  laplacian.laplacian_x += coeff * dv_x;
  laplacian.laplacian_y += coeff * dv_y;
}

return laplacian;
}

// 计算速度散度（使用LSMPS corrective matrix方法）
// 散度 = ∂v_x/∂x + ∂v_y/∂y
double ComputeVelocityDivergence(
    int particle_idx,
    const std::vector<double2>& vector_field,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const Eigen::Matrix<double, 5, 5>& corrective_matrix,
    double smoothing_radius) {
  
  VelocityGradient grad = ComputeVelocityGradient(
      particle_idx, vector_field, fluid_particles, solid_particles,
      corrective_matrix, smoothing_radius);
  
  return grad.grad_xx + grad.grad_yy;  // div(v) = ∂v_x/∂x + ∂v_y/∂y
}

// 生成管道流动测试场景
// 参数：
//   channel_length: 管道长度
//   channel_height: 管道高度
//   particle_spacing: 粒子间距
//   fluid_particles: 输出流体粒子
//   solid_particles: 输出固体粒子（壁面）
void GeneratePipeFlowTest(
    double channel_length,
    double channel_height,
    double particle_spacing,
    FluidParticle& fluid_particles,
    SolidParticle& solid_particles) {
  
  // 计算粒子数量
  // 流体粒子应该在 y = particle_spacing 到 y = channel_height - particle_spacing
  // 壁面在 y = 0 和 y = channel_height
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
  // 流体粒子从 y = particle_spacing 开始，到 y = (ny_fluid) * particle_spacing 结束
  // 确保最后一个粒子距离上壁面为 particle_spacing
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
  // 下壁面
  int nx_bottom = static_cast<int>(channel_length / particle_spacing) + 1;
  // 上壁面
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
    double y = 0;  // 距离底部流体粒子一个粒子间距
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {0.0, 1.0};  // 向上
    ++idx;
  }
  
  // 上壁面（法向量向下）
  for (int i = 0; i < nx_top; ++i) {
    double x = i * particle_spacing;
    double y = channel_height;  // 距离顶部流体粒子一个粒子间距
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {0.0, -1.0};  // 向下
    ++idx;
  }
}

// 计算理论速度分布（Poiseuille流动）
// v_x(y) = v_max * (1 - (2y/h - 1)^2)
// 注意：这个公式假设 y 的范围是 [0, h]，其中 y=0 和 y=h 是壁面（速度=0）
// v_y = 0
// 参数：
//   y: y坐标
//   channel_height: 管道高度
//   v_max: 最大速度（管道中心）
// 返回：速度向量
double2 ComputeTheoreticalVelocity(double y, double channel_height, double v_max) {
  double2 velocity;
  // 确保 y 在有效范围内
  if (y < 0.0) y = 0.0;
  if (y > channel_height) y = channel_height;
  double normalized_y = 2.0 * y / channel_height - 1.0;  // 归一化到[-1, 1]
  velocity.x = v_max * (1.0 - normalized_y * normalized_y);
  velocity.y = 0.0;
  return velocity;
}

// 计算理论速度梯度
// ∂v_x/∂x = 0
// ∂v_x/∂y = -4*v_max*(2y/h - 1)/h
// ∂v_y/∂x = 0
// ∂v_y/∂y = 0
// 注意：在边界 y=0 和 y=h 处，梯度应该考虑边界条件
VelocityGradient ComputeTheoreticalVelocityGradient(
    double y, double channel_height, double v_max) {
  VelocityGradient grad;
  grad.grad_xx = 0.0;
  // 确保 y 在有效范围内
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
// 对于Poiseuille流动：v_x(y) = v_max * (1 - (2y/h - 1)^2)
// ∂^2v_x/∂x^2 = 0
// ∂^2v_x/∂y^2 = -8*v_max/h^2
// ∂^2v_y/∂x^2 = 0
// ∂^2v_y/∂y^2 = 0
VelocityLaplacian ComputeTheoreticalVelocityLaplacian(
    double y, double channel_height, double v_max) {
  VelocityLaplacian laplacian;
  laplacian.laplacian_x = -8.0 * v_max / (channel_height * channel_height);  // ∂^2v_x/∂x^2 + ∂^2v_x/∂y^2 = 0 + (-8*v_max/h^2)
  laplacian.laplacian_y = 0.0;  // ∂^2v_y/∂x^2 + ∂^2v_y/∂y^2 = 0 + 0
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
  file << "Pipe Flow Test - MPS 2D\n";
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
  
  // 写入速度X分量（理论值）
  file << "SCALARS theoretical_velocity_x float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    file << theoretical_velocity[i].x << "\n";
  }
  
  // 写入速度X分量（计算值）
  file << "SCALARS computed_velocity_x float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    file << computed_velocity[i].x << "\n";
  }
  
  // 写入速度Y分量（理论值）
  file << "SCALARS theoretical_velocity_y float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    file << theoretical_velocity[i].y << "\n";
  }
  
  // 写入速度Y分量（计算值）
  file << "SCALARS computed_velocity_y float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    file << computed_velocity[i].y << "\n";
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
  
  // 写入拉普拉斯算子X分量（计算值）
  file << "SCALARS computed_laplacian_x float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    file << computed_laplacian[i].laplacian_x << "\n";
  }
  
  // 写入拉普拉斯算子Y分量（计算值）
  file << "SCALARS computed_laplacian_y float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    file << computed_laplacian[i].laplacian_y << "\n";
  }
  
  // 写入拉普拉斯算子（理论值）- 向量形式
  file << "VECTORS theoretical_velocity_laplacian float\n";
  for (int i = 0; i < num_particles; ++i) {
    file << theoretical_laplacian[i].laplacian_x << " " 
         << theoretical_laplacian[i].laplacian_y << " 0.0\n";
  }
  
  // 写入拉普拉斯算子X分量（理论值）
  file << "SCALARS theoretical_laplacian_x float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    file << theoretical_laplacian[i].laplacian_x << "\n";
  }
  
  // 写入拉普拉斯算子Y分量（理论值）
  file << "SCALARS theoretical_laplacian_y float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    file << theoretical_laplacian[i].laplacian_y << "\n";
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

// 输出壁面粒子数据到VTK文件
void WriteWallParticlesToVTK(
    const std::string& filename,
    const SolidParticle& solid_particles) {
  
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "错误：无法打开文件 " << filename << std::endl;
    return;
  }
  
  int num_particles = solid_particles.particle_num;
  
  file << std::fixed << std::setprecision(15);
  
  // VTK文件头
  file << "# vtk DataFile Version 3.0\n";
  file << "Wall Particles - Pipe Flow Test\n";
  file << "ASCII\n";
  file << "DATASET POLYDATA\n";
  
  // 写入点坐标（2D数据添加z=0）
  file << "POINTS " << num_particles << " float\n";
  for (int i = 0; i < num_particles; ++i) {
    const auto& pos = solid_particles.position[i];
    file << pos.x << " " << pos.y << " 0.0\n";
  }
  
  // 写入顶点（每个点作为一个顶点）
  file << "VERTICES " << num_particles << " " << (num_particles * 2) << "\n";
  for (int i = 0; i < num_particles; ++i) {
    file << "1 " << i << "\n";
  }
  
  // 写入点数据
  file << "POINT_DATA " << num_particles << "\n";
  
  // 写入法向量
  if (solid_particles.normal_vector.size() >= static_cast<size_t>(num_particles)) {
    file << "VECTORS normal_vector float\n";
    for (int i = 0; i < num_particles; ++i) {
      const auto& normal = solid_particles.normal_vector[i];
      file << normal.x << " " << normal.y << " 0.0\n";
    }
  }
  
  // 写入法向量X分量
  if (solid_particles.normal_vector.size() >= static_cast<size_t>(num_particles)) {
    file << "SCALARS normal_vector_x float\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < num_particles; ++i) {
      file << solid_particles.normal_vector[i].x << "\n";
    }
  }
  
  // 写入法向量Y分量
  if (solid_particles.normal_vector.size() >= static_cast<size_t>(num_particles)) {
    file << "SCALARS normal_vector_y float\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < num_particles; ++i) {
      file << solid_particles.normal_vector[i].y << "\n";
    }
  }
  
  // 写入X坐标
  file << "SCALARS x_coordinate float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    file << solid_particles.position[i].x << "\n";
  }
  
  // 写入Y坐标
  file << "SCALARS y_coordinate float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    file << solid_particles.position[i].y << "\n";
  }
  
  // 写入速度（如果有）
  if (solid_particles.velocity.size() >= static_cast<size_t>(num_particles)) {
    file << "VECTORS velocity float\n";
    for (int i = 0; i < num_particles; ++i) {
      const auto& vel = solid_particles.velocity[i];
      file << vel.x << " " << vel.y << " 0.0\n";
    }
  }
  
  file.close();
  std::cout << "  已输出壁面粒子VTK文件: " << filename << std::endl;
}

int main() {
  std::cout << "=== 管道流动测例（LSMPS速度梯度和散度计算） ===" << std::endl;
  
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
  
  // 判定自由面
  std::cout << "判定自由面..." << std::endl;
  SurfaceDetector surface_detector;
  surface_detector.DetectSurfaceParticles(
      fluid_particles, solid_particles,
      smoothing_radius, particle_spacing);
  
  // 设置理论速度场
  std::cout << "\n设置理论速度场..." << std::endl;
  std::vector<double2> theoretical_velocity(fluid_particles.particle_num);
  std::vector<double2> computed_velocity(fluid_particles.particle_num);
  
  // 无滑移边界条件：壁面速度为零，但最靠近壁面的流体粒子速度不为零
  // 流体粒子速度按照理论Poiseuille分布计算
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    double y = fluid_particles.position[i].y;
    theoretical_velocity[i] = ComputeTheoreticalVelocity(y, channel_height, v_max);
    computed_velocity[i] = theoretical_velocity[i];  // 使用理论值作为输入
    fluid_particles.velocity[i] = theoretical_velocity[i];
  }
  
  // 验证边界速度：最靠近壁面的流体粒子速度应该很小但不为零
  double min_y = std::numeric_limits<double>::max();
  double max_y = std::numeric_limits<double>::lowest();
  int min_y_idx = -1, max_y_idx = -1;
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    double y = fluid_particles.position[i].y;
    if (y < min_y) {
      min_y = y;
      min_y_idx = i;
    }
    if (y > max_y) {
      max_y = y;
      max_y_idx = i;
    }
  }
  if (min_y_idx >= 0 && max_y_idx >= 0) {
    std::cout << "  最靠近下壁面的粒子 (y=" << min_y << "): v_x=" 
              << theoretical_velocity[min_y_idx].x << " m/s" << std::endl;
    std::cout << "  最靠近上壁面的粒子 (y=" << max_y << "): v_x=" 
              << theoretical_velocity[max_y_idx].x << " m/s" << std::endl;
    std::cout << "  壁面速度: v_x=0.0 m/s (无滑移边界条件)" << std::endl;
  }
  
  // 计算corrective matrix
  std::cout << "\n计算corrective matrix..." << std::endl;
  CorrectiveMatrix corrective_matrix_calculator;
  std::vector<Eigen::Matrix<double, 5, 5>> matrices(fluid_particles.particle_num);
  
  int num_valid_matrix = 0;
  // 对于速度场，使用第一类边界条件（Dirichlet边界条件，壁面速度为零）
  const bool border_condition = false;  // false表示第一类边界条件
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    matrices[i] = corrective_matrix_calculator.ComputeCorrectiveMatrix(
        i, fluid_particles, solid_particles, smoothing_radius, border_condition);
    
    // 检查是否是单位矩阵
    bool is_identity = true;
    const double tolerance = 1e-6;
    for (int row = 0; row < 5; ++row) {
      for (int col = 0; col < 5; ++col) {
        double expected = (row == col) ? 1.0 : 0.0;
        if (std::abs(matrices[i](row, col) - expected) > tolerance) {
          is_identity = false;
          break;
        }
      }
      if (!is_identity) break;
    }
    if (!is_identity) ++num_valid_matrix;
  }
  std::cout << "  有效corrective matrix: " << num_valid_matrix 
            << " / " << fluid_particles.particle_num << std::endl;
  
  // 计算速度梯度
  std::cout << "\n计算速度梯度..." << std::endl;
  std::vector<VelocityGradient> computed_gradient(fluid_particles.particle_num);
  std::vector<VelocityGradient> theoretical_gradient(fluid_particles.particle_num);
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    computed_gradient[i] = ComputeVelocityGradient(
        i, computed_velocity, fluid_particles, solid_particles,
        matrices[i], smoothing_radius);
    
    double y = fluid_particles.position[i].y;
    theoretical_gradient[i] = ComputeTheoreticalVelocityGradient(
        y, channel_height, v_max);
  }
  
  // 计算散度
  std::cout << "\n计算速度散度..." << std::endl;
  std::vector<double> computed_divergence(fluid_particles.particle_num);
  std::vector<double> theoretical_divergence(fluid_particles.particle_num);
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    computed_divergence[i] = ComputeVelocityDivergence(
        i, computed_velocity, fluid_particles, solid_particles,
        matrices[i], smoothing_radius);
    theoretical_divergence[i] = ComputeTheoreticalDivergence();
  }

  // 计算速度拉普拉斯算子
  std::cout << "\n计算速度拉普拉斯算子..." << std::endl;
  std::vector<VelocityLaplacian> computed_laplacian(fluid_particles.particle_num);
  std::vector<VelocityLaplacian> theoretical_laplacian(fluid_particles.particle_num);
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    computed_laplacian[i] = ComputeVelocityLaplacian(
        i, computed_velocity, fluid_particles, solid_particles,
        matrices[i], smoothing_radius);
    
    double y = fluid_particles.position[i].y;
    theoretical_laplacian[i] = ComputeTheoreticalVelocityLaplacian(
        y, channel_height, v_max);
  }
  
  // 统计误差
  std::cout << "\n速度梯度误差统计:" << std::endl;
  double total_error_xx = 0.0, total_error_xy = 0.0;
  double total_error_yx = 0.0, total_error_yy = 0.0;
  double max_error_xx = 0.0, max_error_xy = 0.0;
  double max_error_yx = 0.0, max_error_yy = 0.0;
  int num_valid = 0;
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    // 只统计内部粒子的误差
    if (fluid_particles.surface_type[i] == SurfaceType::INNER) {
      double error_xx = std::abs(computed_gradient[i].grad_xx - theoretical_gradient[i].grad_xx);
      double error_xy = std::abs(computed_gradient[i].grad_xy - theoretical_gradient[i].grad_xy);
      double error_yx = std::abs(computed_gradient[i].grad_yx - theoretical_gradient[i].grad_yx);
      double error_yy = std::abs(computed_gradient[i].grad_yy - theoretical_gradient[i].grad_yy);
      
      total_error_xx += error_xx;
      total_error_xy += error_xy;
      total_error_yx += error_yx;
      total_error_yy += error_yy;
      
      if (error_xx > max_error_xx) max_error_xx = error_xx;
      if (error_xy > max_error_xy) max_error_xy = error_xy;
      if (error_yx > max_error_yx) max_error_yx = error_yx;
      if (error_yy > max_error_yy) max_error_yy = error_yy;
      
      ++num_valid;
    }
  }
  
  double avg_error_xx = total_error_xx / num_valid;
  double avg_error_xy = total_error_xy / num_valid;
  double avg_error_yx = total_error_yx / num_valid;
  double avg_error_yy = total_error_yy / num_valid;
  
  std::cout << "  统计粒子数（内部粒子）: " << num_valid << std::endl;
  std::cout << "  ∂v_x/∂x 平均误差: " << avg_error_xx << " 1/s" << std::endl;
  std::cout << "  ∂v_x/∂y 平均误差: " << avg_error_xy << " 1/s" << std::endl;
  std::cout << "  ∂v_y/∂x 平均误差: " << avg_error_yx << " 1/s" << std::endl;
  std::cout << "  ∂v_y/∂y 平均误差: " << avg_error_yy << " 1/s" << std::endl;
  std::cout << "  ∂v_x/∂x 最大误差: " << max_error_xx << " 1/s" << std::endl;
  std::cout << "  ∂v_x/∂y 最大误差: " << max_error_xy << " 1/s" << std::endl;
  std::cout << "  ∂v_y/∂x 最大误差: " << max_error_yx << " 1/s" << std::endl;
  std::cout << "  ∂v_y/∂y 最大误差: " << max_error_yy << " 1/s" << std::endl;
  
  std::cout << "position[24]: " << fluid_particles.position[24].x << " " << fluid_particles.position[24].y << std::endl;
  printMatrix(matrices[24], "matrices[24]");
  std::cout << "position[2526]: " << fluid_particles.position[2526].x << " " << fluid_particles.position[2526].y << std::endl;
  printMatrix(matrices[2526], "matrices[2526]");


  
  // 散度误差统计
  std::cout << "\n速度散度误差统计:" << std::endl;
  double total_error_div = 0.0;
  double max_error_div = 0.0;
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    if (fluid_particles.surface_type[i] == SurfaceType::INNER) {
      double error = std::abs(computed_divergence[i] - theoretical_divergence[i]);
      total_error_div += error;
      if (error > max_error_div) max_error_div = error;
    }
  }
  
  double avg_error_div = total_error_div / num_valid;
  std::cout << "  平均散度误差: " << avg_error_div << " 1/s" << std::endl;
  std::cout << "  最大散度误差: " << max_error_div << " 1/s" << std::endl;
  std::cout << "  理论散度: 0.0 1/s（不可压缩流动）" << std::endl;
  
  // 拉普拉斯算子误差统计
  std::cout << "\n速度拉普拉斯算子误差统计:" << std::endl;
  double total_error_lap_x = 0.0, total_error_lap_y = 0.0;
  double max_error_lap_x = 0.0, max_error_lap_y = 0.0;
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    if (fluid_particles.surface_type[i] == SurfaceType::INNER) {
      double error_x = std::abs(computed_laplacian[i].laplacian_x - theoretical_laplacian[i].laplacian_x);
      double error_y = std::abs(computed_laplacian[i].laplacian_y - theoretical_laplacian[i].laplacian_y);
      
      total_error_lap_x += error_x;
      total_error_lap_y += error_y;
      
      if (error_x > max_error_lap_x) max_error_lap_x = error_x;
      if (error_y > max_error_lap_y) max_error_lap_y = error_y;
    }
  }
  
  double avg_error_lap_x = total_error_lap_x / num_valid;
  double avg_error_lap_y = total_error_lap_y / num_valid;
  
  std::cout << "  统计粒子数（内部粒子）: " << num_valid << std::endl;
  std::cout << "  ∇²v_x 平均误差: " << avg_error_lap_x << " 1/(m²·s)" << std::endl;
  std::cout << "  ∇²v_y 平均误差: " << avg_error_lap_y << " 1/(m²·s)" << std::endl;
  std::cout << "  ∇²v_x 最大误差: " << max_error_lap_x << " 1/(m²·s)" << std::endl;
  std::cout << "  ∇²v_y 最大误差: " << max_error_lap_y << " 1/(m²·s)" << std::endl;
  
  // 计算理论拉普拉斯算子值（用于参考）
  double theoretical_lap_x = -8.0 * v_max / (channel_height * channel_height);
  std::cout << "  理论 ∇²v_x: " << theoretical_lap_x << " 1/(m²·s)" << std::endl;
  std::cout << "  理论 ∇²v_y: 0.0 1/(m²·s)" << std::endl;
  
  // 输出到VTK文件
  std::cout << "\n输出VTK文件..." << std::endl;
  WriteVelocityToVTK(
      "pipe_flow.vtk",
      fluid_particles,
      theoretical_velocity,
      computed_velocity,
      computed_gradient,
      theoretical_gradient,
      computed_divergence,
      theoretical_divergence,
      computed_laplacian,
      theoretical_laplacian);
  
  // 输出壁面粒子到VTK文件
  WriteWallParticlesToVTK(
      "pipe_flow_wall_particles.vtk",
      solid_particles);
  
  std::cout << "\n测试完成!" << std::endl;
  
  return 0;
}

