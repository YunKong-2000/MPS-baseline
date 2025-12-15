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
#include <random>

using namespace mps2D;

// 计算标量场的梯度（使用LSMPS corrective matrix方法）
// 参数：
//   particle_idx: 粒子索引
//   scalar_field: 标量场值（如压力）
//   fluid_particles: 流体粒子对象
//   solid_particles: 固体粒子对象
//   corrective_matrix: corrective matrix (5x5)
//   smoothing_radius: 平滑半径
//   rho: 流体密度
//   g: 重力加速度大小
// 返回：梯度向量 (grad_x, grad_y)
double2 ComputeGradient(
    int particle_idx,
    const std::vector<double>& scalar_field,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const Eigen::Matrix<double, 5, 5>& corrective_matrix,
    double smoothing_radius,
    double rho,
    double g) {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  double phi_i = scalar_field[particle_idx];
  
  // 提取corrective matrix的前两行（C1和C2）
  Eigen::Matrix<double, 1, 5> C1 = corrective_matrix.row(0);  // x方向
  Eigen::Matrix<double, 1, 5> C2 = corrective_matrix.row(1);    // y方向
  
  double grad_x = 0.0;
  double grad_y = 0.0;
  
  // 处理流体邻域粒子
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    
    if (dist < 1e-10 || dist > smoothing_radius) continue;  // 避免除零
    
    double phi_j = scalar_field[j];
    double d_ij = (phi_j - phi_i) / dist;  // d_ij = (phi_j - phi_i) / r_ij
    double weight = WeightFunction(dist, smoothing_radius);
    
    // 计算基函数（与 manual 测试保持一致）
    Eigen::Matrix<double, 5, 1> P;
    P << dx / dist,                                    // x/r
         dy / dist,                                    // y/r
         dx * dx / (dist * smoothing_radius),          // x^2/(r*r_e)
         dy * dy / (dist * smoothing_radius),          // y^2/(r*r_e)
         dx * dy / (dist * smoothing_radius);          // x*y/(r*r_e)
    
    // 计算梯度贡献
    grad_x += weight * d_ij * (C1 * P)(0, 0);
    grad_y += weight * d_ij * (C2 * P)(0, 0);
  }
  
  // 处理固体邻域粒子（壁面粒子）
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    const double2& pos_j = solid_particles.position[j];
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    
    if (dist < 1e-10 || dist > smoothing_radius) continue;  // 避免除零
    
    // 获取壁面法向量并归一化
    const double2& normal = solid_particles.normal_vector[j];
    double nn = std::sqrt(normal.x * normal.x + normal.y * normal.y);
    double n_y = (nn > 1e-10) ? normal.y / nn : 0.0;
    
    double d_ij = -rho * g * n_y;
    double weight = WeightFunction(dist, smoothing_radius);
    
    // 壁面基函数
    Eigen::Matrix<double, 5, 1> P;
    P << normal.x,                                    // n_x
         normal.y,                                    // n_y
         2.0 * normal.x * dx / smoothing_radius,      // 2*n_x*x/r_e
         2.0 * normal.y * dy / smoothing_radius,      // 2*n_y*y/r_e
         (normal.y * dx + normal.x * dy) / smoothing_radius;  // (n_x*x + n_y*y)/r_e
    
    // 计算梯度贡献
    grad_x += weight * d_ij * (C1 * P)(0, 0);
    grad_y += weight * d_ij * (C2 * P)(0, 0);
  }
  
  return {grad_x, grad_y};
}

// 计算标量场的拉普拉斯算子（使用LSMPS corrective matrix方法）
// 参数：
//   particle_idx: 粒子索引
//   scalar_field: 标量场值（如压力）
//   fluid_particles: 流体粒子对象
//   solid_particles: 固体粒子对象
//   corrective_matrix: corrective matrix (5x5)
//   smoothing_radius: 平滑半径
//   rho: 流体密度
//   g: 重力加速度大小
// 返回：拉普拉斯算子
double ComputePressureLaplacian(
  int particle_idx,
  const std::vector<double>& scalar_field,
  const FluidParticle& fluid_particles,
  const SolidParticle& solid_particles,
  const Eigen::Matrix<double, 5, 5>& corrective_matrix,
  double smoothing_radius,
  double rho,
  double g 
  ) {

const double2& pos_i = fluid_particles.position[particle_idx];
double phi_i = scalar_field[particle_idx];

// 提取corrective matrix的前两行（C1和C2）
Eigen::Matrix<double, 1, 5> C3 = corrective_matrix.row(2);  // x方向
Eigen::Matrix<double, 1, 5> C4 = corrective_matrix.row(3);  // y方向

double pressure_laplacian = 0.0;

// 处理流体邻域粒子
for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
  const double2& pos_j = fluid_particles.position[j];
  double dx = pos_j.x - pos_i.x;
  double dy = pos_j.y - pos_i.y;
  double dist = ComputeDistance(pos_i, pos_j);
  
  if (dist < 1e-10 || dist > smoothing_radius) continue;  // 避免除零
  
  double phi_j = scalar_field[j];
  double d_ij = (phi_j - phi_i) / dist;  // d_ij = (phi_j - phi_i) / r_ij
  double weight = WeightFunction(dist, smoothing_radius);
  
  // 计算基函数（与 manual 测试保持一致）
  Eigen::Matrix<double, 5, 1> P;
  P << dx / dist,                                    // x/r
       dy / dist,                                    // y/r
       dx * dx / (dist * smoothing_radius),          // x^2/(r*r_e)
       dy * dy / (dist * smoothing_radius),          // y^2/(r*r_e)
       dx * dy / (dist * smoothing_radius);          // x*y/(r*r_e)
  
  // 计算拉普拉斯算子贡献 
  pressure_laplacian += weight * d_ij * (C3 * P)(0, 0);
  pressure_laplacian += weight * d_ij * (C4 * P)(0, 0);
}

// 处理固体邻域粒子（壁面粒子）
for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
  const double2& pos_j = solid_particles.position[j];
  double dx = pos_j.x - pos_i.x;
  double dy = pos_j.y - pos_i.y;
  double dist = ComputeDistance(pos_i, pos_j);
  
  if (dist < 1e-10 || dist > smoothing_radius) continue;  // 避免除零
  
  // 获取壁面法向量并归一化
  const double2& normal = solid_particles.normal_vector[j];
  double nn = std::sqrt(normal.x * normal.x + normal.y * normal.y);
  double n_y = (nn > 1e-10) ? normal.y / nn : 0.0;
  
  double d_ij = -rho * g * n_y;
  double weight = WeightFunction(dist, smoothing_radius);
  
  // 壁面基函数
  Eigen::Matrix<double, 5, 1> P;
  P << normal.x,                                    // n_x
       normal.y,                                    // n_y
       2.0 * normal.x * dx / smoothing_radius,      // 2*n_x*x/r_e
       2.0 * normal.y * dy / smoothing_radius,      // 2*n_y*y/r_e
       (normal.y * dx + normal.x * dy) / smoothing_radius;  // (n_x*x + n_y*y)/r_e
  
  // 计算梯度贡献
  pressure_laplacian += weight * d_ij * (C3 * P)(0, 0);
  pressure_laplacian += weight * d_ij * (C4 * P)(0, 0);
}

return pressure_laplacian;
}

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
  // 壁面粒子与流体粒子的初始距离应该等于粒子间距
  for (int i = 0; i < nx_bottom; ++i) {
    double x = i * particle_spacing;
    double y = -particle_spacing;  // 距离底部流体粒子一个粒子间距
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {0.0, 1.0};  // 向上
    ++idx;
  }
  
  // 左侧壁面（法向量向右）
  // 壁面粒子与流体粒子的初始距离应该等于粒子间距
  for (int j = 0; j < ny_left; ++j) {
    double x = -particle_spacing;  // 距离左侧流体粒子一个粒子间距
    double y = j * particle_spacing;
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {1.0, 0.0};  // 向右
    ++idx;
  }
  
  // 右侧壁面（法向量向左）
  // 壁面粒子与流体粒子的初始距离应该等于粒子间距
  for (int j = 0; j < ny_right; ++j) {
    double x = container_width;  // 距离右侧流体粒子一个粒子间距（容器宽度处是最后一个流体粒子）
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

// 输出压力数据到VTK文件（用于调试）
// 参数：
//   filename: 输出文件名
//   fluid_particles: 流体粒子
//   theoretical_pressure: 理论压力
//   computed_gradient: 计算出的梯度
//   theoretical_gradient: 理论梯度
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
  file << "Hydrostatic Pressure Test - MPS 2D\n";
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
  
  // 写入计算出的压力（如果有）
  if (fluid_particles.pressure.size() >= static_cast<size_t>(num_particles)) {
    file << "SCALARS computed_pressure float\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < num_particles; ++i) {
      file << fluid_particles.pressure[i] << "\n";
    }
  }
  
  // 写入计算出的梯度（压力梯度）
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
  
  // 写入压力梯度大小（计算值）
  file << "SCALARS computed_pressure_gradient_magnitude float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    double grad_mag = std::sqrt(computed_gradient[i].x * computed_gradient[i].x + 
                                 computed_gradient[i].y * computed_gradient[i].y);
    file << grad_mag << "\n";
  }
  
  // 写入压力梯度大小（理论值）
  file << "SCALARS theoretical_pressure_gradient_magnitude float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    double grad_mag = std::sqrt(theoretical_gradient[i].x * theoretical_gradient[i].x + 
                                 theoretical_gradient[i].y * theoretical_gradient[i].y);
    file << grad_mag << "\n";
  }
  
  // 写入压力梯度X分量（计算值）
  file << "SCALARS computed_pressure_gradient_x float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    file << computed_gradient[i].x << "\n";
  }
  
  // 写入压力梯度Y分量（计算值）
  file << "SCALARS computed_pressure_gradient_y float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    file << computed_gradient[i].y << "\n";
  }
  
  // 写入压力梯度X分量（理论值）
  file << "SCALARS theoretical_pressure_gradient_x float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    file << theoretical_gradient[i].x << "\n";
  }

  // 写入压力拉普拉斯算子
  file << "SCALARS pressure_laplacian float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    file << computed_laplacian[i] << "\n";
  }

  // 写入表面类型
  if (fluid_particles.surface_type.size() >= static_cast<size_t>(num_particles)) {
    file << "SCALARS surface_type int\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < num_particles; ++i) {
      file << static_cast<int>(fluid_particles.surface_type[i]) << "\n";
    }
  }
  
  // 写入Y坐标（用于检查压力分布）
  file << "SCALARS y_coordinate float\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; ++i) {
    file << fluid_particles.position[i].y << "\n";
  }
  
  file.close();
  std::cout << "  已输出VTK文件: " << filename << std::endl;
}

// 输出壁面粒子数据到VTK文件
// 参数：
//   filename: 输出文件名
//   solid_particles: 固体粒子（壁面粒子）
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
  file << "Wall Particles - MPS 2D\n";
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
  std::cout << "=== 静水压强测例（LSMPS梯度计算） ===" << std::endl;
  
  // 物理参数
  const double rho = 1000.0;  // 水的密度 (kg/m³)
  const double g = 9.8;       // 重力加速度 (m/s²)
  
  // 几何参数
  const double container_width = 1.0;   // 容器宽度 (m)
  const double container_height = 1.0;  // 容器高度 (m)
  const double water_height = 0.5;      // 水位高度 (m)
  
  // 粒子参数（目标约10000个粒子）
  // 总粒子数 ≈ (width/spacing) * (height/spacing)
  // 10000 ≈ (1.0/spacing) * (0.5/spacing) = 0.5/spacing²
  // spacing² ≈ 0.5/10000 = 0.00005
  // spacing ≈ 0.007
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
  
  // 判定自由面
  std::cout << "判定自由面..." << std::endl;
  SurfaceDetector surface_detector;
  surface_detector.DetectSurfaceParticles(
      fluid_particles, solid_particles,
      smoothing_radius, particle_spacing);
  
  // 统计自由面粒子数
  int num_surface = 0;
  int num_near_surface = 0;
  int num_inner = 0;
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    if (fluid_particles.surface_type[i] == SurfaceType::SURFACE) {
      ++num_surface;
    } else if (fluid_particles.surface_type[i] == SurfaceType::NEAR_SURFACE) {
      ++num_near_surface;
    } else if (fluid_particles.surface_type[i] == SurfaceType::INNER) {
      ++num_inner;
    }
  }
  std::cout << "  自由面粒子: " << num_surface << std::endl;
  std::cout << "  近自由面粒子: " << num_near_surface << std::endl;
  std::cout << "  内部粒子: " << num_inner << std::endl;
  
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
  
  // 计算corrective matrix
  std::cout << "\n计算corrective matrix..." << std::endl;
  CorrectiveMatrix corrective_matrix_calculator;
  std::vector<Eigen::Matrix<double, 5, 5>> matrices(fluid_particles.particle_num);
  
  int num_valid_matrix = 0;
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    matrices[i] = corrective_matrix_calculator.ComputeCorrectiveMatrix(
        i, fluid_particles, solid_particles, smoothing_radius, true);
    
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
  
  // 计算压力梯度
  std::cout << "\n计算压力梯度..." << std::endl;
  std::vector<double2> computed_gradient(fluid_particles.particle_num);
  std::vector<double2> theoretical_gradient(fluid_particles.particle_num);
  std::vector<double> computed_laplacian(fluid_particles.particle_num);
  std::vector<double> theoretical_laplacian(fluid_particles.particle_num);
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    // 计算梯度
    computed_gradient[i] = ComputeGradient(
        i, theoretical_pressure, fluid_particles, solid_particles,
        matrices[i], smoothing_radius, rho, g);
    
    // 理论梯度（静水压强）
    // 压力 p = rho * g * h，其中 h = water_height - y
    // dp/dy = rho * g * d(water_height - y)/dy = rho * g * (-1) = -rho * g
    // 所以y方向梯度应该是负数（压力随y增加而减小）
    theoretical_gradient[i] = {0.0, -rho * g};
    computed_laplacian[i] = ComputePressureLaplacian(
        i, theoretical_pressure, fluid_particles, solid_particles,
        matrices[i], smoothing_radius, rho, g);
    theoretical_laplacian[i] = 0;
  }
  
  // 统计误差
  std::cout << "\n梯度误差统计:" << std::endl;
  double total_error_x = 0.0;
  double total_error_y = 0.0;
  double max_error_x = 0.0;
  double max_error_y = 0.0;
  int num_valid_gradient = 0;
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    // 只统计内部粒子的误差（排除自由面附近）
    if (fluid_particles.surface_type[i] == SurfaceType::INNER) {
      double error_x = std::abs(computed_gradient[i].x - theoretical_gradient[i].x);
      double error_y = std::abs(computed_gradient[i].y - theoretical_gradient[i].y);
      
      total_error_x += error_x;
      total_error_y += error_y;
      
      if (error_x > max_error_x) max_error_x = error_x;
      if (error_y > max_error_y) max_error_y = error_y;
      
      ++num_valid_gradient;
    }
  }
  
  double avg_error_x = total_error_x / num_valid_gradient;
  double avg_error_y = total_error_y / num_valid_gradient;
  double relative_error_y = avg_error_y / (rho * g) * 100.0;
  
  std::cout << "  统计粒子数（内部粒子）: " << num_valid_gradient << std::endl;
  std::cout << "  X方向平均误差: " << avg_error_x << " Pa/m" << std::endl;
  std::cout << "  Y方向平均误差: " << avg_error_y << " Pa/m" << std::endl;
  std::cout << "  Y方向相对误差: " << std::fixed << std::setprecision(2) 
            << relative_error_y << "%" << std::endl;
  std::cout << "  X方向最大误差: " << max_error_x << " Pa/m" << std::endl;
  std::cout << "  Y方向最大误差: " << max_error_y << " Pa/m" << std::endl;
  std::cout << "  理论Y方向梯度: " << (-rho * g) << " Pa/m (压力随y增加而减小)" << std::endl;
  
  // 显示前几个粒子的详细信息（用于调试）
  std::cout << "\n前5个内部粒子的详细信息（用于调试）:" << std::endl;
  std::cout << std::setw(6) << "粒子" 
            << std::setw(12) << "位置Y" 
            << std::setw(12) << "压力"
            << std::setw(12) << "理论GradY"
            << std::setw(12) << "计算GradY"
            << std::setw(12) << "误差Y"
            << std::setw(8) << "邻域数" << std::endl;
  
  int count = 0;
  for (int i = 0; i < fluid_particles.particle_num && count < 5; ++i) {
    if (fluid_particles.surface_type[i] == SurfaceType::INNER) {
      double error_y = computed_gradient[i].y - theoretical_gradient[i].y;
      int num_neighbors = fluid_particles.fluid_neighbour_list[i].size() + 
                         fluid_particles.solid_neighbour_list[i].size();
      
      std::cout << std::setw(6) << i
                << std::setw(12) << std::fixed << std::setprecision(4) 
                << fluid_particles.position[i].y
                << std::setw(12) << theoretical_pressure[i]
                << std::setw(12) << theoretical_gradient[i].y
                << std::setw(12) << computed_gradient[i].y
                << std::setw(12) << error_y
                << std::setw(8) << num_neighbors << std::endl;
      ++count;
    }
  }
  
  // 显示前10个内部粒子的梯度对比
  std::cout << "\n前10个内部粒子的梯度对比:" << std::endl;
  std::cout << std::setw(6) << "粒子" 
            << std::setw(12) << "理论GradX" 
            << std::setw(12) << "计算GradX"
            << std::setw(12) << "误差X"
            << std::setw(12) << "理论GradY"
            << std::setw(12) << "计算GradY"
            << std::setw(12) << "误差Y" << std::endl;
  
  count = 0;
  for (int i = 0; i < fluid_particles.particle_num && count < 10; ++i) {
    if (fluid_particles.surface_type[i] == SurfaceType::INNER) {
      double error_x = computed_gradient[i].x - theoretical_gradient[i].x;
      double error_y = computed_gradient[i].y - theoretical_gradient[i].y;
      
      std::cout << std::setw(6) << i
                << std::setw(12) << std::fixed << std::setprecision(4) << theoretical_gradient[i].x
                << std::setw(12) << computed_gradient[i].x
                << std::setw(12) << error_x
                << std::setw(12) << theoretical_gradient[i].y
                << std::setw(12) << computed_gradient[i].y
                << std::setw(12) << error_y << std::endl;
      ++count;
    }
  }
  
  // 输出压力数据到VTK文件（用于调试）
  std::cout << "\n输出VTK文件..." << std::endl;
  WritePressureToVTK(
      "hydrostatic_pressure.vtk",
      fluid_particles,
      theoretical_pressure,
      computed_gradient,
      theoretical_gradient,
      computed_laplacian);
  
  // 输出壁面粒子到VTK文件
  WriteWallParticlesToVTK(
      "wall_particles.vtk",
      solid_particles);
  
  std::cout << "\n测试完成!" << std::endl;
  
  return 0;
}

