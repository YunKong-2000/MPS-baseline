#include "../src/PPE/PPEMatrixBuilder.hpp"
#include "../src/lsmps/CorrectiveMatrix.hpp"
#include "../src/neighbour_list/NeighborListSearcher.hpp"
#include "../src/core/Particle.hpp"
#include "core/Types.h"
#include "core/MPSUtils.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace mps2D;

// 生成静水压力测试场景
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
      fluid_particles.velocity[idx] = {0.0, 0.0};  // 静止流体
      fluid_particles.pressure[idx] = 0.0;  // 稍后通过PPE求解
      fluid_particles.density[idx] = 1000.0;  // 水的密度
      fluid_particles.surface_type[idx] = SurfaceType::INNER;
      ++idx;
    }
  }
  
  // 生成固体粒子（容器壁面）- 两层壁面粒子
  // 底部壁面
  int nx_bottom = static_cast<int>(container_width / particle_spacing) + 1;
  // 左侧壁面
  int ny_left = static_cast<int>(container_height / particle_spacing) + 1;
  // 右侧壁面
  int ny_right = ny_left;
  
  // 两层壁面：每层都有底部、左侧、右侧
  int num_solid = 2 * (nx_bottom + ny_left + ny_right);
  solid_particles.particle_num = num_solid;
  solid_particles.position.resize(num_solid);
  solid_particles.velocity.resize(num_solid);
  solid_particles.normal_vector.resize(num_solid);
  
  idx = 0;
  
  // 底部壁面 - 第一层（法向量向上）
  for (int i = 0; i < nx_bottom; ++i) {
    double x = i * particle_spacing;
    double y = -particle_spacing;  // 距离底部流体粒子一个粒子间距
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {0.0, 1.0};  // 向上
    ++idx;
  }
  
  // 底部壁面 - 第二层（法向量向上）
  for (int i = 0; i < nx_bottom; ++i) {
    double x = i * particle_spacing;
    double y = -2.0 * particle_spacing;  // 第二层，再向外一个粒子间距
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {0.0, 1.0};  // 向上
    ++idx;
  }
  
  // 左侧壁面 - 第一层（法向量向右）
  for (int j = 0; j < ny_left; ++j) {
    double x = -particle_spacing;
    double y = j * particle_spacing;
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {1.0, 0.0};  // 向右
    ++idx;
  }
  
  // 左侧壁面 - 第二层（法向量向右）
  for (int j = 0; j < ny_left; ++j) {
    double x = -2.0 * particle_spacing;  // 第二层，再向外一个粒子间距
    double y = j * particle_spacing;
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {1.0, 0.0};  // 向右
    ++idx;
  }
  
  // 右侧壁面 - 第一层（法向量向左）
  for (int j = 0; j < ny_right; ++j) {
    double x = container_width;
    double y = j * particle_spacing;
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {-1.0, 0.0};  // 向左
    ++idx;
  }
  
  // 右侧壁面 - 第二层（法向量向左）
  for (int j = 0; j < ny_right; ++j) {
    double x = container_width + particle_spacing;  // 第二层，再向外一个粒子间距
    double y = j * particle_spacing;
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {-1.0, 0.0};  // 向左
    ++idx;
  }
}

int main() {
  std::cout << "=== PPE求解模块静水压力测试 ===" << std::endl;
  
  // 物理参数
  const double rho = 1000.0;  // 水的密度 (kg/m³)
  const double g = 9.8;       // 重力加速度 (m/s²)
  const double gravity_x = 0.0;
  const double gravity_y = -g;
  
  // 几何参数
  const double container_width = 1.0;   // 容器宽度 (m)
  const double container_height = 1.0;  // 容器高度 (m)
  const double water_height = 0.5;      // 水位高度 (m)
  
  // 粒子参数（目标约3000-5000个粒子）
  // 总粒子数 ≈ (width/spacing) * (height/spacing)
  // 3000 ≈ (1.0/spacing) * (0.5/spacing) = 0.5/spacing²
  // spacing² ≈ 0.5/3000 = 0.000167
  // spacing ≈ 0.013
  const double particle_spacing = 0.013;  // 粒子间距 (m)
  const double particle_radius = particle_spacing / 2.0;
  const double smoothing_radius = 2.1 * particle_spacing;
  const double cell_size = 2.0 * smoothing_radius;
  
  // 时间步长
  const double time_step = 0.001;  // 时间步长 (s)
  
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
  std::cout << "  网格尺寸: " << cell_size << " m" << std::endl;
  
  // 生成测试场景
  FluidParticle fluid_particles("fluid");
  SolidParticle solid_particles("solid");
  
  std::cout << "\n生成测试场景..." << std::endl;
  GenerateHydrostaticTest(
      container_width, container_height, water_height,
      particle_spacing, fluid_particles, solid_particles);
  
  std::cout << "  流体粒子数: " << fluid_particles.particle_num << std::endl;
  std::cout << "  固体粒子数: " << solid_particles.particle_num << std::endl;
  std::cout << "  总粒子数: " << fluid_particles.particle_num + solid_particles.particle_num << std::endl;
  
  // 构建邻域列表
  std::cout << "\n构建邻域列表..." << std::endl;
  NeighborListSearcher neighbor_searcher;
  neighbor_searcher.BuildNeighborList(
      fluid_particles, solid_particles,
      particle_radius, smoothing_radius, cell_size);
  
  // 统计邻域信息
  int total_fluid_neighbors = 0;
  int total_solid_neighbors = 0;
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    total_fluid_neighbors += fluid_particles.fluid_neighbour_list[i].size();
    total_solid_neighbors += fluid_particles.solid_neighbour_list[i].size();
  }
  std::cout << "  平均流体邻域粒子数: " 
            << static_cast<double>(total_fluid_neighbors) / fluid_particles.particle_num << std::endl;
  std::cout << "  平均固体邻域粒子数: " 
            << static_cast<double>(total_solid_neighbors) / fluid_particles.particle_num << std::endl;
  
  // 计算corrective matrix
  std::cout << "\n计算corrective matrix..." << std::endl;
  CorrectiveMatrix corrective_matrix_calc;
  std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>> 
      corrective_matrices(fluid_particles.particle_num);
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    corrective_matrices[i] = corrective_matrix_calc.ComputeCorrectiveMatrix(
        i, fluid_particles, solid_particles, smoothing_radius, false);
  }
  std::cout << "  完成" << std::endl;
  
  // 构建PPE系数矩阵和右边项（带调试信息）
  std::cout << "\n构建PPE系数矩阵和右边项..." << std::endl;
  PPEMatrixBuilder matrix_builder;
  Eigen::SparseMatrix<double, Eigen::ColMajor> A;
  Eigen::VectorXd b;
  PPEMatrixBuilder::DebugInfo debug_info;
  
  bool success = matrix_builder.BuildPPEMatrixWithDebug(
      fluid_particles, solid_particles, corrective_matrices,
      smoothing_radius, rho, time_step, gravity_x, gravity_y, A, b, debug_info);
  
  if (!success) {
    std::cerr << "错误：构建PPE矩阵失败" << std::endl;
    return 1;
  }
  
  std::cout << "  矩阵大小: " << A.rows() << " x " << A.cols() << std::endl;
  std::cout << "  非零元素数: " << A.nonZeros() << std::endl;
  std::cout << "  稀疏度: " << (1.0 - static_cast<double>(A.nonZeros()) / (A.rows() * A.cols())) * 100.0 
            << "%" << std::endl;
  
  std::cout << "\nPPE系数矩阵构建完成！" << std::endl;
  
  return 0;
}
