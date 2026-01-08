// 生成静水问题算例的工具程序
// 用途：生成指定粒子数的静水问题初始粒子分布文件

#include "../src/core/Particle.hpp"
#include "../src/core/FileOperator.hpp"
#include "../include/core/Types.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace mps2D;

// 生成静水压力测试场景
// 参数：
//   container_width: 容器宽度 (m)
//   container_height: 容器高度 (m)
//   water_height: 水位高度 (m)
//   particle_spacing: 粒子间距 (m)
//   fluid_particles: 输出流体粒子
//   solid_particles: 输出固体粒子
void GenerateHydrostaticTest(
    double container_width,
    double container_height,
    double water_height,
    double particle_spacing,
    FluidParticle& fluid_particles,
    SolidParticle& solid_particles) {
  
  // 计算x方向的粒子数（从0开始）
  int nx_fluid = static_cast<int>(std::round(container_width / particle_spacing)) + 1;
  
  // 计算y方向的粒子数
  int ny_fluid = static_cast<int>(std::round(water_height / particle_spacing)) + 1;
  
  int num_fluid = nx_fluid * ny_fluid;
  
  // 生成流体粒子（水）
  fluid_particles.particle_num = num_fluid;
  fluid_particles.position.resize(num_fluid);
  fluid_particles.velocity.resize(num_fluid);
  
  int idx = 0;
  for (int j = 0; j < ny_fluid; ++j) {
    for (int i = 0; i < nx_fluid; ++i) {
      double x = i * particle_spacing;
      double y = j * particle_spacing;
      fluid_particles.position[idx] = {x, y};
      fluid_particles.velocity[idx] = {0.0, 0.0};  // 静止流体
      ++idx;
    }
  }
  
  // 生成固体粒子（容器壁面）
  // 底部壁面
  int nx_bottom = static_cast<int>(std::round(container_width / particle_spacing)) + 1;
  
  // 左侧壁面
  int ny_left = static_cast<int>(std::round(container_height / particle_spacing)) + 1;
  
  // 右侧壁面
  int ny_right = ny_left;
  
  int num_solid = nx_bottom + 2 * ny_left;
  solid_particles.particle_num = num_solid;
  solid_particles.position.resize(num_solid);
  solid_particles.velocity.resize(num_solid);
  solid_particles.normal_vector.resize(num_solid);
  
  idx = 0;
  
  // 底部壁面（法向量向上）
  // 壁面粒子与流体粒子的初始距离等于粒子间距
  for (int i = 0; i < nx_bottom; ++i) {
    double x = i * particle_spacing;
    double y = -particle_spacing;  // 距离底部流体粒子一个粒子间距
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {0.0, 1.0};  // 向上
    ++idx;
  }
  
  // 左侧壁面（法向量向右）
  double left_x = -particle_spacing;  // 距离左侧流体粒子一个粒子间距
  for (int j = 0; j < ny_left; ++j) {
    double x = left_x;
    double y = j * particle_spacing;
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {1.0, 0.0};  // 向右
    ++idx;
  }
  
  // 右侧壁面（法向量向左）
  double right_x = container_width + particle_spacing;  // 距离右侧流体粒子一个粒子间距
  for (int j = 0; j < ny_right; ++j) {
    double x = right_x;
    double y = j * particle_spacing;
    solid_particles.position[idx] = {x, y};
    solid_particles.velocity[idx] = {0.0, 0.0};
    solid_particles.normal_vector[idx] = {-1.0, 0.0};  // 向左
    ++idx;
  }
}

int main(int argc, char* argv[]) {
  std::cout << "=== 生成静水问题算例 ===" << std::endl;
  
  // 目标：生成约2万粒子的算例
  // 计算合适的参数
  double container_width = 2.0;   // 容器宽度 (m)
  double container_height = 2.0;  // 容器高度 (m)
  double water_height = 1.0;      // 水位高度 (m)
  double particle_spacing = 0.01; // 粒子间距 (m)，可以根据需要调整
  
  // 如果通过命令行参数指定粒子间距
  if (argc > 1) {
    particle_spacing = std::stod(argv[1]);
  }
  
  // 计算粒子数
  int nx_fluid = static_cast<int>(std::round(container_width / particle_spacing)) + 1;
  int ny_fluid = static_cast<int>(std::round(water_height / particle_spacing)) + 1;
  int num_fluid = nx_fluid * ny_fluid;
  
  int nx_bottom = static_cast<int>(std::round(container_width / particle_spacing)) + 1;
  int ny_wall = static_cast<int>(std::round(container_height / particle_spacing)) + 1;
  int num_solid = nx_bottom + 2 * ny_wall;
  int total_particles = num_fluid + num_solid;
  
  std::cout << "\n几何参数:" << std::endl;
  std::cout << "  容器宽度: " << container_width << " m" << std::endl;
  std::cout << "  容器高度: " << container_height << " m" << std::endl;
  std::cout << "  水位高度: " << water_height << " m" << std::endl;
  std::cout << "  粒子间距: " << particle_spacing << " m" << std::endl;
  
  std::cout << "\n粒子数估算:" << std::endl;
  std::cout << "  流体粒子数: " << num_fluid << std::endl;
  std::cout << "  固体粒子数: " << num_solid << std::endl;
  std::cout << "  总粒子数: " << total_particles << std::endl;
  
  // 生成粒子
  FluidParticle fluid_particles("fluid");
  SolidParticle solid_particles("solid");
  
  std::cout << "\n生成粒子..." << std::endl;
  GenerateHydrostaticTest(
      container_width, container_height, water_height,
      particle_spacing, fluid_particles, solid_particles);
  
  // 输出文件
  FileOperator file_operator;
  
  std::string fluid_file = "data/fluid_particles_2d.txt";
  std::string solid_file = "data/solid_particles_2d.txt";
  
  std::cout << "\n输出粒子文件..." << std::endl;
  
  // 写入流体粒子文件（只包含位置和速度）
  std::ofstream fluid_out(fluid_file);
  if (!fluid_out.is_open()) {
    std::cerr << "错误：无法打开文件 " << fluid_file << std::endl;
    return 1;
  }
  fluid_out << std::fixed << std::setprecision(15);
  fluid_out << "# Fluid Particle Data File\n";
  fluid_out << "# Format: position_x position_y velocity_x velocity_y\n";
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    fluid_out << fluid_particles.position[i].x << " "
              << fluid_particles.position[i].y << " "
              << fluid_particles.velocity[i].x << " "
              << fluid_particles.velocity[i].y << "\n";
  }
  fluid_out.close();
  std::cout << "  已输出流体粒子文件: " << fluid_file << " (" << fluid_particles.particle_num << " 个粒子)" << std::endl;
  
  // 写入固体粒子文件（包含位置、速度和法向向量）
  std::ofstream solid_out(solid_file);
  if (!solid_out.is_open()) {
    std::cerr << "错误：无法打开文件 " << solid_file << std::endl;
    return 1;
  }
  solid_out << std::fixed << std::setprecision(15);
  solid_out << "# Solid Particle Data File\n";
  solid_out << "# Format: position_x position_y velocity_x velocity_y normal_x normal_y\n";
  for (int i = 0; i < solid_particles.particle_num; ++i) {
    solid_out << solid_particles.position[i].x << " "
              << solid_particles.position[i].y << " "
              << solid_particles.velocity[i].x << " "
              << solid_particles.velocity[i].y << " "
              << solid_particles.normal_vector[i].x << " "
              << solid_particles.normal_vector[i].y << "\n";
  }
  solid_out.close();
  std::cout << "  已输出固体粒子文件: " << solid_file << " (" << solid_particles.particle_num << " 个粒子)" << std::endl;
  
  std::cout << "\n算例生成完成！" << std::endl;
  std::cout << "  总粒子数: " << total_particles << std::endl;
  
  return 0;
}

