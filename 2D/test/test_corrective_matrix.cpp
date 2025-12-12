#include "../src/lsmps/CorrectiveMatrix.hpp"
#include "../src/core/Particle.hpp"
#include "core/Types.h"
#include "core/MPSUtils.h"
#include "../src/neighbour_list/NeighborListSearcher.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <algorithm>

using namespace mps2D;

// 测试corrective matrix计算
int main() {
  std::cout << "=== LSMPS Corrective Matrix 测试程序 ===" << std::endl;
  
  // 创建测试粒子
  FluidParticle fluid_particles("fluid");
  SolidParticle solid_particles("solid");
  
  // 读取粒子数据
  std::string fluid_file = "data/fluid_particles_2d.txt";
  std::string solid_file = "data/solid_particles_2d.txt";
  
  int num_fluid = fluid_particles.getParticleFromFile(fluid_file);
  int num_solid = solid_particles.getParticleFromFile(solid_file);
  
  std::cout << "读取流体粒子数: " << num_fluid << std::endl;
  std::cout << "读取固体粒子数: " << num_solid << std::endl;
  
  if (num_fluid <= 0) {
    std::cerr << "错误: 未读取到流体粒子数据 (返回: " << num_fluid << ")" << std::endl;
    std::cerr << "尝试的文件路径: " << fluid_file << std::endl;
    return 1;
  }
  
  // 设置仿真参数
  // 注意：粒子间距约为0.1（从数据文件看），平滑半径应该至少是粒子间距的2-3倍
  // 才能保证每个粒子有足够的邻域粒子（至少5个）来计算corrective matrix
  double particle_spacing = 0.1;  // 粒子间距
  double particle_radius = particle_spacing / 2.0;  // 粒子半径通常为间距的一半
  double smoothing_radius = 2.5 * particle_spacing;  // 平滑半径设为间距的2.5倍，确保覆盖多个粒子
  double cell_size = 2.0 * smoothing_radius;  // 网格尺寸设为平滑半径的2倍
  
  std::cout << "\n仿真参数:" << std::endl;
  std::cout << "  粒子间距: " << particle_spacing << std::endl;
  std::cout << "  粒子半径: " << particle_radius << std::endl;
  std::cout << "  平滑半径: " << smoothing_radius << std::endl;
  std::cout << "  网格尺寸: " << cell_size << std::endl;
  
  // 构建邻居列表
  NeighborListSearcher neighbor_searcher;
  neighbor_searcher.BuildNeighborList(
      fluid_particles, solid_particles,
      particle_radius, smoothing_radius, cell_size);
  
  std::cout << "\n开始计算corrective matrix..." << std::endl;
  std::cout << "平滑半径: " << smoothing_radius << std::endl;
  
  // 检查前几个粒子的邻域情况
  std::cout << "\n前5个粒子的邻域信息:" << std::endl;
  for (int i = 0; i < std::min(5, num_fluid); ++i) {
    int num_fluid_neighbors = fluid_particles.fluid_neighbour_list[i].size();
    int num_solid_neighbors = fluid_particles.solid_neighbour_list[i].size();
    int total_neighbors = num_fluid_neighbors + num_solid_neighbors;
    std::cout << "  粒子 " << i << ": 流体邻域=" << num_fluid_neighbors 
              << ", 固体邻域=" << num_solid_neighbors 
              << ", 总计=" << total_neighbors << std::endl;
  }
  
  // 创建corrective matrix计算器
  CorrectiveMatrix corrective_matrix_calculator;
  
  // 计算所有粒子的corrective matrix（逐个计算）
  std::vector<Eigen::Matrix<double, 5, 5>> matrices(num_fluid);
  
  // 统计失败原因
  int num_insufficient_neighbors = 0;
  int num_singular_matrix = 0;
  int num_success = 0;
  
  auto start_time = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < num_fluid; ++i) {
    int num_fluid_neighbors = fluid_particles.fluid_neighbour_list[i].size();
    int num_solid_neighbors = fluid_particles.solid_neighbour_list[i].size();
    int total_neighbors = num_fluid_neighbors + num_solid_neighbors;
    
    if (total_neighbors < 5) {
      ++num_insufficient_neighbors;
    }
    
    matrices[i] = corrective_matrix_calculator.ComputeCorrectiveMatrix(
        i, fluid_particles, solid_particles, smoothing_radius);
    
    // 检查是否是单位矩阵（表示失败）
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
    
    if (is_identity && total_neighbors >= 5) {
      ++num_singular_matrix;
    } else if (!is_identity) {
      ++num_success;
    }
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      end_time - start_time).count();
  
  std::cout << "\n计算完成，耗时: " << duration << " ms" << std::endl;
  std::cout << "平均每个粒子耗时: " << std::fixed << std::setprecision(4) 
            << (static_cast<double>(duration) / num_fluid) << " ms" << std::endl;
  
  std::cout << "\n失败原因统计:" << std::endl;
  std::cout << "  邻域不足(<5个): " << num_insufficient_neighbors << std::endl;
  std::cout << "  矩阵奇异(>=5个邻域但不可逆): " << num_singular_matrix << std::endl;
  std::cout << "  计算成功: " << num_success << std::endl;
  
  const double tolerance = 1e-6;
  
  // 显示前5个粒子的corrective matrix
  std::cout << "\n前5个粒子的corrective matrix:" << std::endl;
  for (int i = 0; i < std::min(5, num_fluid); ++i) {
    std::cout << "\n粒子 " << i << ":" << std::endl;
    std::cout << matrices[i] << std::endl;
    
    // 检查矩阵是否为奇异矩阵（单位矩阵表示计算失败）
    bool is_identity = true;
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
    
    if (is_identity) {
      std::cout << "  (注意: 此矩阵为单位矩阵，可能表示邻域不足或计算失败)" << std::endl;
    }
  }
  
  // 统计信息
  int num_valid = 0;
  int num_invalid = 0;
  
  for (const auto& matrix : matrices) {
    bool is_identity = true;
    for (int row = 0; row < 5; ++row) {
      for (int col = 0; col < 5; ++col) {
        double expected = (row == col) ? 1.0 : 0.0;
        if (std::abs(matrix(row, col) - expected) > tolerance) {
          is_identity = false;
          break;
        }
      }
      if (!is_identity) break;
    }
    
    if (is_identity) {
      ++num_invalid;
    } else {
      ++num_valid;
    }
  }
  
  std::cout << "\n统计信息:" << std::endl;
  std::cout << "  有效corrective matrix数量: " << num_valid << std::endl;
  std::cout << "  无效corrective matrix数量: " << num_invalid << std::endl;
  std::cout << "  有效率: " << std::fixed << std::setprecision(2) 
            << (100.0 * num_valid / num_fluid) << "%" << std::endl;
  
  std::cout << "\n测试完成!" << std::endl;
  
  return 0;
}

