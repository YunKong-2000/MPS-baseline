#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <cassert>
#include "../src/neighbour_list/NeighborListSearcher.hpp"
#include "../src/core/Particle.hpp"

using namespace mps;

// 辅助函数：打印分隔线
void printSeparator(const std::string& title) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "  " << title << "\n";
    std::cout << std::string(60, '=') << "\n";
}

// 辅助函数：计算两点距离
double computeDistance(const double3& pos1, const double3& pos2) {
    double dx = pos1[0] - pos2[0];
    double dy = pos1[1] - pos2[1];
    double dz = pos1[2] - pos2[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// 辅助函数：验证邻居列表的正确性
bool verifyNeighborList(const FluidParticle& fluid_particles,
                       const SolidParticle& solid_particles,
                       double r_e) {
    bool all_correct = true;
    
    for (int i = 0; i < fluid_particles.particle_num; ++i) {
        // 验证流体邻居列表
        for (int j : fluid_particles.fluid_neighbour_list[i]) {
            if (j < 0 || j >= fluid_particles.particle_num) {
                std::cout << "  错误: 粒子 " << i << " 的流体邻居索引 " << j 
                          << " 超出范围\n";
                all_correct = false;
                continue;
            }
            if (i == j) {
                std::cout << "  错误: 粒子 " << i << " 将自己加入邻居列表\n";
                all_correct = false;
                continue;
            }
            double dist = computeDistance(fluid_particles.position[i],
                                         fluid_particles.position[j]);
            if (dist >= r_e) {
                std::cout << "  错误: 粒子 " << i << " 和 " << j 
                          << " 的距离 " << dist << " 大于 r_e " << r_e << "\n";
                all_correct = false;
            }
        }
        
        // 验证固体邻居列表
        for (int j : fluid_particles.solid_neighbour_list[i]) {
            if (j < 0 || j >= solid_particles.particle_num) {
                std::cout << "  错误: 粒子 " << i << " 的固体邻居索引 " << j 
                          << " 超出范围\n";
                all_correct = false;
                continue;
            }
            double dist = computeDistance(fluid_particles.position[i],
                                         solid_particles.position[j]);
            if (dist >= r_e) {
                std::cout << "  错误: 流体粒子 " << i << " 和固体粒子 " << j 
                          << " 的距离 " << dist << " 大于 r_e " << r_e << "\n";
                all_correct = false;
            }
        }
    }
    
    return all_correct;
}

// 测试1：基本功能测试 - 简单网格场景
void testBasicFunctionality() {
    printSeparator("测试1: 基本功能测试");
    
    NeighborListSearcher searcher;
    
    // 创建简单的测试数据：3x3x3网格，每个网格一个粒子
    FluidParticle fluid("test_fluid");
    SolidParticle solid("test_solid");
    
    double spacing = 1.0;
    double r_e = 1.5;  // 搜索半径
    double r_cell = 1.0;  // 网格尺寸
    
    // 创建9个流体粒子，排列成3x3网格
    fluid.particle_num = 9;
    fluid.position.resize(9);
    fluid.velocity.resize(9);
    fluid.density.resize(9);
    fluid.pressure.resize(9);
    fluid.surface_type.resize(9);
    
    int idx = 0;
    for (int z = 0; z < 3; ++z) {
        for (int y = 0; y < 3; ++y) {
            for (int x = 0; x < 3; ++x) {
                if (idx < 9) {
                    fluid.position[idx] = {x * spacing, y * spacing, z * spacing};
                    fluid.velocity[idx] = {0.0, 0.0, 0.0};
                    fluid.density[idx] = 1000.0;
                    fluid.pressure[idx] = 0.0;
                    fluid.surface_type[idx] = SurfaceType::INNER;
                    ++idx;
                }
            }
        }
    }
    
    // 创建1个固体粒子
    solid.particle_num = 1;
    solid.position.resize(1);
    solid.velocity.resize(1);
    solid.normal_vector.resize(1);
    solid.position[0] = {0.5, 0.5, 0.5};
    solid.velocity[0] = {0.0, 0.0, 0.0};
    solid.normal_vector[0] = {1.0, 0.0, 0.0};
    
    std::cout << "创建了 " << fluid.particle_num << " 个流体粒子和 " 
              << solid.particle_num << " 个固体粒子\n";
    std::cout << "r_e = " << r_e << ", r_cell = " << r_cell << "\n";
    
    // 构建邻居列表
    searcher.BuildNeighborList(fluid, solid, spacing, r_e, r_cell);
    
    // 验证结果
    std::cout << "\n邻居列表统计:\n";
    int total_fluid_neighbors = 0;
    int total_solid_neighbors = 0;
    for (int i = 0; i < fluid.particle_num; ++i) {
        total_fluid_neighbors += fluid.fluid_neighbour_list[i].size();
        total_solid_neighbors += fluid.solid_neighbour_list[i].size();
    }
    std::cout << "  总流体邻居数: " << total_fluid_neighbors << "\n";
    std::cout << "  总固体邻居数: " << total_solid_neighbors << "\n";
    
    // 验证正确性
    if (verifyNeighborList(fluid, solid, r_e)) {
        std::cout << "✓ 邻居列表验证通过\n";
    } else {
        std::cout << "✗ 邻居列表验证失败\n";
    }
}

// 测试2：空粒子列表测试
void testEmptyParticles() {
    printSeparator("测试2: 空粒子列表测试");
    
    NeighborListSearcher searcher;
    
    FluidParticle fluid("empty_fluid");
    SolidParticle solid("empty_solid");
    
    fluid.particle_num = 0;
    solid.particle_num = 0;
    
    double spacing = 1.0;
    double r_e = 1.5;
    double r_cell = 1.0;
    
    std::cout << "测试空粒子列表...\n";
    
    try {
        searcher.BuildNeighborList(fluid, solid, spacing, r_e, r_cell);
        std::cout << "✓ 空粒子列表处理成功\n";
        std::cout << "  流体粒子数: " << fluid.particle_num << "\n";
        std::cout << "  固体粒子数: " << solid.particle_num << "\n";
    } catch (const std::exception& e) {
        std::cout << "✗ 空粒子列表处理失败: " << e.what() << "\n";
    }
}

// 测试3：只有流体粒子
void testFluidOnly() {
    printSeparator("测试3: 只有流体粒子");
    
    NeighborListSearcher searcher;
    
    FluidParticle fluid("fluid_only");
    SolidParticle solid("empty_solid");
    
    double spacing = 1.0;
    double r_e = 1.5;
    double r_cell = 1.0;
    
    // 创建5个流体粒子，排列成一条线
    fluid.particle_num = 5;
    fluid.position.resize(5);
    fluid.velocity.resize(5);
    fluid.density.resize(5);
    fluid.pressure.resize(5);
    fluid.surface_type.resize(5);
    
    for (int i = 0; i < 5; ++i) {
        fluid.position[i] = {i * spacing, 0.0, 0.0};
        fluid.velocity[i] = {0.0, 0.0, 0.0};
        fluid.density[i] = 1000.0;
        fluid.pressure[i] = 0.0;
        fluid.surface_type[i] = SurfaceType::INNER;
    }
    
    solid.particle_num = 0;
    
    std::cout << "创建了 " << fluid.particle_num << " 个流体粒子\n";
    
    searcher.BuildNeighborList(fluid, solid, spacing, r_e, r_cell);
    
    std::cout << "\n邻居列表统计:\n";
    for (int i = 0; i < fluid.particle_num; ++i) {
        std::cout << "  粒子 " << i << ": " 
                  << fluid.fluid_neighbour_list[i].size() << " 个流体邻居, "
                  << fluid.solid_neighbour_list[i].size() << " 个固体邻居\n";
    }
    
    if (verifyNeighborList(fluid, solid, r_e)) {
        std::cout << "✓ 邻居列表验证通过\n";
    } else {
        std::cout << "✗ 邻居列表验证失败\n";
    }
}

// 测试4：密集粒子场景
void testDenseParticles() {
    printSeparator("测试4: 密集粒子场景");
    
    NeighborListSearcher searcher;
    
    FluidParticle fluid("dense_fluid");
    SolidParticle solid("dense_solid");
    
    double spacing = 0.5;
    double r_e = 1.0;
    double r_cell = 0.6;
    
    // 创建27个流体粒子（3x3x3）
    fluid.particle_num = 27;
    fluid.position.resize(27);
    fluid.velocity.resize(27);
    fluid.density.resize(27);
    fluid.pressure.resize(27);
    fluid.surface_type.resize(27);
    
    int idx = 0;
    for (int z = 0; z < 3; ++z) {
        for (int y = 0; y < 3; ++y) {
            for (int x = 0; x < 3; ++x) {
                fluid.position[idx] = {x * spacing, y * spacing, z * spacing};
                fluid.velocity[idx] = {0.0, 0.0, 0.0};
                fluid.density[idx] = 1000.0;
                fluid.pressure[idx] = 0.0;
                fluid.surface_type[idx] = SurfaceType::INNER;
                ++idx;
            }
        }
    }
    
    // 创建8个固体粒子，位于角落
    solid.particle_num = 8;
    solid.position.resize(8);
    solid.velocity.resize(8);
    solid.normal_vector.resize(8);
    
    idx = 0;
    for (int z = 0; z < 2; ++z) {
        for (int y = 0; y < 2; ++y) {
            for (int x = 0; x < 2; ++x) {
                solid.position[idx] = {x * spacing - 0.2, y * spacing - 0.2, z * spacing - 0.2};
                solid.velocity[idx] = {0.0, 0.0, 0.0};
                solid.normal_vector[idx] = {1.0, 0.0, 0.0};
                ++idx;
            }
        }
    }
    
    std::cout << "创建了 " << fluid.particle_num << " 个流体粒子和 " 
              << solid.particle_num << " 个固体粒子\n";
    std::cout << "r_e = " << r_e << ", r_cell = " << r_cell << "\n";
    
    searcher.BuildNeighborList(fluid, solid, spacing, r_e, r_cell);
    
    std::cout << "\n邻居列表统计:\n";
    int total_fluid_neighbors = 0;
    int total_solid_neighbors = 0;
    int max_fluid_neighbors = 0;
    int max_solid_neighbors = 0;
    
    for (int i = 0; i < fluid.particle_num; ++i) {
        int fluid_count = fluid.fluid_neighbour_list[i].size();
        int solid_count = fluid.solid_neighbour_list[i].size();
        total_fluid_neighbors += fluid_count;
        total_solid_neighbors += solid_count;
        max_fluid_neighbors = std::max(max_fluid_neighbors, fluid_count);
        max_solid_neighbors = std::max(max_solid_neighbors, solid_count);
    }
    
    std::cout << "  总流体邻居数: " << total_fluid_neighbors << "\n";
    std::cout << "  总固体邻居数: " << total_solid_neighbors << "\n";
    std::cout << "  最大流体邻居数: " << max_fluid_neighbors << "\n";
    std::cout << "  最大固体邻居数: " << max_solid_neighbors << "\n";
    
    if (verifyNeighborList(fluid, solid, r_e)) {
        std::cout << "✓ 邻居列表验证通过\n";
    } else {
        std::cout << "✗ 邻居列表验证失败\n";
    }
}

// 测试5：粒子排序验证
void testParticleSorting() {
    printSeparator("测试5: 粒子排序验证");
    
    NeighborListSearcher searcher;
    
    FluidParticle fluid("sort_test");
    SolidParticle solid("sort_solid");
    
    double spacing = 1.0;
    double r_e = 1.5;
    double r_cell = 1.0;
    
    // 创建粒子，故意打乱顺序
    fluid.particle_num = 5;
    fluid.position.resize(5);
    fluid.velocity.resize(5);
    fluid.density.resize(5);
    fluid.pressure.resize(5);
    fluid.surface_type.resize(5);
    
    // 打乱顺序：z=2, z=0, z=1, z=2, z=0
    fluid.position[0] = {0.0, 0.0, 2.0};
    fluid.position[1] = {0.0, 0.0, 0.0};
    fluid.position[2] = {0.0, 0.0, 1.0};
    fluid.position[3] = {1.0, 0.0, 2.0};
    fluid.position[4] = {1.0, 0.0, 0.0};
    
    for (int i = 0; i < 5; ++i) {
        fluid.velocity[i] = {static_cast<double>(i), 0.0, 0.0};
        fluid.density[i] = 1000.0 + i;
        fluid.pressure[i] = static_cast<double>(i);
        fluid.surface_type[i] = SurfaceType::INNER;
    }
    
    solid.particle_num = 0;
    
    std::cout << "排序前粒子位置:\n";
    for (int i = 0; i < fluid.particle_num; ++i) {
        std::cout << "  粒子 " << i << ": (" 
                  << fluid.position[i][0] << ", "
                  << fluid.position[i][1] << ", "
                  << fluid.position[i][2] << ")\n";
    }
    
    // 保存原始数据用于验证
    std::vector<double3> original_positions = fluid.position;
    std::vector<double3> original_velocities = fluid.velocity;
    std::vector<double> original_densities = fluid.density;
    
    searcher.BuildNeighborList(fluid, solid, spacing, r_e, r_cell);
    
    std::cout << "\n排序后粒子位置:\n";
    for (int i = 0; i < fluid.particle_num; ++i) {
        std::cout << "  粒子 " << i << ": (" 
                  << fluid.position[i][0] << ", "
                  << fluid.position[i][1] << ", "
                  << fluid.position[i][2] << ")\n";
    }
    
    // 验证属性是否一起移动
    bool sorting_correct = true;
    for (int i = 0; i < fluid.particle_num; ++i) {
        // 找到原始位置对应的索引
        bool found = false;
        for (int j = 0; j < fluid.particle_num; ++j) {
            if (std::abs(fluid.position[i][0] - original_positions[j][0]) < 1e-9 &&
                std::abs(fluid.position[i][1] - original_positions[j][1]) < 1e-9 &&
                std::abs(fluid.position[i][2] - original_positions[j][2]) < 1e-9) {
                // 检查其他属性是否匹配
                if (std::abs(fluid.velocity[i][0] - original_velocities[j][0]) > 1e-9 ||
                    std::abs(fluid.density[i] - original_densities[j]) > 1e-9) {
                    std::cout << "  错误: 粒子 " << i << " 的属性不匹配\n";
                    sorting_correct = false;
                }
                found = true;
                break;
            }
        }
        if (!found) {
            std::cout << "  错误: 找不到粒子 " << i << " 的原始位置\n";
            sorting_correct = false;
        }
    }
    
    if (sorting_correct) {
        std::cout << "✓ 粒子排序验证通过（属性一起移动）\n";
    } else {
        std::cout << "✗ 粒子排序验证失败\n";
    }
}

// 测试6：边界情况 - 粒子在边界上
void testBoundaryParticles() {
    printSeparator("测试6: 边界粒子测试");
    
    NeighborListSearcher searcher;
    
    FluidParticle fluid("boundary_fluid");
    SolidParticle solid("boundary_solid");
    
    double spacing = 1.0;
    double r_e = 1.2;
    double r_cell = 1.0;
    
    // 创建粒子，一些在边界上
    fluid.particle_num = 4;
    fluid.position.resize(4);
    fluid.velocity.resize(4);
    fluid.density.resize(4);
    fluid.pressure.resize(4);
    fluid.surface_type.resize(4);
    
    fluid.position[0] = {0.0, 0.0, 0.0};  // 原点
    fluid.position[1] = {10.0, 0.0, 0.0};  // 远离
    fluid.position[2] = {0.1, 0.1, 0.1};  // 接近原点
    fluid.position[3] = {9.9, 0.0, 0.0};  // 接近第二个粒子
    
    for (int i = 0; i < 4; ++i) {
        fluid.velocity[i] = {0.0, 0.0, 0.0};
        fluid.density[i] = 1000.0;
        fluid.pressure[i] = 0.0;
        fluid.surface_type[i] = SurfaceType::INNER;
    }
    
    solid.particle_num = 0;
    
    std::cout << "创建了 " << fluid.particle_num << " 个边界粒子\n";
    
    searcher.BuildNeighborList(fluid, solid, spacing, r_e, r_cell);
    
    std::cout << "\n邻居列表:\n";
    for (int i = 0; i < fluid.particle_num; ++i) {
        std::cout << "  粒子 " << i << " (" 
                  << fluid.position[i][0] << ", "
                  << fluid.position[i][1] << ", "
                  << fluid.position[i][2] << "): "
                  << fluid.fluid_neighbour_list[i].size() << " 个邻居\n";
        for (int j : fluid.fluid_neighbour_list[i]) {
            double dist = computeDistance(fluid.position[i], fluid.position[j]);
            std::cout << "    -> 粒子 " << j << ", 距离: " << dist << "\n";
        }
    }
    
    if (verifyNeighborList(fluid, solid, r_e)) {
        std::cout << "✓ 边界粒子测试通过\n";
    } else {
        std::cout << "✗ 边界粒子测试失败\n";
    }
}

// 测试7：大网格尺寸测试
void testLargeCellSize() {
    printSeparator("测试7: 大网格尺寸测试");
    
    NeighborListSearcher searcher;
    
    FluidParticle fluid("large_cell");
    SolidParticle solid("large_cell_solid");
    
    double spacing = 1.0;
    double r_e = 1.5;
    double r_cell = 5.0;  // 大网格
    
    // 创建粒子
    fluid.particle_num = 5;
    fluid.position.resize(5);
    fluid.velocity.resize(5);
    fluid.density.resize(5);
    fluid.pressure.resize(5);
    fluid.surface_type.resize(5);
    
    for (int i = 0; i < 5; ++i) {
        fluid.position[i] = {i * spacing, 0.0, 0.0};
        fluid.velocity[i] = {0.0, 0.0, 0.0};
        fluid.density[i] = 1000.0;
        fluid.pressure[i] = 0.0;
        fluid.surface_type[i] = SurfaceType::INNER;
    }
    
    solid.particle_num = 0;
    
    std::cout << "创建了 " << fluid.particle_num << " 个粒子\n";
    std::cout << "r_cell = " << r_cell << " (大网格)\n";
    
    searcher.BuildNeighborList(fluid, solid, spacing, r_e, r_cell);
    
    std::cout << "\n邻居列表统计:\n";
    for (int i = 0; i < fluid.particle_num; ++i) {
        std::cout << "  粒子 " << i << ": " 
                  << fluid.fluid_neighbour_list[i].size() << " 个邻居\n";
    }
    
    if (verifyNeighborList(fluid, solid, r_e)) {
        std::cout << "✓ 大网格尺寸测试通过\n";
    } else {
        std::cout << "✗ 大网格尺寸测试失败\n";
    }
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "MPS Baseline - NeighborListSearcher 功能测试\n";
    std::cout << std::string(60, '=') << "\n";
    
    try {
        // 运行所有测试
        testBasicFunctionality();
        testEmptyParticles();
        testFluidOnly();
        testDenseParticles();
        testParticleSorting();
        testBoundaryParticles();
        testLargeCellSize();
        
        printSeparator("所有测试完成");
        std::cout << "\n测试总结:\n";
        std::cout << "  - 基本功能测试\n";
        std::cout << "  - 空粒子列表测试\n";
        std::cout << "  - 只有流体粒子测试\n";
        std::cout << "  - 密集粒子场景测试\n";
        std::cout << "  - 粒子排序验证\n";
        std::cout << "  - 边界粒子测试\n";
        std::cout << "  - 大网格尺寸测试\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\n错误: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}

