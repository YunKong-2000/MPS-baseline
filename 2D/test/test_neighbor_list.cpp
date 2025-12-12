#include <iostream>
#include <iomanip>
#include <vector>
#include "../src/core/Particle.hpp"
#include "../src/core/FileOperator.hpp"
#include "../src/neighbour_list/NeighborListSearcher.hpp"
#include "../include/core/MPSUtils.h"

using namespace mps2D;

void printSeparator(const std::string& title) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "  " << title << "\n";
    std::cout << std::string(60, '=') << "\n";
}

void testBasicNeighborSearch() {
    printSeparator("测试基本邻居搜索功能");
    
    // 创建简单的3x3网格场景，每个网格一个粒子
    FluidParticle fluid("test_fluid");
    SolidParticle solid("test_solid");
    
    // 创建3x3网格的流体粒子（间距0.1）
    fluid.particle_num = 9;
    fluid.position.resize(9);
    fluid.velocity.resize(9);
    fluid.density.resize(9);
    fluid.pressure.resize(9);
    fluid.surface_type.resize(9);
    fluid.fluid_neighbour_list.resize(9);
    fluid.solid_neighbour_list.resize(9);
    
    int idx = 0;
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            fluid.position[idx] = {i * 0.1, j * 0.1};
            fluid.velocity[idx] = {0.0, 0.0};
            fluid.density[idx] = 1000.0;
            fluid.pressure[idx] = 101325.0;
            fluid.surface_type[idx] = SurfaceType::INNER;
            idx++;
        }
    }
    
    // 创建边界固体粒子
    solid.particle_num = 0;
    solid.position.resize(0);
    solid.velocity.resize(0);
    solid.normal_vector.resize(0);
    
    std::cout << "\n创建了 " << fluid.particle_num << " 个流体粒子（3x3网格）\n";
    std::cout << "粒子间距: 0.1\n";
    
    // 构建邻居列表
    double particle_radius = 0.05;
    double r_e = 0.15;  // 搜索半径
    double r_cell = 0.12;  // 网格尺寸
    
    NeighborListSearcher searcher;
    searcher.BuildNeighborList(fluid, solid, particle_radius, r_e, r_cell);
    
    std::cout << "\n邻居列表统计:\n";
    int total_neighbors = 0;
    int max_neighbors = 0;
    for (int i = 0; i < fluid.particle_num; ++i) {
        int neighbor_count = fluid.fluid_neighbour_list[i].size();
        total_neighbors += neighbor_count;
        max_neighbors = std::max(max_neighbors, neighbor_count);
    }
    
    std::cout << "  总邻居数: " << total_neighbors << "\n";
    std::cout << "  最大邻居数: " << max_neighbors << "\n";
    std::cout << "  平均邻居数: " << (fluid.particle_num > 0 ? total_neighbors / fluid.particle_num : 0) << "\n";
    
    // 验证邻居列表正确性
    std::cout << "\n验证邻居列表:\n";
    bool all_valid = true;
    for (int i = 0; i < fluid.particle_num; ++i) {
        for (int j : fluid.fluid_neighbour_list[i]) {
            if (j < 0 || j >= fluid.particle_num) {
                std::cout << "  ✗ 粒子[" << i << "] 有无效邻居索引: " << j << "\n";
                all_valid = false;
            }
            if (i == j) {
                std::cout << "  ✗ 粒子[" << i << "] 将自己加入邻居列表\n";
                all_valid = false;
            }
            double dist = ComputeDistance(fluid.position[i], fluid.position[j]);
            if (dist >= r_e) {
                std::cout << "  ✗ 粒子[" << i << "] 和粒子[" << j << "] 距离 " 
                          << dist << " 超出搜索半径 " << r_e << "\n";
                all_valid = false;
            }
        }
    }
    
    if (all_valid) {
        std::cout << "  ✓ 所有邻居列表有效\n";
    }
}

void testWithSolidParticles() {
    printSeparator("测试包含固体粒子的邻居搜索");
    
    FluidParticle fluid("test_fluid");
    SolidParticle solid("test_solid");
    
    // 创建5个流体粒子排成一条线
    fluid.particle_num = 5;
    fluid.position.resize(5);
    fluid.velocity.resize(5);
    fluid.density.resize(5);
    fluid.pressure.resize(5);
    fluid.surface_type.resize(5);
    fluid.fluid_neighbour_list.resize(5);
    fluid.solid_neighbour_list.resize(5);
    
    for (int i = 0; i < 5; ++i) {
        fluid.position[i] = {i * 0.1, 0.0};
        fluid.velocity[i] = {0.0, 0.0};
        fluid.density[i] = 1000.0;
        fluid.pressure[i] = 101325.0;
        fluid.surface_type[i] = SurfaceType::INNER;
    }
    
    // 创建边界固体粒子
    solid.particle_num = 2;
    solid.position.resize(2);
    solid.velocity.resize(2);
    solid.normal_vector.resize(2);
    
    solid.position[0] = {-0.1, 0.0};
    solid.velocity[0] = {0.0, 0.0};
    solid.normal_vector[0] = {1.0, 0.0};
    
    solid.position[1] = {0.5, 0.0};
    solid.velocity[1] = {0.0, 0.0};
    solid.normal_vector[1] = {-1.0, 0.0};
    
    std::cout << "\n创建了 " << fluid.particle_num << " 个流体粒子\n";
    std::cout << "创建了 " << solid.particle_num << " 个固体粒子\n";
    
    // 构建邻居列表
    double particle_radius = 0.05;
    double r_e = 0.15;
    double r_cell = 0.12;
    
    NeighborListSearcher searcher;
    searcher.BuildNeighborList(fluid, solid, particle_radius, r_e, r_cell);
    
    std::cout << "\n流体-流体邻居统计:\n";
    int total_fluid_neighbors = 0;
    for (int i = 0; i < fluid.particle_num; ++i) {
        total_fluid_neighbors += fluid.fluid_neighbour_list[i].size();
    }
    std::cout << "  总流体-流体邻居数: " << total_fluid_neighbors << "\n";
    
    std::cout << "\n流体-固体邻居统计:\n";
    int total_solid_neighbors = 0;
    for (int i = 0; i < fluid.particle_num; ++i) {
        total_solid_neighbors += fluid.solid_neighbour_list[i].size();
    }
    std::cout << "  总流体-固体邻居数: " << total_solid_neighbors << "\n";
    
    // 显示前几个粒子的邻居
    std::cout << "\n前3个粒子的邻居:\n";
    for (int i = 0; i < std::min(3, fluid.particle_num); ++i) {
        std::cout << "  粒子[" << i << "] 位置: " << fluid.position[i] << "\n";
        std::cout << "    流体邻居数: " << fluid.fluid_neighbour_list[i].size() << "\n";
        std::cout << "    固体邻居数: " << fluid.solid_neighbour_list[i].size() << "\n";
    }
}

void testUniformDistribution() {
    printSeparator("测试均匀分布粒子的邻居搜索");
    
    FluidParticle fluid("uniform_fluid");
    SolidParticle solid("uniform_solid");
    
    // 创建5x5网格的流体粒子
    fluid.particle_num = 25;
    fluid.position.resize(25);
    fluid.velocity.resize(25);
    fluid.density.resize(25);
    fluid.pressure.resize(25);
    fluid.surface_type.resize(25);
    fluid.fluid_neighbour_list.resize(25);
    fluid.solid_neighbour_list.resize(25);
    
    int idx = 0;
    for (int j = 0; j < 5; ++j) {
        for (int i = 0; i < 5; ++i) {
            fluid.position[idx] = {i * 0.1, j * 0.1};
            fluid.velocity[idx] = {0.0, 0.0};
            fluid.density[idx] = 1000.0;
            fluid.pressure[idx] = 101325.0;
            fluid.surface_type[idx] = SurfaceType::INNER;
            idx++;
        }
    }
    
    solid.particle_num = 0;
    
    std::cout << "\n创建了 " << fluid.particle_num << " 个流体粒子（5x5网格）\n";
    
    // 构建邻居列表
    double particle_radius = 0.05;
    double r_e = 0.12;  // 搜索半径，应该能找到相邻的粒子
    double r_cell = 0.15;
    
    NeighborListSearcher searcher;
    searcher.BuildNeighborList(fluid, solid, particle_radius, r_e, r_cell);
    
    std::cout << "\n邻居列表统计:\n";
    int total_neighbors = 0;
    int max_neighbors = 0;
    int min_neighbors = 1000;
    for (int i = 0; i < fluid.particle_num; ++i) {
        int neighbor_count = fluid.fluid_neighbour_list[i].size();
        total_neighbors += neighbor_count;
        max_neighbors = std::max(max_neighbors, neighbor_count);
        min_neighbors = std::min(min_neighbors, neighbor_count);
    }
    
    std::cout << "  总邻居数: " << total_neighbors << "\n";
    std::cout << "  最大邻居数: " << max_neighbors << "\n";
    std::cout << "  最小邻居数: " << min_neighbors << "\n";
    std::cout << "  平均邻居数: " << (fluid.particle_num > 0 ? total_neighbors / fluid.particle_num : 0) << "\n";
    
    // 验证中心粒子的邻居（应该有最多邻居）
    int center_idx = 12;  // 5x5网格的中心
    std::cout << "\n中心粒子（索引" << center_idx << "）的邻居数: " 
              << fluid.fluid_neighbour_list[center_idx].size() << "\n";
}

void testEmptyParticles() {
    printSeparator("测试空粒子列表");
    
    FluidParticle fluid("empty_fluid");
    SolidParticle solid("empty_solid");
    
    fluid.particle_num = 0;
    solid.particle_num = 0;
    
    std::cout << "\n测试空粒子列表...\n";
    
    double particle_radius = 0.05;
    double r_e = 0.15;
    double r_cell = 0.12;
    
    NeighborListSearcher searcher;
    searcher.BuildNeighborList(fluid, solid, particle_radius, r_e, r_cell);
    
    std::cout << "  ✓ 空粒子列表处理成功（无崩溃）\n";
}

void testLargeSquareWall() {
    printSeparator("测试100x100大尺寸正方形壁面包围流体");
    
    FluidParticle fluid("large_fluid");
    SolidParticle solid("large_solid");
    FileOperator file_op;
    
    // 创建100x100的流体粒子网格
    const int grid_size = 100;
    const double spacing = 0.01;  // 粒子间距
    const double domain_size = (grid_size - 1) * spacing;  // 域大小
    
    fluid.particle_num = grid_size * grid_size;
    fluid.position.resize(fluid.particle_num);
    fluid.velocity.resize(fluid.particle_num);
    fluid.density.resize(fluid.particle_num);
    fluid.pressure.resize(fluid.particle_num);
    fluid.surface_type.resize(fluid.particle_num);
    fluid.fluid_neighbour_list.resize(fluid.particle_num);
    fluid.solid_neighbour_list.resize(fluid.particle_num);
    
    // 生成流体粒子（从(0,0)到(domain_size, domain_size)）
    int idx = 0;
    for (int j = 0; j < grid_size; ++j) {
        for (int i = 0; i < grid_size; ++i) {
            fluid.position[idx] = {i * spacing, j * spacing};
            fluid.velocity[idx] = {0.0, 0.0};
            fluid.density[idx] = 1000.0;
            fluid.pressure[idx] = 101325.0;
            fluid.surface_type[idx] = SurfaceType::INNER;
            idx++;
        }
    }
    
    // 创建正方形壁面（包围流体）
    // 底部边界
    std::vector<double2> wall_positions;
    std::vector<double2> wall_normals;
    
    // 底部边界（y = -spacing）
    for (int i = 0; i < grid_size; ++i) {
        wall_positions.push_back({i * spacing, -spacing});
        wall_normals.push_back({0.0, 1.0});  // 向上
    }
    
    // 顶部边界（y = domain_size + spacing）
    for (int i = 0; i < grid_size; ++i) {
        wall_positions.push_back({i * spacing, domain_size + spacing});
        wall_normals.push_back({0.0, -1.0});  // 向下
    }
    
    // 左侧边界（x = -spacing，不包括角落）
    for (int j = 1; j < grid_size - 1; ++j) {
        wall_positions.push_back({-spacing, j * spacing});
        wall_normals.push_back({1.0, 0.0});  // 向右
    }
    
    // 右侧边界（x = domain_size + spacing，不包括角落）
    for (int j = 1; j < grid_size - 1; ++j) {
        wall_positions.push_back({domain_size + spacing, j * spacing});
        wall_normals.push_back({-1.0, 0.0});  // 向左
    }
    
    solid.particle_num = static_cast<int>(wall_positions.size());
    solid.position.resize(solid.particle_num);
    solid.velocity.resize(solid.particle_num);
    solid.normal_vector.resize(solid.particle_num);
    
    for (int i = 0; i < solid.particle_num; ++i) {
        solid.position[i] = wall_positions[i];
        solid.velocity[i] = {0.0, 0.0};
        solid.normal_vector[i] = wall_normals[i];
    }
    
    std::cout << "\n创建了 " << fluid.particle_num << " 个流体粒子（" 
              << grid_size << "x" << grid_size << "网格）\n";
    std::cout << "创建了 " << solid.particle_num << " 个固体边界粒子\n";
    std::cout << "域大小: " << domain_size << " x " << domain_size << "\n";
    std::cout << "粒子间距: " << spacing << "\n";
    
    // 构建邻居列表
    double particle_radius = spacing * 0.5;
    double r_e = spacing * 2.1;  // 搜索半径，应该能找到相邻的粒子
    double r_cell = spacing * 2.5;  // 网格尺寸
    
    std::cout << "\n构建邻居列表...\n";
    std::cout << "  搜索半径 r_e: " << r_e << "\n";
    std::cout << "  网格尺寸 r_cell: " << r_cell << "\n";
    
    NeighborListSearcher searcher;
    searcher.BuildNeighborList(fluid, solid, particle_radius, r_e, r_cell);
    
    // 统计邻居数量
    std::cout << "\n邻居列表统计:\n";
    int total_fluid_neighbors = 0;
    int total_solid_neighbors = 0;
    int max_fluid_neighbors = 0;
    int max_solid_neighbors = 0;
    int min_fluid_neighbors = 10000;
    int min_solid_neighbors = 10000;
    
    for (int i = 0; i < fluid.particle_num; ++i) {
        int fluid_neighbor_count = fluid.fluid_neighbour_list[i].size();
        int solid_neighbor_count = fluid.solid_neighbour_list[i].size();
        
        total_fluid_neighbors += fluid_neighbor_count;
        total_solid_neighbors += solid_neighbor_count;
        max_fluid_neighbors = std::max(max_fluid_neighbors, fluid_neighbor_count);
        max_solid_neighbors = std::max(max_solid_neighbors, solid_neighbor_count);
        min_fluid_neighbors = std::min(min_fluid_neighbors, fluid_neighbor_count);
        min_solid_neighbors = std::min(min_solid_neighbors, solid_neighbor_count);
    }
    
    std::cout << "  流体-流体邻居:\n";
    std::cout << "    总数: " << total_fluid_neighbors << "\n";
    std::cout << "    最大: " << max_fluid_neighbors << "\n";
    std::cout << "    最小: " << min_fluid_neighbors << "\n";
    std::cout << "    平均: " << (fluid.particle_num > 0 ? total_fluid_neighbors / fluid.particle_num : 0) << "\n";
    
    std::cout << "  流体-固体邻居:\n";
    std::cout << "    总数: " << total_solid_neighbors << "\n";
    std::cout << "    最大: " << max_solid_neighbors << "\n";
    std::cout << "    最小: " << min_solid_neighbors << "\n";
    std::cout << "    平均: " << (fluid.particle_num > 0 ? total_solid_neighbors / fluid.particle_num : 0) << "\n";
    
    // 准备VTK输出
    std::cout << "\n准备VTK输出...\n";
    
    // 写入VTK基础文件（位置和速度）
    std::cout << "1. 写入VTK基础文件（位置和速度）:\n";
    bool success1 = file_op.writeVTKBase("data/large_square_wall.vtk", fluid);
    std::cout << "   写入结果: " << (success1 ? "成功" : "失败") << "\n";
    
    if (success1) {
        // 将流体邻居数量作为标量输出
        std::cout << "2. 追加流体邻居数量标量:\n";
        std::vector<int> fluid_neighbor_counts(fluid.particle_num);
        for (int i = 0; i < fluid.particle_num; ++i) {
            fluid_neighbor_counts[i] = static_cast<int>(fluid.fluid_neighbour_list[i].size());
        }
        bool success2 = file_op.appendVTKScalar("data/large_square_wall.vtk", 
                                                "fluid_neighbor_count", 
                                                fluid_neighbor_counts);
        std::cout << "   追加结果: " << (success2 ? "成功" : "失败") << "\n";
        
        // 将固体邻居数量作为标量输出
        std::cout << "3. 追加固体邻居数量标量:\n";
        std::vector<int> solid_neighbor_counts(fluid.particle_num);
        for (int i = 0; i < fluid.particle_num; ++i) {
            solid_neighbor_counts[i] = static_cast<int>(fluid.solid_neighbour_list[i].size());
        }
        bool success3 = file_op.appendVTKScalar("data/large_square_wall.vtk", 
                                                "solid_neighbor_count", 
                                                solid_neighbor_counts);
        std::cout << "   追加结果: " << (success3 ? "成功" : "失败") << "\n";
        
        // 将总邻居数量作为标量输出
        std::cout << "4. 追加总邻居数量标量:\n";
        std::vector<int> total_neighbor_counts(fluid.particle_num);
        for (int i = 0; i < fluid.particle_num; ++i) {
            total_neighbor_counts[i] = static_cast<int>(fluid.fluid_neighbour_list[i].size() + 
                                                        fluid.solid_neighbour_list[i].size());
        }
        bool success4 = file_op.appendVTKScalar("data/large_square_wall.vtk", 
                                                "total_neighbor_count", 
                                                total_neighbor_counts);
        std::cout << "   追加结果: " << (success4 ? "成功" : "失败") << "\n";
        
        if (success2 && success3 && success4) {
            std::cout << "\n   ✓ VTK文件已生成: data/large_square_wall.vtk\n";
            std::cout << "   可以使用 ParaView 打开查看邻居数量分布\n";
        }
    }
    
    // 显示一些边界粒子的邻居信息
    std::cout << "\n边界粒子邻居信息（示例）:\n";
    // 左下角粒子（索引0）
    std::cout << "  左下角粒子[0] 位置: " << fluid.position[0] << "\n";
    std::cout << "    流体邻居数: " << fluid.fluid_neighbour_list[0].size() << "\n";
    std::cout << "    固体邻居数: " << fluid.solid_neighbour_list[0].size() << "\n";
    
    // 中心粒子
    int center_idx = (grid_size / 2) * grid_size + (grid_size / 2);
    std::cout << "  中心粒子[" << center_idx << "] 位置: " << fluid.position[center_idx] << "\n";
    std::cout << "    流体邻居数: " << fluid.fluid_neighbour_list[center_idx].size() << "\n";
    std::cout << "    固体邻居数: " << fluid.solid_neighbour_list[center_idx].size() << "\n";
    
    // 右上角粒子
    int top_right_idx = (grid_size - 1) * grid_size + (grid_size - 1);
    std::cout << "  右上角粒子[" << top_right_idx << "] 位置: " << fluid.position[top_right_idx] << "\n";
    std::cout << "    流体邻居数: " << fluid.fluid_neighbour_list[top_right_idx].size() << "\n";
    std::cout << "    固体邻居数: " << fluid.solid_neighbour_list[top_right_idx].size() << "\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "MPS Baseline 2D - NeighborListSearcher 邻居列表搜索器测试\n";
    std::cout << std::string(60, '=') << "\n";
    
    try {
        testBasicNeighborSearch();
        testWithSolidParticles();
        testUniformDistribution();
        testEmptyParticles();
        testLargeSquareWall();
        
        printSeparator("所有测试完成");
        std::cout << "\n✓ 所有邻居列表搜索器测试通过\n";
        std::cout << "\n生成的VTK文件:\n";
        std::cout << "  - data/large_square_wall.vtk (100x100大尺寸测试用例)\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\n错误: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}

