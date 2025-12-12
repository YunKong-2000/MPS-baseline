#include <iostream>
#include <iomanip>
#include <vector>
#include "../src/core/Particle.hpp"
#include "../src/core/FileOperator.hpp"
#include "../src/neighbour_list/NeighborListSearcher.hpp"
#include "../src/surface_detection/SurfaceDetector.hpp"
#include "../include/core/MPSUtils.h"

using namespace mps2D;

void printSeparator(const std::string& title) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "  " << title << "\n";
    std::cout << std::string(60, '=') << "\n";
}

void testSquareContainer() {
    printSeparator("测试1：方形容器盛装液体");
    
    FluidParticle fluid("square_container_fluid");
    SolidParticle solid("square_container_solid");
    FileOperator file_op;
    
    // 创建方形容器：50x50的流体粒子网格，容器壁面包围
    const int grid_size = 50;
    const double spacing = 0.02;
    const double domain_size = (grid_size - 1) * spacing;
    
    fluid.particle_num = grid_size * grid_size;
    fluid.position.resize(fluid.particle_num);
    fluid.velocity.resize(fluid.particle_num);
    fluid.density.resize(fluid.particle_num);
    fluid.pressure.resize(fluid.particle_num);
    fluid.surface_type.resize(fluid.particle_num);
    fluid.fluid_neighbour_list.resize(fluid.particle_num);
    fluid.solid_neighbour_list.resize(fluid.particle_num);
    
    int idx = 0;
    for (int j = 0; j < grid_size; ++j) {
        for (int i = 0; i < grid_size; ++i) {
            double2 pos = {i * spacing, j * spacing};
            fluid.position[idx] = pos;
            fluid.velocity[idx] = {0.0, 0.0};
            fluid.density[idx] = 1000.0;
            fluid.pressure[idx] = 101325.0;
            fluid.surface_type[idx] = SurfaceType::INNER;
            idx++;
        }
    }
    
    // 创建方形容器壁面（顶部开口，其他三面封闭）
    std::vector<double2> wall_positions;
    std::vector<double2> wall_normals;
    
    // 底部边界（完整，包括角落）
    for (int i = 0; i < grid_size; ++i) {
        wall_positions.push_back({i * spacing, -spacing});
        wall_normals.push_back({0.0, 1.0});
    }
    
    // 左侧边界（完整，从底部到顶部，包括角落）
    for (int j = 0; j < grid_size; ++j) {
        wall_positions.push_back({-spacing, j * spacing});
        wall_normals.push_back({1.0, 0.0});
    }
    
    // 右侧边界（完整，从底部到顶部，包括角落）
    for (int j = 0; j < grid_size; ++j) {
        wall_positions.push_back({domain_size + spacing, j * spacing});
        wall_normals.push_back({-1.0, 0.0});
    }
    
    // 注意：顶部不添加壁面，形成开口
    
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
    std::cout << "创建了 " << solid.particle_num << " 个固体容器壁面粒子\n";
    
    // 构建邻居列表
    double particle_radius = spacing * 0.5;
    double r_e = spacing * 2.1;
    double r_cell = spacing * 2.5;
    
    std::cout << "\n构建邻居列表...\n";
    NeighborListSearcher searcher;
    searcher.BuildNeighborList(fluid, solid, particle_radius, r_e, r_cell);
    
    // 检测自由面
    double smoothing_radius = r_e;
    std::cout << "检测自由面...\n";
    SurfaceDetector detector;
    detector.DetectSurfaceParticles(fluid, solid, smoothing_radius, spacing);
    
    // 统计结果
    std::cout << "\n自由面检测结果统计:\n";
    int inner_count = 0, surface_count = 0, near_surface_count = 0, splash_count = 0;
    for (int i = 0; i < fluid.particle_num; ++i) {
        switch (fluid.surface_type[i]) {
            case SurfaceType::INNER:
                inner_count++;
                break;
            case SurfaceType::SURFACE:
                surface_count++;
                break;
            case SurfaceType::NEAR_SURFACE:
                near_surface_count++;
                break;
            case SurfaceType::SPLASH:
                splash_count++;
                break;
        }
    }
    
    std::cout << "  内部粒子: " << inner_count << "\n";
    std::cout << "  自由面粒子: " << surface_count << "\n";
    std::cout << "  近自由面粒子: " << near_surface_count << "\n";
    std::cout << "  飞溅粒子: " << splash_count << "\n";
    
    // 计算阴影面积比例
    std::vector<double> shadow_ratios(fluid.particle_num);
    for (int i = 0; i < fluid.particle_num; ++i) {
        shadow_ratios[i] = detector.ComputeShadowAreaRatio(
            i, fluid, solid, smoothing_radius, spacing);
    }
    
    // 写入VTK文件
    std::cout << "\n写入VTK文件...\n";
    bool success1 = file_op.writeVTKBase("data/square_container.vtk", fluid);
    if (success1) {
        std::vector<int> surface_types(fluid.particle_num);
        for (int i = 0; i < fluid.particle_num; ++i) {
            surface_types[i] = static_cast<int>(fluid.surface_type[i]);
        }
        file_op.appendVTKScalar("data/square_container.vtk", "surface_type", surface_types);
        file_op.appendVTKScalar("data/square_container.vtk", "shadow_area_ratio", shadow_ratios);
        std::cout << "  ✓ 流体粒子VTK文件已生成: data/square_container.vtk\n";
    }
    
    // 输出壁面（固体粒子）VTK文件
    if (solid.particle_num > 0) {
        bool success2 = file_op.writeVTKBase("data/square_container_wall.vtk", solid);
        if (success2) {
            file_op.appendVTKVector("data/square_container_wall.vtk", "normal_vector", solid.normal_vector);
            std::cout << "  ✓ 壁面VTK文件已生成: data/square_container_wall.vtk\n";
        }
    }
}

void testComplexContainer() {
    printSeparator("测试2：复杂外形容器装水");
    
    FluidParticle fluid("complex_container_fluid");
    SolidParticle solid("complex_container_solid");
    FileOperator file_op;
    
    // 创建圆形/复杂形状容器：使用圆形区域内的流体粒子
    const double spacing = 0.02;
    const double center_x = 0.5;
    const double center_y = 0.5;
    const double radius = 0.4;
    const double wall_thickness = spacing * 2.0;
    
    // 生成圆形区域内的流体粒子
    std::vector<double2> fluid_positions;
    const int grid_size = 50;
    
    for (int j = 0; j < grid_size; ++j) {
        for (int i = 0; i < grid_size; ++i) {
            double x = i * spacing;
            double y = j * spacing;
            double dx = x - center_x;
            double dy = y - center_y;
            double dist = std::sqrt(dx * dx + dy * dy);
            
            // 在圆形区域内
            if (dist < radius) {
                fluid_positions.push_back({x, y});
            }
        }
    }
    
    fluid.particle_num = static_cast<int>(fluid_positions.size());
    fluid.position.resize(fluid.particle_num);
    fluid.velocity.resize(fluid.particle_num);
    fluid.density.resize(fluid.particle_num);
    fluid.pressure.resize(fluid.particle_num);
    fluid.surface_type.resize(fluid.particle_num);
    fluid.fluid_neighbour_list.resize(fluid.particle_num);
    fluid.solid_neighbour_list.resize(fluid.particle_num);
    
    for (int i = 0; i < fluid.particle_num; ++i) {
        fluid.position[i] = fluid_positions[i];
        fluid.velocity[i] = {0.0, 0.0};
        fluid.density[i] = 1000.0;
        fluid.pressure[i] = 101325.0;
        fluid.surface_type[i] = SurfaceType::INNER;
    }
    
    // 创建圆形容器壁面
    std::vector<double2> wall_positions;
    std::vector<double2> wall_normals;
    
    const int num_wall_particles = static_cast<int>(2.0 * M_PI * (radius + wall_thickness) / spacing);
    for (int i = 0; i < num_wall_particles; ++i) {
        double angle = 2.0 * M_PI * i / num_wall_particles;
        double wall_radius = radius + wall_thickness;
        double x = center_x + wall_radius * std::cos(angle);
        double y = center_y + wall_radius * std::sin(angle);
        
        // 计算法向向量（指向圆心）
        double2 normal = {-std::cos(angle), -std::sin(angle)};
        double norm = std::sqrt(normal.x * normal.x + normal.y * normal.y);
        if (norm > 1e-10) {
            normal.x /= norm;
            normal.y /= norm;
        }
        
        wall_positions.push_back({x, y});
        wall_normals.push_back(normal);
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
    
    std::cout << "\n创建了 " << fluid.particle_num << " 个流体粒子（圆形区域）\n";
    std::cout << "创建了 " << solid.particle_num << " 个固体容器壁面粒子（圆形）\n";
    std::cout << "容器半径: " << radius << "\n";
    
    // 构建邻居列表
    double particle_radius = spacing * 0.5;
    double r_e = spacing * 2.1;
    double r_cell = spacing * 2.5;
    
    std::cout << "\n构建邻居列表...\n";
    NeighborListSearcher searcher;
    searcher.BuildNeighborList(fluid, solid, particle_radius, r_e, r_cell);
    
    // 检测自由面
    double smoothing_radius = r_e;
    std::cout << "检测自由面...\n";
    SurfaceDetector detector;
    detector.DetectSurfaceParticles(fluid, solid, smoothing_radius, spacing);
    
    // 统计结果
    std::cout << "\n自由面检测结果统计:\n";
    int inner_count = 0, surface_count = 0, near_surface_count = 0, splash_count = 0;
    for (int i = 0; i < fluid.particle_num; ++i) {
        switch (fluid.surface_type[i]) {
            case SurfaceType::INNER:
                inner_count++;
                break;
            case SurfaceType::SURFACE:
                surface_count++;
                break;
            case SurfaceType::NEAR_SURFACE:
                near_surface_count++;
                break;
            case SurfaceType::SPLASH:
                splash_count++;
                break;
        }
    }
    
    std::cout << "  内部粒子: " << inner_count << "\n";
    std::cout << "  自由面粒子: " << surface_count << "\n";
    std::cout << "  近自由面粒子: " << near_surface_count << "\n";
    std::cout << "  飞溅粒子: " << splash_count << "\n";
    
    // 计算阴影面积比例
    std::vector<double> shadow_ratios(fluid.particle_num);
    for (int i = 0; i < fluid.particle_num; ++i) {
        shadow_ratios[i] = detector.ComputeShadowAreaRatio(
            i, fluid, solid, smoothing_radius, spacing);
    }
    
    // 写入VTK文件
    std::cout << "\n写入VTK文件...\n";
    bool success1 = file_op.writeVTKBase("data/complex_container.vtk", fluid);
    if (success1) {
        std::vector<int> surface_types(fluid.particle_num);
        for (int i = 0; i < fluid.particle_num; ++i) {
            surface_types[i] = static_cast<int>(fluid.surface_type[i]);
        }
        file_op.appendVTKScalar("data/complex_container.vtk", "surface_type", surface_types);
        file_op.appendVTKScalar("data/complex_container.vtk", "shadow_area_ratio", shadow_ratios);
        std::cout << "  ✓ 流体粒子VTK文件已生成: data/complex_container.vtk\n";
    }
    
    // 输出壁面（固体粒子）VTK文件
    if (solid.particle_num > 0) {
        bool success2 = file_op.writeVTKBase("data/complex_container_wall.vtk", solid);
        if (success2) {
            file_op.appendVTKVector("data/complex_container_wall.vtk", "normal_vector", solid.normal_vector);
            std::cout << "  ✓ 壁面VTK文件已生成: data/complex_container_wall.vtk\n";
        }
    }
}

void testLargeDeformationSurface() {
    printSeparator("测试3：具有大变形液面的场景");
    
    FluidParticle fluid("large_deformation_fluid");
    SolidParticle solid("large_deformation_solid");
    FileOperator file_op;
    
    // 创建具有波浪/大变形液面的场景
    const double spacing = 0.015;
    const int grid_x = 60;
    const int grid_y = 40;
    const double domain_width = (grid_x - 1) * spacing;
    const double domain_height = (grid_y - 1) * spacing;
    
    // 波浪参数
    const double wave_amplitude = 0.15;  // 波浪幅度
    const double wave_frequency = 2.0 * M_PI / domain_width * 2.0;  // 波浪频率
    const double liquid_level = domain_height * 0.6;  // 基础液面高度
    
    // 生成具有波浪形状的流体粒子
    std::vector<double2> fluid_positions;
    
    for (int j = 0; j < grid_y; ++j) {
        for (int i = 0; i < grid_x; ++i) {
            double x = i * spacing;
            double base_y = j * spacing;
            
            // 创建波浪形状的自由面
            // 使用正弦波创建大变形液面
            double wave_y = base_y + wave_amplitude * std::sin(wave_frequency * x);
            
            // 添加第二个波浪分量（创建更复杂的变形）
            double wave_y2 = wave_amplitude * 0.5 * std::sin(wave_frequency * x * 2.0);
            wave_y += wave_y2;
            
            // 只保留在波浪下方的粒子（模拟液体）
            if (wave_y < liquid_level) {
                fluid_positions.push_back({x, base_y});
            }
        }
    }
    
    fluid.particle_num = static_cast<int>(fluid_positions.size());
    fluid.position.resize(fluid.particle_num);
    fluid.velocity.resize(fluid.particle_num);
    fluid.density.resize(fluid.particle_num);
    fluid.pressure.resize(fluid.particle_num);
    fluid.surface_type.resize(fluid.particle_num);
    fluid.fluid_neighbour_list.resize(fluid.particle_num);
    fluid.solid_neighbour_list.resize(fluid.particle_num);
    
    for (int i = 0; i < fluid.particle_num; ++i) {
        fluid.position[i] = fluid_positions[i];
        fluid.velocity[i] = {0.0, 0.0};
        fluid.density[i] = 1000.0;
        fluid.pressure[i] = 101325.0;
        fluid.surface_type[i] = SurfaceType::INNER;
    }
    
    // 创建底部和侧面的容器壁面
    std::vector<double2> wall_positions;
    std::vector<double2> wall_normals;
    
    // 底部边界
    for (int i = 0; i < grid_x; ++i) {
        wall_positions.push_back({i * spacing, -spacing});
        wall_normals.push_back({0.0, 1.0});
    }
    
    // 左侧边界
    for (int j = 0; j < grid_y; ++j) {
        wall_positions.push_back({-spacing, j * spacing});
        wall_normals.push_back({1.0, 0.0});
    }
    
    // 右侧边界
    for (int j = 0; j < grid_y; ++j) {
        wall_positions.push_back({domain_width + spacing, j * spacing});
        wall_normals.push_back({-1.0, 0.0});
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
    
    std::cout << "\n创建了 " << fluid.particle_num << " 个流体粒子（波浪形状）\n";
    std::cout << "创建了 " << solid.particle_num << " 个固体容器壁面粒子\n";
    std::cout << "波浪幅度: " << wave_amplitude << "\n";
    std::cout << "基础液面高度: " << liquid_level << "\n";
    
    // 构建邻居列表
    double particle_radius = spacing * 0.5;
    double r_e = spacing * 2.1;
    double r_cell = spacing * 2.5;
    
    std::cout << "\n构建邻居列表...\n";
    NeighborListSearcher searcher;
    searcher.BuildNeighborList(fluid, solid, particle_radius, r_e, r_cell);
    
    // 检测自由面
    double smoothing_radius = r_e;
    std::cout << "检测自由面...\n";
    SurfaceDetector detector;
    detector.DetectSurfaceParticles(fluid, solid, smoothing_radius, spacing);
    
    // 统计结果
    std::cout << "\n自由面检测结果统计:\n";
    int inner_count = 0, surface_count = 0, near_surface_count = 0, splash_count = 0;
    for (int i = 0; i < fluid.particle_num; ++i) {
        switch (fluid.surface_type[i]) {
            case SurfaceType::INNER:
                inner_count++;
                break;
            case SurfaceType::SURFACE:
                surface_count++;
                break;
            case SurfaceType::NEAR_SURFACE:
                near_surface_count++;
                break;
            case SurfaceType::SPLASH:
                splash_count++;
                break;
        }
    }
    
    std::cout << "  内部粒子: " << inner_count << "\n";
    std::cout << "  自由面粒子: " << surface_count << "\n";
    std::cout << "  近自由面粒子: " << near_surface_count << "\n";
    std::cout << "  飞溅粒子: " << splash_count << "\n";
    
    // 计算阴影面积比例
    std::vector<double> shadow_ratios(fluid.particle_num);
    for (int i = 0; i < fluid.particle_num; ++i) {
        shadow_ratios[i] = detector.ComputeShadowAreaRatio(
            i, fluid, solid, smoothing_radius, spacing);
    }
    
    // 写入VTK文件
    std::cout << "\n写入VTK文件...\n";
    bool success1 = file_op.writeVTKBase("data/large_deformation_surface.vtk", fluid);
    if (success1) {
        std::vector<int> surface_types(fluid.particle_num);
        for (int i = 0; i < fluid.particle_num; ++i) {
            surface_types[i] = static_cast<int>(fluid.surface_type[i]);
        }
        file_op.appendVTKScalar("data/large_deformation_surface.vtk", "surface_type", surface_types);
        file_op.appendVTKScalar("data/large_deformation_surface.vtk", "shadow_area_ratio", shadow_ratios);
        std::cout << "  ✓ 流体粒子VTK文件已生成: data/large_deformation_surface.vtk\n";
    }
    
    // 输出壁面（固体粒子）VTK文件
    if (solid.particle_num > 0) {
        bool success2 = file_op.writeVTKBase("data/large_deformation_surface_wall.vtk", solid);
        if (success2) {
            file_op.appendVTKVector("data/large_deformation_surface_wall.vtk", "normal_vector", solid.normal_vector);
            std::cout << "  ✓ 壁面VTK文件已生成: data/large_deformation_surface_wall.vtk\n";
        }
    }
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "MPS Baseline 2D - SurfaceDetector 自由面检测测试\n";
    std::cout << std::string(60, '=') << "\n";
    
    try {
        testSquareContainer();
        testComplexContainer();
        testLargeDeformationSurface();
        
        printSeparator("所有测试完成");
        std::cout << "\n✓ 所有自由面检测测试通过\n";
        std::cout << "\n生成的VTK文件:\n";
        std::cout << "  流体粒子:\n";
        std::cout << "    - data/square_container.vtk (方形容器测试用例)\n";
        std::cout << "    - data/complex_container.vtk (复杂外形容器测试用例)\n";
        std::cout << "    - data/large_deformation_surface.vtk (大变形液面测试用例)\n";
        std::cout << "  壁面（固体粒子）:\n";
        std::cout << "    - data/square_container_wall.vtk (方形容器壁面)\n";
        std::cout << "    - data/complex_container_wall.vtk (复杂外形容器壁面)\n";
        std::cout << "    - data/large_deformation_surface_wall.vtk (大变形液面壁面)\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\n错误: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}

