#include <iostream>
#include <iomanip>
#include <set>
#include "../src/core/Particle.hpp"
#include "../src/core/FileOperator.hpp"
#include "../include/core/MPSUtils.h"

using namespace mps2D;

void printSeparator(const std::string& title) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "  " << title << "\n";
    std::cout << std::string(60, '=') << "\n";
}

void testFluidParticleRead() {
    printSeparator("测试 FluidParticle 读取功能");
    
    // 测试从TXT文件读取
    std::cout << "\n1. 从TXT文件读取 FluidParticle:\n";
    FluidParticle fluid_txt("fluid_txt");
    int count_txt = fluid_txt.getParticleFromFile("data/fluid_particles_2d.txt");
    std::cout << "   读取的粒子数量: " << count_txt << "\n";
    
    if (count_txt > 0) {
        std::cout << "   前3个粒子位置:\n";
        for (int i = 0; i < std::min(3, count_txt); ++i) {
            std::cout << "     粒子[" << i << "] 位置: " << fluid_txt.position[i] << "\n";
        }
    }
    
    // 测试从CSV文件读取
    std::cout << "\n2. 从CSV文件读取 FluidParticle:\n";
    FluidParticle fluid_csv("fluid_csv");
    int count_csv = fluid_csv.getParticleFromFile("data/fluid_particles_2d.csv");
    std::cout << "   读取的粒子数量: " << count_csv << "\n";
    
    if (count_csv > 0) {
        std::cout << "   前3个粒子位置:\n";
        for (int i = 0; i < std::min(3, count_csv); ++i) {
            std::cout << "     粒子[" << i << "] 位置: " << fluid_csv.position[i] << "\n";
        }
    }
    
    // 验证读取的数据是否一致
    if (count_txt == count_csv && count_txt > 0) {
        bool data_match = true;
        for (int i = 0; i < count_txt; ++i) {
            double dist = ComputeDistance(fluid_txt.position[i], fluid_csv.position[i]);
            if (dist > 1e-10) {
                data_match = false;
                break;
            }
        }
        if (data_match) {
            std::cout << "\n   ✓ TXT和CSV文件读取的数据一致\n";
        } else {
            std::cout << "\n   ✗ TXT和CSV文件读取的数据不一致\n";
        }
    }
}

void testSolidParticleRead() {
    printSeparator("测试 SolidParticle 读取功能");
    
    // 测试从TXT文件读取
    std::cout << "\n1. 从TXT文件读取 SolidParticle:\n";
    SolidParticle solid_txt("solid_txt");
    int count = solid_txt.getParticleFromFile("data/solid_particles_2d.txt");
    std::cout << "   读取的粒子数量: " << count << "\n";
    
    if (count > 0 && count <= static_cast<int>(solid_txt.position.size()) && 
        count <= static_cast<int>(solid_txt.normal_vector.size())) {
        std::cout << "   前3个粒子信息:\n";
        int show_count = std::min(3, count);
        for (int i = 0; i < show_count; ++i) {
            std::cout << "     粒子[" << i << "] 位置: " << solid_txt.position[i] << "\n";
            std::cout << "     粒子[" << i << "] 法向向量: " << solid_txt.normal_vector[i] << "\n";
        }
    } else if (count > 0) {
        std::cout << "   ⚠ 数据大小不匹配\n";
    }
}

void testParticleWrite() {
    printSeparator("测试粒子数据写入功能");
    
    FileOperator file_op;
    
    // 读取数据
    FluidParticle fluid("fluid_write");
    fluid.getParticleFromFile("data/fluid_particles_2d.txt");
    std::cout << "读取了 " << fluid.particle_num << " 个流体粒子\n";
    
    SolidParticle solid("solid_write");
    solid.getParticleFromFile("data/solid_particles_2d.txt");
    std::cout << "读取了 " << solid.particle_num << " 个固体粒子\n";
    
    // 写入TXT文件
    std::cout << "\n1. 写入 FluidParticle 到 TXT 文件:\n";
    bool success1 = file_op.writeParticleToFile("data/output_fluid_2d.txt", fluid);
    std::cout << "   写入结果: " << (success1 ? "成功" : "失败") << "\n";
    
    // 写入CSV文件
    std::cout << "\n2. 写入 FluidParticle 到 CSV 文件:\n";
    bool success2 = file_op.writeParticleToFile("data/output_fluid_2d.csv", fluid);
    std::cout << "   写入结果: " << (success2 ? "成功" : "失败") << "\n";
    
    // 写入SolidParticle到TXT文件（包含法向向量）
    std::cout << "\n3. 写入 SolidParticle 到 TXT 文件（包含法向向量）:\n";
    std::set<std::string> custom_fields = {"normal_vector"};
    bool success3 = file_op.writeParticleToFile("data/output_solid_2d.txt", solid, custom_fields);
    std::cout << "   写入结果: " << (success3 ? "成功" : "失败") << "\n";
    
    // 写入SolidParticle到CSV文件（包含法向向量）
    std::cout << "\n4. 写入 SolidParticle 到 CSV 文件（包含法向向量）:\n";
    bool success4 = file_op.writeParticleToFile("data/output_solid_2d.csv", solid, custom_fields);
    std::cout << "   写入结果: " << (success4 ? "成功" : "失败") << "\n";
}

void testUniformDistributionReadAndVTKWrite() {
    printSeparator("测试均匀分布粒子读取和VTK写入");
    
    FileOperator file_op;
    
    // 读取均匀分布的流体粒子
    std::cout << "\n1. 读取均匀分布的流体粒子数据:\n";
    FluidParticle fluid("uniform_fluid");
    int count = fluid.getParticleFromFile("data/fluid_particles_2d.txt");
    std::cout << "   读取的粒子数量: " << count << "\n";
    
    if (count == 0) {
        std::cout << "   ✗ 未能读取粒子数据\n";
        return;
    }
    
    // 验证粒子分布
    std::cout << "\n2. 验证粒子分布（5x5网格）:\n";
    if (count == 25) {
        std::cout << "   ✓ 粒子数量正确（5x5=25）\n";
    } else {
        std::cout << "   ⚠ 粒子数量: " << count << "（预期25）\n";
    }
    
    // 检查粒子位置范围
    double min_x = 1e10, max_x = -1e10;
    double min_y = 1e10, max_y = -1e10;
    for (int i = 0; i < count; ++i) {
        if (fluid.position[i].x < min_x) min_x = fluid.position[i].x;
        if (fluid.position[i].x > max_x) max_x = fluid.position[i].x;
        if (fluid.position[i].y < min_y) min_y = fluid.position[i].y;
        if (fluid.position[i].y > max_y) max_y = fluid.position[i].y;
    }
    std::cout << "   位置范围: x=[" << min_x << ", " << max_x << "], y=[" 
              << min_y << ", " << max_y << "]\n";
    
    // 设置一些测试数据（密度、压力等）
    std::cout << "\n3. 设置测试数据:\n";
    fluid.density.resize(count);
    fluid.pressure.resize(count);
    fluid.surface_type.resize(count);
    
    for (int i = 0; i < count; ++i) {
        fluid.density[i] = 1000.0;
        fluid.pressure[i] = 101325.0 + i * 10.0;
        fluid.surface_type[i] = (i < 5 || i >= 20 || i % 5 == 0 || i % 5 == 4) 
                                ? SurfaceType::SURFACE 
                                : SurfaceType::INNER;
    }
    std::cout << "   已设置密度、压力和表面类型\n";
    
    // 写入VTK基础文件（位置和速度）
    std::cout << "\n4. 写入VTK基础文件（位置和速度）:\n";
    bool success1 = file_op.writeVTKBase("data/uniform_fluid_2d.vtk", fluid);
    std::cout << "   写入结果: " << (success1 ? "成功" : "失败") << "\n";
    
    if (success1) {
        // 追加密度标量
        std::cout << "\n5. 追加密度标量到VTK文件:\n";
        bool success2 = file_op.appendVTKScalar("data/uniform_fluid_2d.vtk", "density", fluid.density);
        std::cout << "   追加结果: " << (success2 ? "成功" : "失败") << "\n";
        
        // 追加压力标量
        std::cout << "\n6. 追加压力标量到VTK文件:\n";
        bool success3 = file_op.appendVTKScalar("data/uniform_fluid_2d.vtk", "pressure", fluid.pressure);
        std::cout << "   追加结果: " << (success3 ? "成功" : "失败") << "\n";
        
        // 追加表面类型标量
        std::cout << "\n7. 追加表面类型标量到VTK文件:\n";
        std::vector<int> surface_types(count);
        for (int i = 0; i < count; ++i) {
            surface_types[i] = static_cast<int>(fluid.surface_type[i]);
        }
        bool success4 = file_op.appendVTKScalar("data/uniform_fluid_2d.vtk", "surface_type", surface_types);
        std::cout << "   追加结果: " << (success4 ? "成功" : "失败") << "\n";
        
        std::cout << "\n   ✓ VTK文件已生成: data/uniform_fluid_2d.vtk\n";
        std::cout << "   可以使用 ParaView 或其他VTK可视化工具打开查看\n";
    }
}

void testVectorWrite() {
    printSeparator("测试向量写入功能（Debug）");
    
    FileOperator file_op;
    
    // 测试写入double向量
    std::vector<double> double_vec = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::cout << "\n1. 写入 double 向量到 TXT 文件:\n";
    bool success1 = file_op.writeVectorToFile("data/debug_double_2d.txt", double_vec);
    std::cout << "   写入结果: " << (success1 ? "成功" : "失败") << "\n";
    
    // 测试写入double2向量
    std::vector<double2> double2_vec = {
        {1.0, 2.0},
        {3.0, 4.0},
        {5.0, 6.0}
    };
    std::cout << "\n2. 写入 double2 向量到 TXT 文件:\n";
    bool success2 = file_op.writeVectorToFile("data/debug_double2_2d.txt", double2_vec);
    std::cout << "   写入结果: " << (success2 ? "成功" : "失败") << "\n";
    
    // 测试写入double2向量到CSV文件
    std::cout << "\n3. 写入 double2 向量到 CSV 文件:\n";
    bool success3 = file_op.writeVectorToFile("data/debug_double2_2d.csv", double2_vec);
    std::cout << "   写入结果: " << (success3 ? "成功" : "失败") << "\n";
}

void testVTKWriteWithSolidParticle() {
    printSeparator("测试 SolidParticle VTK 写入功能");
    
    FileOperator file_op;
    
    // 读取固体粒子数据
    SolidParticle solid("solid_vtk");
    solid.getParticleFromFile("data/solid_particles_2d.txt");
    std::cout << "读取了 " << solid.particle_num << " 个固体粒子\n";
    
    if (solid.particle_num == 0) {
        std::cout << "   ✗ 未能读取粒子数据\n";
        return;
    }
    
    // 写入VTK基础文件
    std::cout << "\n1. 写入 SolidParticle 基础 VTK 文件（位置和速度）:\n";
    bool success1 = file_op.writeVTKBase("data/solid_particles_2d.vtk", solid);
    std::cout << "   写入结果: " << (success1 ? "成功" : "失败") << "\n";
    
    if (success1) {
        // 追加法向向量
        std::cout << "\n2. 追加法向向量到VTK文件:\n";
        bool success2 = file_op.appendVTKVector("data/solid_particles_2d.vtk", "normal_vector", solid.normal_vector);
        std::cout << "   追加结果: " << (success2 ? "成功" : "失败") << "\n";
        
        std::cout << "\n   ✓ VTK文件已生成: data/solid_particles_2d.vtk\n";
    }
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "MPS Baseline 2D - FileOperator 文件读写模块测试\n";
    std::cout << std::string(60, '=') << "\n";
    
    try {
        // 测试文件读取功能
        testFluidParticleRead();
        testSolidParticleRead();
        
        // 测试文件写入功能
        testParticleWrite();
        
        // 测试均匀分布粒子读取和VTK写入（主要测试用例）
        testUniformDistributionReadAndVTKWrite();
        
        // 测试向量写入功能
        testVectorWrite();
        
        // 测试SolidParticle VTK写入
        testVTKWriteWithSolidParticle();
        
        printSeparator("所有测试完成");
        std::cout << "\n测试结果文件已生成在 data/ 目录下:\n";
        std::cout << "  - output_fluid_2d.txt / output_fluid_2d.csv\n";
        std::cout << "  - output_solid_2d.txt / output_solid_2d.csv\n";
        std::cout << "  - uniform_fluid_2d.vtk (均匀分布粒子VTK文件)\n";
        std::cout << "  - solid_particles_2d.vtk (固体粒子VTK文件)\n";
        std::cout << "  - debug_double_2d.txt / debug_double2_2d.txt / debug_double2_2d.csv\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\n错误: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}

