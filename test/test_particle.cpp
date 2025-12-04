#include <iostream>
#include <iomanip>
#include <set>
#include "../src/core/Particle.hpp"
#include "../src/core/FileOperator.hpp"

using namespace mps;

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
    fluid_txt.getParticleFromFile("data/fluid_particles.txt");
    std::cout << "   读取的粒子数量: " << fluid_txt.particle_num << "\n";
    
    // 测试从CSV文件读取
    std::cout << "\n2. 从CSV文件读取 FluidParticle:\n";
    FluidParticle fluid_csv("fluid_csv");
    fluid_csv.getParticleFromFile("data/fluid_particles.csv");
    std::cout << "   读取的粒子数量: " << fluid_csv.particle_num << "\n";
    
    // 测试手动读取
    std::cout << "\n3. 手动调用 getParticleFromFile:\n";
    FluidParticle fluid_manual("fluid_manual");
    fluid_manual.getParticleFromFile("data/fluid_particles.txt");
    std::cout << "   读取的粒子数量: " << fluid_manual.particle_num << "\n";
}

void testSolidParticleRead() {
    printSeparator("测试 SolidParticle 读取功能");
    
    // 测试从TXT文件读取
    std::cout << "\n1. 从TXT文件读取 SolidParticle:\n";
    SolidParticle solid_txt("solid_txt");
    solid_txt.getParticleFromFile("data/solid_particles.txt");
    std::cout << "   读取的粒子数量: " << solid_txt.particle_num << "\n";
    
    // 测试从CSV文件读取
    std::cout << "\n2. 从CSV文件读取 SolidParticle:\n";
    SolidParticle solid_csv("solid_csv");
    solid_csv.getParticleFromFile("data/solid_particles.csv");
    std::cout << "   读取的粒子数量: " << solid_csv.particle_num << "\n";
    
    // 测试手动读取
    std::cout << "\n3. 手动调用 getParticleFromFile:\n";
    SolidParticle solid_manual("solid_manual");
    solid_manual.getParticleFromFile("data/solid_particles.txt");
    std::cout << "   读取的粒子数量: " << solid_manual.particle_num << "\n";
}

void testParticleWrite() {
    printSeparator("测试粒子数据写入功能");
    
    FileOperator file_op;
    
    // 读取数据
    FluidParticle fluid("fluid_write");
    fluid.getParticleFromFile("data/fluid_particles.txt");
    std::cout << "读取了 " << fluid.particle_num << " 个流体粒子\n";
    
    SolidParticle solid("solid_write");
    solid.getParticleFromFile("data/solid_particles.txt");
    std::cout << "读取了 " << solid.particle_num << " 个固体粒子\n";
    
    // 写入TXT文件
    std::cout << "\n1. 写入 FluidParticle 到 TXT 文件:\n";
    bool success1 = file_op.writeParticleToFile("data/output_fluid.txt", fluid);
    std::cout << "   写入结果: " << (success1 ? "成功" : "失败") << "\n";
    
    // 写入CSV文件
    std::cout << "\n2. 写入 FluidParticle 到 CSV 文件:\n";
    bool success2 = file_op.writeParticleToFile("data/output_fluid.csv", fluid);
    std::cout << "   写入结果: " << (success2 ? "成功" : "失败") << "\n";
    
    // 写入SolidParticle到TXT文件（包含法向向量）
    std::cout << "\n3. 写入 SolidParticle 到 TXT 文件（包含法向向量）:\n";
    std::set<std::string> custom_fields = {"normal_vector"};
    bool success3 = file_op.writeParticleToFile("data/output_solid.txt", solid, custom_fields);
    std::cout << "   写入结果: " << (success3 ? "成功" : "失败") << "\n";
    
    // 写入SolidParticle到CSV文件（包含法向向量）
    std::cout << "\n4. 写入 SolidParticle 到 CSV 文件（包含法向向量）:\n";
    bool success4 = file_op.writeParticleToFile("data/output_solid.csv", solid, custom_fields);
    std::cout << "   写入结果: " << (success4 ? "成功" : "失败") << "\n";
}

void testParticleCopy() {
    printSeparator("测试粒子数据拷贝功能");
    
    // 读取源粒子
    FluidParticle fluid_source("fluid_source");
    fluid_source.getParticleFromFile("data/fluid_particles.txt");
    std::cout << "源粒子数量: " << fluid_source.particle_num << "\n";
    
    // 创建目标粒子（初始为空）
    FluidParticle fluid_target("fluid_target");
    std::cout << "目标粒子初始数量: " << fluid_target.particle_num << "\n";
    
    // 执行拷贝
    std::cout << "\n执行拷贝操作...\n";
    fluid_target.copyParticle(fluid_source);
    std::cout << "拷贝后目标粒子数量: " << fluid_target.particle_num << "\n";
    
    // 验证拷贝结果
    if (fluid_target.particle_num == fluid_source.particle_num) {
        std::cout << "✓ 粒子数量匹配\n";
    } else {
        std::cout << "✗ 粒子数量不匹配\n";
    }
    
    // 测试SolidParticle拷贝
    std::cout << "\n测试 SolidParticle 拷贝:\n";
    SolidParticle solid_source("solid_source");
    solid_source.getParticleFromFile("data/solid_particles.txt");
    std::cout << "源粒子数量: " << solid_source.particle_num << "\n";
    
    SolidParticle solid_target("solid_target");
    solid_target.copyParticle(solid_source);
    std::cout << "拷贝后目标粒子数量: " << solid_target.particle_num << "\n";
    
    if (solid_target.particle_num == solid_source.particle_num) {
        std::cout << "✓ 粒子数量匹配\n";
    } else {
        std::cout << "✗ 粒子数量不匹配\n";
    }
}

void testVectorWrite() {
    printSeparator("测试向量写入功能（Debug）");
    
    FileOperator file_op;
    
    // 测试写入double向量
    std::vector<double> double_vec = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::cout << "\n1. 写入 double 向量到 TXT 文件:\n";
    bool success1 = file_op.writeVectorToFile("data/debug_double.txt", double_vec);
    std::cout << "   写入结果: " << (success1 ? "成功" : "失败") << "\n";
    
    // 测试写入double3向量
    std::vector<double3> double3_vec = {
        {1.0, 2.0, 3.0},
        {4.0, 5.0, 6.0},
        {7.0, 8.0, 9.0}
    };
    std::cout << "\n2. 写入 double3 向量到 TXT 文件:\n";
    bool success2 = file_op.writeVectorToFile("data/debug_double3.txt", double3_vec);
    std::cout << "   写入结果: " << (success2 ? "成功" : "失败") << "\n";
    
    // 测试写入double3向量到CSV文件
    std::cout << "\n3. 写入 double3 向量到 CSV 文件:\n";
    bool success3 = file_op.writeVectorToFile("data/debug_double3.csv", double3_vec);
    std::cout << "   写入结果: " << (success3 ? "成功" : "失败") << "\n";
}

void testFileTypeDetection() {
    printSeparator("测试文件类型自动检测");
    
    FileOperator file_op;
    
    // 测试不同扩展名的文件
    std::cout << "\n测试文件类型检测:\n";
    std::cout << "  data/fluid_particles.txt -> 应识别为TXT\n";
    std::cout << "  data/fluid_particles.csv -> 应识别为CSV\n";
    std::cout << "  data/solid_particles.txt -> 应识别为TXT\n";
    std::cout << "  data/solid_particles.csv -> 应识别为CSV\n";
    
    // 实际读取测试
    FluidParticle fluid1("test1");
    int count1 = fluid1.getParticleFromFile("data/fluid_particles.txt");
    std::cout << "\n  读取 .txt 文件: " << count1 << " 个粒子\n";
    
    FluidParticle fluid2("test2");
    int count2 = fluid2.getParticleFromFile("data/fluid_particles.csv");
    std::cout << "  读取 .csv 文件: " << count2 << " 个粒子\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "MPS Baseline - Particle 类功能测试\n";
    std::cout << std::string(60, '=') << "\n";
    
    try {
        // 测试文件读取功能
        testFluidParticleRead();
        testSolidParticleRead();
        
        // 测试文件写入功能
        testParticleWrite();
        
        // 测试粒子拷贝功能
        testParticleCopy();
        
        // 测试向量写入功能
        testVectorWrite();
        
        // 测试文件类型检测
        testFileTypeDetection();
        
        printSeparator("所有测试完成");
        std::cout << "\n测试结果文件已生成在 data/ 目录下:\n";
        std::cout << "  - output_fluid.txt / output_fluid.csv\n";
        std::cout << "  - output_solid.txt / output_solid.csv\n";
        std::cout << "  - debug_double.txt / debug_double3.txt / debug_double3.csv\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\n错误: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}

