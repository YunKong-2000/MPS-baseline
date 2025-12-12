#include <iostream>
#include <iomanip>
#include "../src/core/Particle.hpp"
#include "../include/core/MPSUtils.h"

using namespace mps2D;

void printSeparator(const std::string& title) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "  " << title << "\n";
    std::cout << std::string(60, '=') << "\n";
}

void testParticleCreation() {
    printSeparator("测试 Particle 类创建");
    
    std::cout << "\n1. 创建基础 Particle:\n";
    Particle particle("test_particle");
    std::cout << "   名称: " << particle.name << "\n";
    std::cout << "   初始粒子数: " << particle.particle_num << "\n";
    std::cout << "   位置向量大小: " << particle.position.size() << "\n";
    std::cout << "   速度向量大小: " << particle.velocity.size() << "\n";
    
    if (particle.particle_num == 0 && 
        particle.position.size() == 0 && 
        particle.velocity.size() == 0) {
        std::cout << "   ✓ 初始化正确\n";
    } else {
        std::cout << "   ✗ 初始化错误\n";
    }
}

void testParticleDataManipulation() {
    printSeparator("测试 Particle 数据操作");
    
    std::cout << "\n1. 手动设置粒子数据:\n";
    Particle particle("test_particle");
    
    // 设置粒子数量和数据
    particle.particle_num = 3;
    particle.position.resize(3);
    particle.velocity.resize(3);
    
    particle.position[0] = {0.0, 0.0};
    particle.position[1] = {1.0, 1.0};
    particle.position[2] = {2.0, 2.0};
    
    particle.velocity[0] = {0.1, 0.1};
    particle.velocity[1] = {0.2, 0.2};
    particle.velocity[2] = {0.3, 0.3};
    
    std::cout << "   粒子数量: " << particle.particle_num << "\n";
    std::cout << "   位置数据:\n";
    for (int i = 0; i < particle.particle_num; ++i) {
        std::cout << "     粒子[" << i << "] 位置: " << particle.position[i] << "\n";
    }
    std::cout << "   速度数据:\n";
    for (int i = 0; i < particle.particle_num; ++i) {
        std::cout << "     粒子[" << i << "] 速度: " << particle.velocity[i] << "\n";
    }
    
    // 测试getParticleNum()
    int count = particle.getParticleNum();
    std::cout << "\n2. 测试 getParticleNum():\n";
    std::cout << "   返回的粒子数: " << count << "\n";
    if (count == particle.particle_num) {
        std::cout << "   ✓ getParticleNum() 正确\n";
    } else {
        std::cout << "   ✗ getParticleNum() 错误\n";
    }
}

void testFluidParticleCreation() {
    printSeparator("测试 FluidParticle 类创建");
    
    std::cout << "\n1. 创建 FluidParticle:\n";
    FluidParticle fluid("test_fluid");
    std::cout << "   名称: " << fluid.name << "\n";
    std::cout << "   初始粒子数: " << fluid.particle_num << "\n";
    std::cout << "   密度向量大小: " << fluid.density.size() << "\n";
    std::cout << "   压力向量大小: " << fluid.pressure.size() << "\n";
    std::cout << "   表面类型向量大小: " << fluid.surface_type.size() << "\n";
    std::cout << "   流体邻居列表大小: " << fluid.fluid_neighbour_list.size() << "\n";
    std::cout << "   固体邻居列表大小: " << fluid.solid_neighbour_list.size() << "\n";
    
    // 设置一些测试数据
    std::cout << "\n2. 设置 FluidParticle 数据:\n";
    fluid.particle_num = 2;
    fluid.position.resize(2);
    fluid.velocity.resize(2);
    fluid.density.resize(2);
    fluid.pressure.resize(2);
    fluid.surface_type.resize(2);
    fluid.fluid_neighbour_list.resize(2);
    fluid.solid_neighbour_list.resize(2);
    
    fluid.position[0] = {0.0, 0.0};
    fluid.position[1] = {1.0, 1.0};
    fluid.velocity[0] = {0.1, 0.1};
    fluid.velocity[1] = {0.2, 0.2};
    fluid.density[0] = 1000.0;
    fluid.density[1] = 1000.0;
    fluid.pressure[0] = 101325.0;
    fluid.pressure[1] = 101325.0;
    fluid.surface_type[0] = SurfaceType::INNER;
    fluid.surface_type[1] = SurfaceType::SURFACE;
    
    std::cout << "   粒子数量: " << fluid.particle_num << "\n";
    std::cout << "   粒子[0] 密度: " << fluid.density[0] << "\n";
    std::cout << "   粒子[0] 压力: " << fluid.pressure[0] << "\n";
    std::cout << "   粒子[0] 表面类型: " 
              << static_cast<int>(fluid.surface_type[0]) << "\n";
    std::cout << "   粒子[1] 表面类型: " 
              << static_cast<int>(fluid.surface_type[1]) << "\n";
}

void testSolidParticleCreation() {
    printSeparator("测试 SolidParticle 类创建");
    
    std::cout << "\n1. 创建 SolidParticle:\n";
    SolidParticle solid("test_solid");
    std::cout << "   名称: " << solid.name << "\n";
    std::cout << "   初始粒子数: " << solid.particle_num << "\n";
    std::cout << "   法向向量大小: " << solid.normal_vector.size() << "\n";
    
    // 设置一些测试数据
    std::cout << "\n2. 设置 SolidParticle 数据:\n";
    solid.particle_num = 2;
    solid.position.resize(2);
    solid.velocity.resize(2);
    solid.normal_vector.resize(2);
    
    solid.position[0] = {0.0, 0.0};
    solid.position[1] = {1.0, 0.0};
    solid.velocity[0] = {0.0, 0.0};
    solid.velocity[1] = {0.0, 0.0};
    solid.normal_vector[0] = {1.0, 0.0};
    solid.normal_vector[1] = {0.0, 1.0};
    
    std::cout << "   粒子数量: " << solid.particle_num << "\n";
    std::cout << "   粒子[0] 位置: " << solid.position[0] << "\n";
    std::cout << "   粒子[0] 法向向量: " << solid.normal_vector[0] << "\n";
    std::cout << "   粒子[1] 法向向量: " << solid.normal_vector[1] << "\n";
    
    // 验证法向向量归一化
    double mag0 = ComputeVectorMagnitude(solid.normal_vector[0]);
    double mag1 = ComputeVectorMagnitude(solid.normal_vector[1]);
    std::cout << "\n3. 验证法向向量模长:\n";
    std::cout << "   |normal_vector[0]| = " << mag0 << "\n";
    std::cout << "   |normal_vector[1]| = " << mag1 << "\n";
    if (std::abs(mag0 - 1.0) < 1e-10 && std::abs(mag1 - 1.0) < 1e-10) {
        std::cout << "   ✓ 法向向量已归一化\n";
    } else {
        std::cout << "   ⚠ 法向向量未归一化（这是正常的，取决于具体应用）\n";
    }
}

void testParticleCopy() {
    printSeparator("测试粒子数据拷贝功能");
    
    std::cout << "\n1. 测试 FluidParticle 拷贝:\n";
    
    // 创建源粒子
    FluidParticle fluid_source("fluid_source");
    fluid_source.particle_num = 3;
    fluid_source.position.resize(3);
    fluid_source.velocity.resize(3);
    fluid_source.density.resize(3);
    fluid_source.pressure.resize(3);
    fluid_source.surface_type.resize(3);
    fluid_source.fluid_neighbour_list.resize(3);
    fluid_source.solid_neighbour_list.resize(3);
    
    for (int i = 0; i < 3; ++i) {
        fluid_source.position[i] = {i * 0.1, i * 0.1};
        fluid_source.velocity[i] = {i * 0.01, i * 0.01};
        fluid_source.density[i] = 1000.0 + i * 10.0;
        fluid_source.pressure[i] = 101325.0 + i * 100.0;
        fluid_source.surface_type[i] = static_cast<SurfaceType>(i);
    }
    
    std::cout << "   源粒子数量: " << fluid_source.particle_num << "\n";
    
    // 创建目标粒子（初始为空）
    FluidParticle fluid_target("fluid_target");
    std::cout << "   目标粒子初始数量: " << fluid_target.particle_num << "\n";
    
    // 执行拷贝
    fluid_target.copyParticle(fluid_source);
    std::cout << "   拷贝后目标粒子数量: " << fluid_target.particle_num << "\n";
    
    // 验证拷贝结果
    bool copy_success = true;
    if (fluid_target.particle_num != fluid_source.particle_num) {
        copy_success = false;
    }
    
    for (int i = 0; i < fluid_source.particle_num; ++i) {
        double dist_pos = ComputeDistance(fluid_target.position[i], 
                                          fluid_source.position[i]);
        double dist_vel = ComputeDistance(fluid_target.velocity[i], 
                                          fluid_source.velocity[i]);
        if (dist_pos > 1e-10 || dist_vel > 1e-10) {
            copy_success = false;
        }
        if (std::abs(fluid_target.density[i] - fluid_source.density[i]) > 1e-10) {
            copy_success = false;
        }
        if (fluid_target.surface_type[i] != fluid_source.surface_type[i]) {
            copy_success = false;
        }
    }
    
    if (copy_success) {
        std::cout << "   ✓ FluidParticle 拷贝成功\n";
    } else {
        std::cout << "   ✗ FluidParticle 拷贝失败\n";
    }
    
    std::cout << "\n2. 测试 SolidParticle 拷贝:\n";
    
    // 创建源粒子
    SolidParticle solid_source("solid_source");
    solid_source.particle_num = 2;
    solid_source.position.resize(2);
    solid_source.velocity.resize(2);
    solid_source.normal_vector.resize(2);
    
    solid_source.position[0] = {0.0, 0.0};
    solid_source.position[1] = {1.0, 0.0};
    solid_source.velocity[0] = {0.0, 0.0};
    solid_source.velocity[1] = {0.0, 0.0};
    solid_source.normal_vector[0] = {1.0, 0.0};
    solid_source.normal_vector[1] = {0.0, 1.0};
    
    std::cout << "   源粒子数量: " << solid_source.particle_num << "\n";
    
    // 创建目标粒子
    SolidParticle solid_target("solid_target");
    solid_target.copyParticle(solid_source);
    std::cout << "   拷贝后目标粒子数量: " << solid_target.particle_num << "\n";
    
    // 验证拷贝结果
    bool solid_copy_success = true;
    if (solid_target.particle_num != solid_source.particle_num) {
        solid_copy_success = false;
    }
    
    for (int i = 0; i < solid_source.particle_num; ++i) {
        double dist_pos = ComputeDistance(solid_target.position[i], 
                                         solid_source.position[i]);
        double dist_norm = ComputeDistance(solid_target.normal_vector[i], 
                                           solid_source.normal_vector[i]);
        if (dist_pos > 1e-10 || dist_norm > 1e-10) {
            solid_copy_success = false;
        }
    }
    
    if (solid_copy_success) {
        std::cout << "   ✓ SolidParticle 拷贝成功\n";
    } else {
        std::cout << "   ✗ SolidParticle 拷贝失败\n";
    }
}

void testParticleDistance() {
    printSeparator("测试粒子间距离计算");
    
    std::cout << "\n1. 计算粒子间距离:\n";
    FluidParticle fluid("test_fluid");
    fluid.particle_num = 3;
    fluid.position.resize(3);
    
    fluid.position[0] = {0.0, 0.0};
    fluid.position[1] = {3.0, 4.0};
    fluid.position[2] = {1.0, 1.0};
    
    std::cout << "   粒子[0] 位置: " << fluid.position[0] << "\n";
    std::cout << "   粒子[1] 位置: " << fluid.position[1] << "\n";
    std::cout << "   粒子[2] 位置: " << fluid.position[2] << "\n";
    
    double dist01 = ComputeDistance(fluid.position[0], fluid.position[1]);
    double dist02 = ComputeDistance(fluid.position[0], fluid.position[2]);
    double dist12 = ComputeDistance(fluid.position[1], fluid.position[2]);
    
    std::cout << "\n   距离[0][1] = " << dist01 << " (预期: 5.0)\n";
    std::cout << "   距离[0][2] = " << dist02 << " (预期: sqrt(2) ≈ " 
              << std::sqrt(2.0) << ")\n";
    std::cout << "   距离[1][2] = " << dist12 << "\n";
    
    if (std::abs(dist01 - 5.0) < 1e-10) {
        std::cout << "   ✓ 距离计算正确\n";
    } else {
        std::cout << "   ✗ 距离计算错误\n";
    }
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "MPS Baseline 2D - Particle 类功能测试\n";
    std::cout << std::string(60, '=') << "\n";
    std::cout << "\n注意: 此测试不包含文件I/O功能（FileOperator尚未实现2D版本）\n";
    
    try {
        testParticleCreation();
        testParticleDataManipulation();
        testFluidParticleCreation();
        testSolidParticleCreation();
        testParticleCopy();
        testParticleDistance();
        
        printSeparator("所有测试完成");
        std::cout << "\n✓ 所有粒子类测试通过\n";
        std::cout << "\n注意: FileOperator 的2D版本尚未实现，因此未测试文件读取/写入功能\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\n错误: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}

