#include <iostream>
#include <iomanip>
#include "../src/config/MPSConfig.hpp"
#include "ini/SimpleIni.h"

using namespace mps;

void printSeparator(const std::string& title) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "  " << title << "\n";
    std::cout << std::string(60, '=') << "\n";
}

void testSimpleIni() {
    printSeparator("测试 SimpleIni 基本功能");
    
    SimpleIni ini;
    
    // 加载配置文件
    std::cout << "\n1. 加载配置文件:\n";
    bool loaded = ini.LoadFile("config.ini");
    std::cout << "   加载结果: " << (loaded ? "成功" : "失败") << "\n";
    
    if (!loaded) {
        std::cout << "   警告: 配置文件未找到，使用默认值\n";
        return;
    }
    
    // 测试读取字符串值
    std::cout << "\n2. 读取字符串值:\n";
    std::string input_dir = ini.GetValue("File", "InputDir", "data");
    std::cout << "   InputDir: " << input_dir << "\n";
    std::string fluid_file = ini.GetValue("File", "FluidParticleFile", "");
    std::cout << "   FluidParticleFile: " << fluid_file << "\n";
    
    // 测试读取整数值
    std::cout << "\n3. 读取整数值:\n";
    int max_iter = ini.GetIntValue("Simulation", "MaxIterations", 1000);
    std::cout << "   MaxIterations: " << max_iter << "\n";
    int particle_count = ini.GetIntValue("Particle", "ParticleCount", 1000);
    std::cout << "   ParticleCount: " << particle_count << "\n";
    
    // 测试读取浮点数值
    std::cout << "\n4. 读取浮点数值:\n";
    double time_step = ini.GetDoubleValue("Simulation", "TimeStep", 0.001);
    std::cout << "   TimeStep: " << std::fixed << std::setprecision(6) << time_step << "\n";
    double density = ini.GetDoubleValue("Simulation", "Density", 1000.0);
    std::cout << "   Density: " << density << "\n";
    double gravity_z = ini.GetDoubleValue("Simulation", "GravityZ", -9.8);
    std::cout << "   GravityZ: " << gravity_z << "\n";
    
    // 测试读取布尔值
    std::cout << "\n5. 测试键是否存在:\n";
    bool has_key1 = ini.HasKey("File", "InputDir");
    std::cout << "   File.InputDir 存在: " << (has_key1 ? "是" : "否") << "\n";
    bool has_key2 = ini.HasKey("File", "NonExistentKey");
    std::cout << "   File.NonExistentKey 存在: " << (has_key2 ? "是" : "否") << "\n";
}

void testMPSConfig() {
    printSeparator("测试 MPSConfig 参数类");
    
    // 创建 MPSConfig 对象
    MPSConfig mps_config;
    
    // 从 SimpleIni 加载参数
    std::cout << "\n1. 从配置文件加载参数:\n";
    SimpleIni ini;
    bool loaded = ini.LoadFile("config.ini");
    if (loaded) {
        bool success = mps_config.LoadFromConfig(ini);
        std::cout << "   加载结果: " << (success ? "成功" : "失败") << "\n";
    } else {
        std::cout << "   警告: 配置文件未找到，使用默认配置\n";
    }
    
    // 验证参数
    std::cout << "\n2. 验证参数有效性:\n";
    bool valid = mps_config.Validate();
    std::cout << "   参数有效性: " << (valid ? "有效" : "无效") << "\n";
    
    // 显示文件配置
    std::cout << "\n3. 文件配置:\n";
    const auto& file_config = mps_config.GetFileConfig();
    std::cout << "   InputDir: " << file_config.input_dir << "\n";
    std::cout << "   OutputDir: " << file_config.output_dir << "\n";
    std::cout << "   FluidParticleFile: " << file_config.fluid_particle_file << "\n";
    std::cout << "   SolidParticleFile: " << file_config.solid_particle_file << "\n";
    
    // 显示仿真配置
    std::cout << "\n4. 仿真配置:\n";
    const auto& sim_config = mps_config.GetSimulationConfig();
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "   TimeStep: " << sim_config.time_step << "\n";
    std::cout << "   TotalTime: " << sim_config.total_time << "\n";
    std::cout << "   MaxIterations: " << sim_config.max_iterations << "\n";
    std::cout << "   Density: " << sim_config.density << "\n";
    std::cout << "   Viscosity: " << sim_config.viscosity << "\n";
    std::cout << "   Gravity: (" << sim_config.gravity_x << ", " 
              << sim_config.gravity_y << ", " << sim_config.gravity_z << ")\n";
    
    // 显示粒子配置
    std::cout << "\n5. 粒子配置:\n";
    const auto& particle_config = mps_config.GetParticleConfig();
    std::cout << "   ParticleRadius: " << particle_config.particle_radius << "\n";
    std::cout << "   SmoothingRadius: " << particle_config.smoothing_radius << "\n";
    std::cout << "   ParticleCount: " << particle_config.particle_count << "\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "MPS Baseline - 配置文件读取测试\n";
    std::cout << std::string(60, '=') << "\n";
    
    try {
        testSimpleIni();
        testMPSConfig();
        
        printSeparator("所有测试完成");
        
    } catch (const std::exception& e) {
        std::cerr << "\n错误: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}

