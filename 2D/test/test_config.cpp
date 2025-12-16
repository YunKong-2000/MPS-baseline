#include <iostream>
#include <iomanip>
#include "../src/config/MPSConfig2D.hpp"
#include "../../third_party/ini/SimpleIni.h"
#include "core/ErrorHandling.h"

using namespace mps2D;

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
    
    // 测试读取浮点数值
    std::cout << "\n4. 读取浮点数值:\n";
    double time_step = ini.GetDoubleValue("Simulation", "TimeStep", 0.001);
    std::cout << "   TimeStep: " << std::fixed << std::setprecision(6) << time_step << "\n";
    double density = ini.GetDoubleValue("Simulation", "Density", 1000.0);
    std::cout << "   Density: " << density << "\n";
    double gravity_y = ini.GetDoubleValue("Simulation", "GravityY", -9.8);
    std::cout << "   GravityY: " << gravity_y << "\n";
}

void testMPSConfig2D() {
    printSeparator("测试 MPSConfig2D 配置类");
    
    MPSConfig2D config;
    SimpleIni ini;
    
    // 加载配置文件
    std::cout << "\n1. 加载配置文件:\n";
    bool loaded = ini.LoadFile("config.ini");
    std::cout << "   加载结果: " << (loaded ? "成功" : "失败") << "\n";
    
    if (!loaded) {
        std::cout << "   警告: 配置文件未找到，使用默认值\n";
    }
    
    // 从配置文件加载参数
    std::cout << "\n2. 从配置文件加载参数:\n";
    bool success = config.LoadFromConfig(ini);
    std::cout << "   加载结果: " << (success ? "成功" : "失败") << "\n";
    
    if (!success) {
        std::cout << "   错误: 参数验证失败\n";
        return;
    }
    
    // 显示文件配置
    std::cout << "\n3. 文件配置:\n";
    const auto& file_config = config.GetFileConfig();
    std::cout << "   InputDir: " << file_config.input_dir << "\n";
    std::cout << "   OutputDir: " << file_config.output_dir << "\n";
    std::cout << "   FluidParticleFile: " << file_config.fluid_particle_file << "\n";
    std::cout << "   SolidParticleFile: " << file_config.solid_particle_file << "\n";
    std::cout << "   ResultFile: " << file_config.result_file << "\n";
    
    // 显示仿真配置
    std::cout << "\n4. 仿真配置:\n";
    const auto& sim_config = config.GetSimulationConfig();
    std::cout << "   TimeStep: " << std::fixed << std::setprecision(6) << sim_config.time_step << " s\n";
    std::cout << "   TotalTime: " << sim_config.total_time << " s\n";
    std::cout << "   MaxIterations: " << sim_config.max_iterations << "\n";
    std::cout << "   Density: " << sim_config.density << " kg/m³\n";
    std::cout << "   KinematicViscosity: " << sim_config.kinematic_viscosity << " m²/s\n";
    std::cout << "   Gravity: (" << sim_config.gravity_x << ", " 
              << sim_config.gravity_y << ") m/s²\n";
    
    // 显示粒子配置
    std::cout << "\n5. 粒子配置:\n";
    const auto& particle_config = config.GetParticleConfig();
    std::cout << "   ParticleSpacing: " << particle_config.particle_spacing << " m\n";
    std::cout << "   ParticleRadius: " << particle_config.particle_radius << " m\n";
    std::cout << "   SmoothingRadius: " << particle_config.smoothing_radius << " m\n";
    std::cout << "   CellSize: " << particle_config.cell_size << " m\n";
}

void testDefaultValues() {
    printSeparator("测试默认值");
    
    MPSConfig2D config;
    
    std::cout << "\n使用默认配置（不加载文件）:\n";
    
    const auto& file_config = config.GetFileConfig();
    std::cout << "   InputDir: " << file_config.input_dir << "\n";
    
    const auto& sim_config = config.GetSimulationConfig();
    std::cout << "   TimeStep: " << sim_config.time_step << " s\n";
    std::cout << "   GravityY: " << sim_config.gravity_y << " m/s²\n";
    
    const auto& particle_config = config.GetParticleConfig();
    std::cout << "   ParticleSpacing: " << particle_config.particle_spacing << " m\n";
}

void testParameterValidation() {
    printSeparator("测试参数验证（使用CHECK宏）");
    
    MPSConfig2D config;
    
    std::cout << "\n测试无效参数值:\n";
    
    // 测试无效的时间步长
    std::cout << "\n1. 测试无效的时间步长 (TimeStep = -0.001):\n";
    auto sim_config = config.GetSimulationConfig();
    sim_config.time_step = -0.001;
    config.SetSimulationConfig(sim_config);
    bool valid = config.Validate();
    std::cout << "   验证结果: " << (valid ? "通过" : "失败") << "\n";
    std::cout << "   预期: 失败（时间步长必须为正数）\n";
    
    // 测试无效的密度
    std::cout << "\n2. 测试无效的密度 (Density = 0):\n";
    sim_config = config.GetSimulationConfig();
    sim_config.time_step = 0.001;  // 恢复有效值
    sim_config.density = 0.0;
    config.SetSimulationConfig(sim_config);
    valid = config.Validate();
    std::cout << "   验证结果: " << (valid ? "通过" : "失败") << "\n";
    std::cout << "   预期: 失败（密度必须为正数）\n";
    
    // 测试无效的动力学粘性系数
    std::cout << "\n3. 测试无效的动力学粘性系数 (KinematicViscosity = -0.001):\n";
    sim_config = config.GetSimulationConfig();
    sim_config.density = 1000.0;  // 恢复有效值
    sim_config.kinematic_viscosity = -0.001;
    config.SetSimulationConfig(sim_config);
    valid = config.Validate();
    std::cout << "   验证结果: " << (valid ? "通过" : "失败") << "\n";
    std::cout << "   预期: 失败（动力学粘性系数必须为非负数）\n";
    
    // 测试无效的粒子间距
    std::cout << "\n4. 测试无效的粒子间距 (ParticleSpacing = 0):\n";
    auto particle_config = config.GetParticleConfig();
    particle_config.particle_spacing = 0.0;
    config.SetParticleConfig(particle_config);
    valid = config.Validate();
    std::cout << "   验证结果: " << (valid ? "通过" : "失败") << "\n";
    std::cout << "   预期: 失败（粒子间距必须为正数）\n";
    
    // 测试所有参数都有效
    std::cout << "\n5. 测试所有参数都有效:\n";
    sim_config = config.GetSimulationConfig();
    sim_config.kinematic_viscosity = 0.001;  // 恢复有效值
    config.SetSimulationConfig(sim_config);
    particle_config = config.GetParticleConfig();
    particle_config.particle_spacing = 0.02;  // 恢复有效值
    config.SetParticleConfig(particle_config);
    valid = config.Validate();
    std::cout << "   验证结果: " << (valid ? "通过" : "失败") << "\n";
    std::cout << "   预期: 通过\n";
    
    // 测试CHECK宏直接抛出异常
    std::cout << "\n6. 测试CHECK宏直接抛出异常:\n";
    try {
        CHECK_POSITIVE(-0.001, "测试值");
        std::cout << "   错误: 应该抛出异常但没有抛出\n";
    } catch (const MPSException& e) {
        std::cout << "   捕获到异常: " << e.what() << "\n";
        std::cout << "   预期: 包含'测试值'和'必须为正数'\n";
    }
    
    // 测试CHECK_NON_NEGATIVE
    std::cout << "\n7. 测试CHECK_NON_NEGATIVE宏:\n";
    try {
        CHECK_NON_NEGATIVE(-0.001, "测试值");
        std::cout << "   错误: 应该抛出异常但没有抛出\n";
    } catch (const MPSException& e) {
        std::cout << "   捕获到异常: " << e.what() << "\n";
        std::cout << "   预期: 包含'测试值'和'必须为非负数'\n";
    }
}

int main() {
    std::cout << "=== MPS 2D 配置模块测试 ===\n";
    
    testSimpleIni();
    testMPSConfig2D();
    testDefaultValues();
    testParameterValidation();
    
    std::cout << "\n测试完成!\n";
    return 0;
}
