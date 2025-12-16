#include "core/TimeStepManager.h"
#include "../src/config/MPSConfig2D.hpp"
#include "../src/core/Particle.hpp"
#include "core/MPSUtils.h"
#include "core/ErrorHandling.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>

using namespace mps2D;

void printSeparator(const std::string& title) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "  " << title << "\n";
    std::cout << std::string(60, '=') << "\n";
}

// 创建测试用的配置对象
MPSConfig2D CreateTestConfig() {
    MPSConfig2D config;
    auto sim_config = config.GetSimulationConfig();
    
    // 设置时间步管理参数
    sim_config.time_step = 0.001;        // 初始时间步
    sim_config.total_time = 1.0;         // 总时间
    sim_config.min_time_step = 1e-6;    // 最小时间步
    sim_config.max_time_step = 0.01;     // 最大时间步
    sim_config.max_cfl = 0.5;            // 最大CFL数
    sim_config.output_interval = 0.1;    // 输出间隔
    
    config.SetSimulationConfig(sim_config);
    return config;
}

// 创建测试用的流体粒子
FluidParticle CreateTestFluidParticles(int num_particles, double max_velocity) {
    FluidParticle particles("test_fluid");
    particles.particle_num = num_particles;
    particles.position.resize(num_particles);
    particles.velocity.resize(num_particles);
    particles.density.resize(num_particles);
    particles.pressure.resize(num_particles);
    particles.surface_type.resize(num_particles);
    particles.fluid_neighbour_list.resize(num_particles);
    particles.solid_neighbour_list.resize(num_particles);
    
    // 设置粒子位置和速度
    for (int i = 0; i < num_particles; ++i) {
        particles.position[i] = {static_cast<double>(i) * 0.02, 0.0};
        // 设置速度，使最大速度为max_velocity
        double velocity_magnitude = max_velocity * (1.0 - static_cast<double>(i) / num_particles);
        particles.velocity[i] = {velocity_magnitude, 0.0};
        particles.density[i] = 1000.0;
        particles.pressure[i] = 0.0;
        particles.surface_type[i] = SurfaceType::INNER;
    }
    
    return particles;
}

void testBasicFunctionality() {
    printSeparator("测试基本功能");
    
    MPSConfig2D config = CreateTestConfig();
    double particle_spacing = 0.02;
    
    std::cout << "\n1. 创建TimeStepManager:\n";
    TimeStepManager manager(config, particle_spacing);
    
    std::cout << "   初始时间步: " << std::fixed << std::setprecision(6) 
              << manager.GetTimeStep() << " s\n";
    std::cout << "   当前时间: " << manager.GetCurrentTime() << " s\n";
    std::cout << "   总时间: " << manager.GetTotalTime() << " s\n";
    std::cout << "   模拟是否完成: " << (manager.IsSimulationFinished() ? "是" : "否") << "\n";
    
    // 验证初始值
    assert(std::abs(manager.GetTimeStep() - 0.001) < 1e-10);
    assert(std::abs(manager.GetCurrentTime() - 0.0) < 1e-10);
    assert(std::abs(manager.GetTotalTime() - 1.0) < 1e-10);
    assert(!manager.IsSimulationFinished());
    
    std::cout << "   ✓ 基本功能测试通过\n";
}

void testTimeAdvancement() {
    printSeparator("测试时间推进");
    
    MPSConfig2D config = CreateTestConfig();
    TimeStepManager manager(config, 0.02);
    
    std::cout << "\n1. 推进时间:\n";
    double initial_time = manager.GetCurrentTime();
    double time_step = manager.GetTimeStep();
    
    manager.AdvanceTime();
    double new_time = manager.GetCurrentTime();
    
    std::cout << "   初始时间: " << initial_time << " s\n";
    std::cout << "   时间步: " << time_step << " s\n";
    std::cout << "   推进后时间: " << new_time << " s\n";
    
    assert(std::abs(new_time - (initial_time + time_step)) < 1e-10);
    std::cout << "   ✓ 时间推进测试通过\n";
    
    std::cout << "\n2. 多次推进时间:\n";
    int steps = 0;
    while (!manager.IsSimulationFinished() && steps < 1000) {
        manager.AdvanceTime();
        steps++;
    }
    
    std::cout << "   推进步数: " << steps << "\n";
    std::cout << "   最终时间: " << manager.GetCurrentTime() << " s\n";
    std::cout << "   模拟完成: " << (manager.IsSimulationFinished() ? "是" : "否") << "\n";
    
    assert(manager.IsSimulationFinished() || manager.GetCurrentTime() >= manager.GetTotalTime());
    std::cout << "   ✓ 多次推进测试通过\n";
}

void testOutputInterval() {
    printSeparator("测试输出文件间隔");
    
    MPSConfig2D config = CreateTestConfig();
    TimeStepManager manager(config, 0.02);
    
    std::cout << "\n1. 初始状态:\n";
    std::cout << "   当前时间: " << manager.GetCurrentTime() << " s\n";
    std::cout << "   输出间隔: 0.1 s\n";
    std::cout << "   是否需要输出: " << (manager.ShouldOutput() ? "是" : "否") << "\n";
    
    assert(!manager.ShouldOutput());  // 初始时不应该输出
    std::cout << "   ✓ 初始状态测试通过\n";
    
    std::cout << "\n2. 推进到输出时间:\n";
    int output_count = 0;
    double current_time = 0.0;
    
    while (!manager.IsSimulationFinished() && output_count < 5) {
        // 调整时间步（使用零速度粒子，时间步不会变化太大）
        FluidParticle particles = CreateTestFluidParticles(10, 0.0);
        manager.AdjustTimeStep(particles);
        
        manager.AdvanceTime();
        current_time = manager.GetCurrentTime();
        
        if (manager.ShouldOutput()) {
            output_count++;
            std::cout << "   输出 #" << output_count << " 在时间: " 
                      << std::fixed << std::setprecision(6) << current_time << " s\n";
            manager.MarkOutputDone();
        }
    }
    
    assert(output_count > 0);
    std::cout << "   ✓ 输出间隔测试通过\n";
}

void testTimeStepAdjustment() {
    printSeparator("测试时间步调整");
    
    MPSConfig2D config = CreateTestConfig();
    double particle_spacing = 0.02;
    TimeStepManager manager(config, particle_spacing);
    
    std::cout << "\n1. 低速度情况（应该增加时间步）:\n";
    FluidParticle low_velocity_particles = CreateTestFluidParticles(10, 0.1);  // 最大速度0.1 m/s
    double initial_dt = manager.GetTimeStep();
    
    double adjusted_dt = manager.AdjustTimeStep(low_velocity_particles);
    std::cout << "   初始时间步: " << initial_dt << " s\n";
    std::cout << "   调整后时间步: " << adjusted_dt << " s\n";
    
    // 计算CFL
    double max_velocity = 0.1;
    double cfl = manager.GetCFL(max_velocity);
    std::cout << "   当前CFL数: " << cfl << "\n";
    
    // 低速度时，时间步应该可以增加（在限制范围内）
    assert(adjusted_dt >= initial_dt * 0.8);  // 至少不会大幅减少
    std::cout << "   ✓ 低速度测试通过\n";
    
    std::cout << "\n2. 高速度情况（应该减少时间步）:\n";
    FluidParticle high_velocity_particles = CreateTestFluidParticles(10, 10.0);  // 最大速度10 m/s
    manager.Reset();  // 重置到初始状态
    initial_dt = manager.GetTimeStep();
    
    adjusted_dt = manager.AdjustTimeStep(high_velocity_particles);
    std::cout << "   初始时间步: " << initial_dt << " s\n";
    std::cout << "   调整后时间步: " << adjusted_dt << " s\n";
    
    max_velocity = 10.0;
    cfl = manager.GetCFL(max_velocity);
    std::cout << "   当前CFL数: " << cfl << "\n";
    
    // 高速度时，时间步应该减少以满足CFL条件
    assert(adjusted_dt <= initial_dt * 1.2);  // 不会大幅增加
    assert(cfl <= config.GetSimulationConfig().max_cfl + 0.1);  // CFL应该在合理范围内
    std::cout << "   ✓ 高速度测试通过\n";
    
    std::cout << "\n3. 零速度情况:\n";
    FluidParticle zero_velocity_particles = CreateTestFluidParticles(10, 0.0);
    manager.Reset();
    initial_dt = manager.GetTimeStep();
    
    adjusted_dt = manager.AdjustTimeStep(zero_velocity_particles);
    std::cout << "   初始时间步: " << initial_dt << " s\n";
    std::cout << "   调整后时间步: " << adjusted_dt << " s\n";
    
    // 零速度时，应该使用最大时间步
    assert(adjusted_dt <= config.GetSimulationConfig().max_time_step + 1e-10);
    std::cout << "   ✓ 零速度测试通过\n";
}

void testCFLCalculation() {
    printSeparator("测试CFL计算");
    
    MPSConfig2D config = CreateTestConfig();
    double particle_spacing = 0.02;
    TimeStepManager manager(config, particle_spacing);
    
    std::cout << "\n1. CFL计算公式验证:\n";
    double max_velocity = 5.0;  // m/s
    double time_step = manager.GetTimeStep();
    double expected_cfl = (max_velocity * time_step) / particle_spacing;
    double actual_cfl = manager.GetCFL(max_velocity);
    
    std::cout << "   最大速度: " << max_velocity << " m/s\n";
    std::cout << "   时间步: " << time_step << " s\n";
    std::cout << "   粒子间距: " << particle_spacing << " m\n";
    std::cout << "   预期CFL: " << expected_cfl << "\n";
    std::cout << "   实际CFL: " << actual_cfl << "\n";
    
    assert(std::abs(actual_cfl - expected_cfl) < 1e-10);
    std::cout << "   ✓ CFL计算测试通过\n";
    
    std::cout << "\n2. 零速度CFL:\n";
    double zero_cfl = manager.GetCFL(0.0);
    std::cout << "   零速度CFL: " << zero_cfl << "\n";
    assert(std::abs(zero_cfl) < 1e-10);
    std::cout << "   ✓ 零速度CFL测试通过\n";
}

void testReset() {
    printSeparator("测试重置功能");
    
    MPSConfig2D config = CreateTestConfig();
    TimeStepManager manager(config, 0.02);
    
    std::cout << "\n1. 初始状态:\n";
    double initial_time = manager.GetCurrentTime();
    double initial_dt = manager.GetTimeStep();
    std::cout << "   初始时间: " << initial_time << " s\n";
    std::cout << "   初始时间步: " << initial_dt << " s\n";
    
    std::cout << "\n2. 推进和调整:\n";
    // 推进时间
    manager.AdvanceTime();
    manager.AdvanceTime();
    
    // 调整时间步
    FluidParticle particles = CreateTestFluidParticles(10, 5.0);
    manager.AdjustTimeStep(particles);
    
    double modified_time = manager.GetCurrentTime();
    double modified_dt = manager.GetTimeStep();
    std::cout << "   修改后时间: " << modified_time << " s\n";
    std::cout << "   修改后时间步: " << modified_dt << " s\n";
    
    assert(modified_time > initial_time);
    
    std::cout << "\n3. 重置:\n";
    manager.Reset();
    double reset_time = manager.GetCurrentTime();
    double reset_dt = manager.GetTimeStep();
    std::cout << "   重置后时间: " << reset_time << " s\n";
    std::cout << "   重置后时间步: " << reset_dt << " s\n";
    
    assert(std::abs(reset_time - initial_time) < 1e-10);
    assert(std::abs(reset_dt - initial_dt) < 1e-10);
    std::cout << "   ✓ 重置功能测试通过\n";
}

void testTimeStepLimits() {
    printSeparator("测试时间步限制");
    
    MPSConfig2D config = CreateTestConfig();
    double particle_spacing = 0.02;
    TimeStepManager manager(config, particle_spacing);
    
    std::cout << "\n1. 测试最小时间步限制:\n";
    // 使用极高的速度，应该触发最小时间步限制
    FluidParticle extreme_velocity_particles = CreateTestFluidParticles(10, 1000.0);  // 1000 m/s
    double adjusted_dt = manager.AdjustTimeStep(extreme_velocity_particles);
    
    std::cout << "   调整后时间步: " << adjusted_dt << " s\n";
    std::cout << "   最小时间步限制: " << config.GetSimulationConfig().min_time_step << " s\n";
    
    assert(adjusted_dt >= config.GetSimulationConfig().min_time_step - 1e-10);
    std::cout << "   ✓ 最小时间步限制测试通过\n";
    
    std::cout << "\n2. 测试最大时间步限制:\n";
    manager.Reset();
    // 使用零速度，应该使用最大时间步
    FluidParticle zero_particles = CreateTestFluidParticles(10, 0.0);
    adjusted_dt = manager.AdjustTimeStep(zero_particles);
    
    std::cout << "   调整后时间步: " << adjusted_dt << " s\n";
    std::cout << "   最大时间步限制: " << config.GetSimulationConfig().max_time_step << " s\n";
    
    assert(adjusted_dt <= config.GetSimulationConfig().max_time_step + 1e-10);
    std::cout << "   ✓ 最大时间步限制测试通过\n";
}

void testSimulationCompletion() {
    printSeparator("测试模拟完成判断");
    
    MPSConfig2D config = CreateTestConfig();
    TimeStepManager manager(config, 0.02);
    
    std::cout << "\n1. 初始状态:\n";
    std::cout << "   当前时间: " << manager.GetCurrentTime() << " s\n";
    std::cout << "   总时间: " << manager.GetTotalTime() << " s\n";
    std::cout << "   模拟完成: " << (manager.IsSimulationFinished() ? "是" : "否") << "\n";
    
    assert(!manager.IsSimulationFinished());
    
    std::cout << "\n2. 推进到完成:\n";
    int steps = 0;
    while (!manager.IsSimulationFinished() && steps < 2000) {
        FluidParticle particles = CreateTestFluidParticles(10, 1.0);
        manager.AdjustTimeStep(particles);
        manager.AdvanceTime();
        steps++;
        
        if (steps % 100 == 0) {
            std::cout << "   步数: " << steps 
                      << ", 时间: " << std::fixed << std::setprecision(6) 
                      << manager.GetCurrentTime() << " s\n";
        }
    }
    
    std::cout << "   总步数: " << steps << "\n";
    std::cout << "   最终时间: " << manager.GetCurrentTime() << " s\n";
    std::cout << "   模拟完成: " << (manager.IsSimulationFinished() ? "是" : "否") << "\n";
    
    assert(manager.IsSimulationFinished());
    assert(manager.GetCurrentTime() >= manager.GetTotalTime() - 1e-6);
    std::cout << "   ✓ 模拟完成判断测试通过\n";
}

int main() {
    std::cout << "=== MPS 2D 时间步管理模块测试 ===\n";
    
    try {
        testBasicFunctionality();
        testTimeAdvancement();
        testOutputInterval();
        testTimeStepAdjustment();
        testCFLCalculation();
        testReset();
        testTimeStepLimits();
        testSimulationCompletion();
        
        std::cout << "\n" << std::string(60, '=') << "\n";
        std::cout << "  所有测试通过！\n";
        std::cout << std::string(60, '=') << "\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\n错误: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
