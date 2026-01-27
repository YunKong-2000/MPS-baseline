#include "MPSConfig2D.hpp"
#include "core/ErrorHandling.h"
#include <iostream>

namespace mps2D {

bool MPSConfig2D::LoadFromConfig(const SimpleIni& config) {
    // 加载文件配置
    file_config_.input_dir = config.GetValue("File", "InputDir", file_config_.input_dir);
    file_config_.output_dir = config.GetValue("File", "OutputDir", file_config_.output_dir);
    file_config_.fluid_particle_file = config.GetValue("File", "FluidParticleFile", 
                                                       file_config_.fluid_particle_file);
    file_config_.solid_particle_file = config.GetValue("File", "SolidParticleFile", 
                                                       file_config_.solid_particle_file);
    file_config_.result_file = config.GetValue("File", "ResultFile", file_config_.result_file);

    // 加载仿真配置
    simulation_config_.time_step = config.GetDoubleValue("Simulation", "TimeStep", 
                                                         simulation_config_.time_step);
    simulation_config_.total_time = config.GetDoubleValue("Simulation", "TotalTime", 
                                                          simulation_config_.total_time);
    simulation_config_.max_iterations = config.GetIntValue("Simulation", "MaxIterations", 
                                                           simulation_config_.max_iterations);
    simulation_config_.density = config.GetDoubleValue("Simulation", "Density", 
                                                       simulation_config_.density);
    simulation_config_.kinematic_viscosity = config.GetDoubleValue("Simulation", "KinematicViscosity", 
                                                                   simulation_config_.kinematic_viscosity);
    simulation_config_.gravity_x = config.GetDoubleValue("Simulation", "GravityX", 
                                                         simulation_config_.gravity_x);
    simulation_config_.gravity_y = config.GetDoubleValue("Simulation", "GravityY", 
                                                         simulation_config_.gravity_y);
    
    // 加载时间步管理参数
    simulation_config_.min_time_step = config.GetDoubleValue("Simulation", "MinTimeStep", 
                                                             simulation_config_.min_time_step);
    simulation_config_.max_time_step = config.GetDoubleValue("Simulation", "MaxTimeStep", 
                                                             simulation_config_.max_time_step);
    simulation_config_.max_cfl = config.GetDoubleValue("Simulation", "MaxCFL", 
                                                       simulation_config_.max_cfl);
    simulation_config_.output_interval = config.GetDoubleValue("Simulation", "OutputInterval", 
                                                               simulation_config_.output_interval);

    // 加载粒子配置
    particle_config_.particle_spacing = config.GetDoubleValue("Particle", "ParticleSpacing", 
                                                              particle_config_.particle_spacing);
    particle_config_.particle_radius = config.GetDoubleValue("Particle", "ParticleRadius", 
                                                              particle_config_.particle_radius);
    particle_config_.smoothing_radius = config.GetDoubleValue("Particle", "SmoothingRadius", 
                                                               particle_config_.smoothing_radius);
    particle_config_.cell_size = config.GetDoubleValue("Particle", "CellSize", 
                                                       particle_config_.cell_size);

    return Validate();
}

bool MPSConfig2D::Validate() const {
    try {
        // 验证时间步长
        CHECK_POSITIVE(simulation_config_.time_step, "时间步长 (TimeStep)");

        // 验证总时间
        CHECK_POSITIVE(simulation_config_.total_time, "总时间 (TotalTime)");

        // 验证最大迭代次数
        CHECK_POSITIVE(simulation_config_.max_iterations, "最大迭代次数 (MaxIterations)");

        // 验证密度
        CHECK_POSITIVE(simulation_config_.density, "密度 (Density)");

        // 验证动力学粘性系数
        CHECK_NON_NEGATIVE(simulation_config_.kinematic_viscosity, "动力学粘性系数 (KinematicViscosity)");

        // 验证粒子间距
        CHECK_POSITIVE(particle_config_.particle_spacing, "粒子间距 (ParticleSpacing)");

        // 验证粒子半径
        CHECK_POSITIVE(particle_config_.particle_radius, "粒子半径 (ParticleRadius)");

        // 验证平滑半径
        CHECK_POSITIVE(particle_config_.smoothing_radius, "平滑半径 (SmoothingRadius)");

        // 验证网格单元大小
        CHECK_POSITIVE(particle_config_.cell_size, "网格单元大小 (CellSize)");
        
        // 验证时间步管理参数
        CHECK_POSITIVE(simulation_config_.min_time_step, "最小时间步 (MinTimeStep)");
        CHECK_POSITIVE(simulation_config_.max_time_step, "最大时间步 (MaxTimeStep)");
        CHECK_POSITIVE(simulation_config_.max_cfl, "最大CFL数 (MaxCFL)");
        CHECK_POSITIVE(simulation_config_.output_interval, "输出文件间隔 (OutputInterval)");
        
        // 验证时间步范围合理性
        if (simulation_config_.min_time_step > simulation_config_.max_time_step) {
            throw MPSException("最小时间步必须小于最大时间步");
        }
        if (simulation_config_.time_step < simulation_config_.min_time_step || 
            simulation_config_.time_step > simulation_config_.max_time_step) {
            throw MPSException("初始时间步必须在最小和最大时间步范围内");
        }

        return true;
    } catch (const MPSException& e) {
        // 参数验证失败，输出错误信息
        std::cerr << "配置验证错误: " << e.what() << std::endl;
        return false;
    }
}

} // namespace mps2D
