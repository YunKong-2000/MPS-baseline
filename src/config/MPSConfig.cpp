#include "MPSConfig.hpp"

namespace mps {

bool MPSConfig::LoadFromConfig(const SimpleIni& config) {
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
    simulation_config_.viscosity = config.GetDoubleValue("Simulation", "Viscosity", 
                                                         simulation_config_.viscosity);
    simulation_config_.gravity_x = config.GetDoubleValue("Simulation", "GravityX", 
                                                         simulation_config_.gravity_x);
    simulation_config_.gravity_y = config.GetDoubleValue("Simulation", "GravityY", 
                                                         simulation_config_.gravity_y);
    simulation_config_.gravity_z = config.GetDoubleValue("Simulation", "GravityZ", 
                                                         simulation_config_.gravity_z);

    // 加载粒子配置
    particle_config_.particle_radius = config.GetDoubleValue("Particle", "ParticleRadius", 
                                                             particle_config_.particle_radius);
    particle_config_.smoothing_radius = config.GetDoubleValue("Particle", "SmoothingRadius", 
                                                               particle_config_.smoothing_radius);
    particle_config_.particle_count = config.GetIntValue("Particle", "ParticleCount", 
                                                        particle_config_.particle_count);

    return Validate();
}

bool MPSConfig::Validate() const {
    // 验证时间步长
    if (simulation_config_.time_step <= 0.0) {
        return false;
    }

    // 验证总时间
    if (simulation_config_.total_time <= 0.0) {
        return false;
    }

    // 验证最大迭代次数
    if (simulation_config_.max_iterations <= 0) {
        return false;
    }

    // 验证密度
    if (simulation_config_.density <= 0.0) {
        return false;
    }

    // 验证粒子半径
    if (particle_config_.particle_radius <= 0.0) {
        return false;
    }

    // 验证平滑半径
    if (particle_config_.smoothing_radius <= 0.0) {
        return false;
    }

    // 验证粒子数量
    if (particle_config_.particle_count <= 0) {
        return false;
    }

    return true;
}

} // namespace mps

