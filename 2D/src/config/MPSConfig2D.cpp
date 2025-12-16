#include "MPSConfig2D.hpp"

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

    // 验证动力学粘性系数
    if (simulation_config_.kinematic_viscosity < 0.0) {
        return false;
    }

    // 验证粒子间距
    if (particle_config_.particle_spacing <= 0.0) {
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

    // 验证网格单元大小
    if (particle_config_.cell_size <= 0.0) {
        return false;
    }

    return true;
}

} // namespace mps2D
