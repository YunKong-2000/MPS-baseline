#pragma once
#include "ConfigParams.hpp"
#include "core/Types.h"
#include <string>

namespace mps {

// MPS 仿真配置参数类
class MPSConfig : public ConfigParams {
public:
    // 文件相关配置
    struct FileConfig {
        std::string input_dir = "data";
        std::string output_dir = "output";
        std::string fluid_particle_file = "fluid_particles.txt";
        std::string solid_particle_file = "solid_particles.txt";
        std::string result_file = "result.txt";
    };

    // 仿真参数配置
    struct SimulationConfig {
        double time_step = 0.001;
        double total_time = 1.0;
        int max_iterations = 1000;
        double density = 1000.0;
        double viscosity = 0.001;
        double gravity_x = 0.0;
        double gravity_y = 0.0;
        double gravity_z = -9.8;
    };

    // 粒子参数配置
    struct ParticleConfig {
        double particle_radius = 0.01;
        double smoothing_radius = 3;
        double cell_size = 4;
        int particle_count = 1000;
    };

    MPSConfig() = default;
    ~MPSConfig() override = default;

    // 从配置管理器加载参数
    bool LoadFromConfig(const SimpleIni& config) override;

    // 验证参数有效性
    bool Validate() const override;

    // 获取配置
    const FileConfig& GetFileConfig() const { return file_config_; }
    const SimulationConfig& GetSimulationConfig() const { return simulation_config_; }
    const ParticleConfig& GetParticleConfig() const { return particle_config_; }

    // 设置配置（用于测试或手动设置）
    void SetFileConfig(const FileConfig& config) { file_config_ = config; }
    void SetSimulationConfig(const SimulationConfig& config) { simulation_config_ = config; }
    void SetParticleConfig(const ParticleConfig& config) { particle_config_ = config; }

private:
    FileConfig file_config_;
    SimulationConfig simulation_config_;
    ParticleConfig particle_config_;
};

} // namespace mps

