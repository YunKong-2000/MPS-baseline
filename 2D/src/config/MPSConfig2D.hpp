#pragma once
#include "ConfigParams.hpp"
#include "../core/Particle.hpp"
#include <string>

namespace mps2D {

// MPS 2D 仿真配置参数类
class MPSConfig2D : public ConfigParams {
public:
    // 文件相关配置
    struct FileConfig {
        std::string input_dir = "data";
        std::string output_dir = "output";
        std::string fluid_particle_file = "fluid_particles_2d.txt";
        std::string solid_particle_file = "solid_particles_2d.txt";
        std::string result_file = "result_2d.txt";
    };

    // 仿真参数配置
    struct SimulationConfig {
        double time_step = 0.001;
        double total_time = 1.0;
        int max_iterations = 1000;
        double density = 1000.0;
        double kinematic_viscosity = 0.001;  // 动力学粘性系数 (m²/s)
        double gravity_x = 0.0;
        double gravity_y = -9.8;
    };

    // 粒子参数配置
    struct ParticleConfig {
        double particle_spacing = 0.02;      // 粒子间距 (m)
        double particle_radius = 0.01;       // 粒子半径 (m)
        double smoothing_radius = 0.042;    // 平滑半径 (m)，通常为2.1倍粒子间距
        double cell_size = 0.084;            // 网格单元大小 (m)，通常为2倍平滑半径
    };

    MPSConfig2D() = default;
    ~MPSConfig2D() override = default;

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

} // namespace mps2D
