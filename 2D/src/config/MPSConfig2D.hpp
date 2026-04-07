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
        double time_step = 0.001;              // 初始时间步 (s)
        double total_time = 1.0;               // 总仿真时间 (s)
        int max_iterations = 1000;             // 最大迭代次数
        double density = 1000.0;              // 流体密度 (kg/m³)
        double kinematic_viscosity = 0.001;   // 动力学粘性系数 (m²/s)
        double gravity_x = 0.0;               // 重力加速度 X方向 (m/s²)
        double gravity_y = -9.8;              // 重力加速度 Y方向 (m/s²)
        double ppe_penalty_mu = 0.0;          // PPE罚函数系数μ（<=0时使用程序默认尺度）
        
        // 时间步管理参数
        double min_time_step = 1e-6;          // 最小时间步 (s)
        double max_time_step = 0.01;          // 最大时间步 (s)
        double max_cfl = 0.5;                 // 最大CFL数
        double output_interval = 0.1;         // 输出文件间隔 (s)
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
