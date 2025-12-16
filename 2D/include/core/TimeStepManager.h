#pragma once
#include "../../src/config/MPSConfig2D.hpp"
#include "../../src/core/Particle.hpp"
#include "core/MPSUtils.h"
#include <cmath>

namespace mps2D {

/**
 * @brief 时间步管理类
 * 
 * 负责管理模拟过程中的时间步动态调整、模拟时间记录和输出文件间隔控制。
 * 根据CFL条件和粒子最大速度自动调整时间步大小。
 */
class TimeStepManager {
public:
    /**
     * @brief 构造函数
     * @param config MPS配置对象，包含时间步管理相关参数
     * @param particle_spacing 粒子间距，用于CFL计算
     */
    TimeStepManager(const MPSConfig2D& config, double particle_spacing);
    
    /**
     * @brief 析构函数
     */
    ~TimeStepManager() = default;
    
    /**
     * @brief 获取当前时间步
     * @return 当前时间步大小 (s)
     */
    double GetTimeStep() const { return current_time_step_; }
    
    /**
     * @brief 获取当前模拟时间
     * @return 当前模拟时间 (s)
     */
    double GetCurrentTime() const { return current_time_; }
    
    /**
     * @brief 获取总模拟时间
     * @return 总模拟时间 (s)
     */
    double GetTotalTime() const { return total_time_; }
    
    /**
     * @brief 检查模拟是否完成
     * @return true表示模拟已完成，false表示继续
     */
    bool IsSimulationFinished() const { return current_time_ >= total_time_; }
    
    /**
     * @brief 检查是否需要输出文件
     * @return true表示需要输出，false表示不需要
     */
    bool ShouldOutput() const;
    
    /**
     * @brief 根据粒子最大速度动态调整时间步
     * @param fluid_particles 流体粒子对象，用于计算最大速度
     * @return 调整后的时间步大小 (s)
     */
    double AdjustTimeStep(const FluidParticle& fluid_particles);
    
    /**
     * @brief 更新模拟时间（在完成一个时间步后调用）
     */
    void AdvanceTime();
    
    /**
     * @brief 标记已输出文件（在输出文件后调用）
     */
    void MarkOutputDone();
    
    /**
     * @brief 重置时间步管理器（用于重新开始模拟）
     */
    void Reset();
    
    /**
     * @brief 获取当前CFL数
     * @param max_velocity 粒子最大速度
     * @return CFL数
     */
    double GetCFL(double max_velocity) const;

private:
    // 配置参数
    double initial_time_step_;      // 初始时间步
    double min_time_step_;          // 最小时间步
    double max_time_step_;          // 最大时间步
    double max_cfl_;                // 最大CFL数
    double output_interval_;        // 输出文件间隔
    double total_time_;             // 总模拟时间
    double particle_spacing_;       // 粒子间距（用于CFL计算）
    
    // 当前状态
    double current_time_step_;      // 当前时间步
    double current_time_;           // 当前模拟时间
    double next_output_time_;       // 下次输出时间
    
    /**
     * @brief 计算粒子最大速度
     * @param fluid_particles 流体粒子对象
     * @return 最大速度大小 (m/s)
     */
    double ComputeMaxVelocity(const FluidParticle& fluid_particles) const;
    
    /**
     * @brief 根据CFL条件计算时间步
     * @param max_velocity 最大速度
     * @return 建议的时间步大小
     */
    double ComputeTimeStepFromCFL(double max_velocity) const;
};

} // namespace mps2D
