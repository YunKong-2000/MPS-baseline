#include "core/TimeStepManager.h"
#include "core/ErrorHandling.h"

namespace mps2D {

TimeStepManager::TimeStepManager(const MPSConfig2D& config, double particle_spacing)
    : particle_spacing_(particle_spacing) {
    const auto& sim_config = config.GetSimulationConfig();
    
    initial_time_step_ = sim_config.time_step;
    min_time_step_ = sim_config.min_time_step;
    max_time_step_ = sim_config.max_time_step;
    max_cfl_ = sim_config.max_cfl;
    output_interval_ = sim_config.output_interval;
    total_time_ = sim_config.total_time;
    
    // 初始化当前状态
    current_time_step_ = initial_time_step_;
    current_time_ = 0.0;
    next_output_time_ = output_interval_;
    
    // 参数验证
    CHECK_POSITIVE(particle_spacing_, "粒子间距");
    CHECK_POSITIVE(initial_time_step_, "初始时间步");
    CHECK_POSITIVE(min_time_step_, "最小时间步");
    CHECK_POSITIVE(max_time_step_, "最大时间步");
    CHECK_POSITIVE(max_cfl_, "最大CFL数");
    CHECK_POSITIVE(output_interval_, "输出文件间隔");
    CHECK_POSITIVE(total_time_, "总模拟时间");
    
    if (min_time_step_ >= max_time_step_) {
        throw MPSException("最小时间步必须小于最大时间步");
    }
    if (initial_time_step_ < min_time_step_ || initial_time_step_ > max_time_step_) {
        throw MPSException("初始时间步必须在最小和最大时间步范围内");
    }
}

double TimeStepManager::ComputeMaxVelocity(const FluidParticle& fluid_particles) const {
    double max_velocity = 0.0;
    int num_particles = fluid_particles.particle_num;
    
    for (int i = 0; i < num_particles; ++i) {
        if (i < static_cast<int>(fluid_particles.velocity.size())) {
            double velocity_magnitude = ComputeVectorMagnitude(fluid_particles.velocity[i]);
            if (velocity_magnitude > max_velocity) {
                max_velocity = velocity_magnitude;
            }
        }
    }
    
    return max_velocity;
}

double TimeStepManager::ComputeTimeStepFromCFL(double max_velocity) const {
    if (max_velocity <= 1e-10) {
        // 如果速度很小，使用最大时间步
        return max_time_step_;
    }
    
    // CFL = (max_velocity * dt) / characteristic_length
    // dt = (CFL * characteristic_length) / max_velocity
    // 特征长度使用粒子间距
    double suggested_dt = (max_cfl_ * particle_spacing_) / max_velocity;
    
    // 限制在最小和最大时间步范围内
    if (suggested_dt < min_time_step_) {
        return min_time_step_;
    } else if (suggested_dt > max_time_step_) {
        return max_time_step_;
    }
    
    return suggested_dt;
}

double TimeStepManager::GetCFL(double max_velocity) const {
    if (max_velocity <= 1e-10) {
        return 0.0;
    }
    return (max_velocity * current_time_step_) / particle_spacing_;
}

double TimeStepManager::AdjustTimeStep(const FluidParticle& fluid_particles) {
    // 计算当前最大速度
    double max_velocity = ComputeMaxVelocity(fluid_particles);
    
    // 根据CFL条件计算建议的时间步
    double suggested_dt = ComputeTimeStepFromCFL(max_velocity);
    
    // 更新当前时间步
    // 为了稳定性，可以限制时间步的变化率（例如不超过20%的变化）
    double max_change_factor = 1.2;  // 最多增加20%
    double min_change_factor = 0.8;  // 最多减少20%
    
    if (suggested_dt > current_time_step_ * max_change_factor) {
        current_time_step_ = current_time_step_ * max_change_factor;
    } else if (suggested_dt < current_time_step_ * min_change_factor) {
        current_time_step_ = current_time_step_ * min_change_factor;
    } else {
        current_time_step_ = suggested_dt;
    }
    
    // 确保在允许范围内
    if (current_time_step_ < min_time_step_) {
        current_time_step_ = min_time_step_;
    } else if (current_time_step_ > max_time_step_) {
        current_time_step_ = max_time_step_;
    }
    
    return current_time_step_;
}

bool TimeStepManager::ShouldOutput() const {
    // 检查是否到达输出时间
    return current_time_ >= next_output_time_ - 1e-10;  // 使用小的容差避免浮点误差
}

void TimeStepManager::AdvanceTime() {
    current_time_ += current_time_step_;
    
    // 确保不超过总时间
    if (current_time_ > total_time_) {
        current_time_ = total_time_;
    }
}

void TimeStepManager::MarkOutputDone() {
    // 更新下次输出时间
    next_output_time_ += output_interval_;
}

void TimeStepManager::Reset() {
    current_time_step_ = initial_time_step_;
    current_time_ = 0.0;
    next_output_time_ = output_interval_;
}

} // namespace mps2D
