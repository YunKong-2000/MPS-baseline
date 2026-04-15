#pragma once

#include "../core/Particle.hpp"

namespace mps2D {

// 压力结果后处理模块
// 用于在PPE求解后统一执行压力裁剪规则。
class PressurePostProcessor {
 public:
  PressurePostProcessor() = default;
  ~PressurePostProcessor() = default;

  // 对压力场执行后处理：
  // 1) splash 粒子压力强制为 0；
  // 2) 近自由面粒子若出现负压，裁剪为 0。
  void ProcessPressureField(FluidParticle& fluid_particles) const;
};

}  // namespace mps2D

