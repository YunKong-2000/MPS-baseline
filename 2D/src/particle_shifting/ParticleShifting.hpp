#pragma once

#include "../core/Particle.hpp"
#include "core/Types.h"
#include <vector>

namespace mps2D {

class ParticleShifting {
 public:
  // 基于可变截断半径PS算法计算所有流体粒子的shifting位移
  static std::vector<double2> ComputeShiftingDisplacement(
      const FluidParticle& fluid_particles,
      const SolidParticle& solid_particles,
      const std::vector<double2>& velocity_after_pressure,
      double smoothing_radius,
      double particle_spacing);
};

}  // namespace mps2D

