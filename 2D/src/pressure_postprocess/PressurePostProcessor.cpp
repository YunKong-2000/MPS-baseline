#include "PressurePostProcessor.hpp"

#include <algorithm>

namespace mps2D {

void PressurePostProcessor::ProcessPressureField(
    FluidParticle& fluid_particles) const {
  const int num_particles = fluid_particles.particle_num;
  if (num_particles <= 0) {
    return;
  }

  const int pressure_size = static_cast<int>(fluid_particles.pressure.size());
  const int surface_type_size =
      static_cast<int>(fluid_particles.surface_type.size());
  const int limit = std::min(num_particles, std::min(pressure_size, surface_type_size));

  for (int i = 0; i < limit; ++i) {
    const SurfaceType surface_type = fluid_particles.surface_type[i];
    if (surface_type == SurfaceType::SPLASH) {
      fluid_particles.pressure[i] = 0.0;
      continue;
    }
    // if (surface_type == SurfaceType::NEAR_SURFACE &&
    //     fluid_particles.pressure[i] < 0.0) {
    //   fluid_particles.pressure[i] = 0.0;
    // }
  }
}

}  // namespace mps2D

