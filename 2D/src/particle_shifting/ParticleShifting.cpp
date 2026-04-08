#include "ParticleShifting.hpp"

#include "core/MPSUtils.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace mps2D {

namespace {

constexpr double kSmallEps = 1e-10;
constexpr double kShiftLambda = 0.1;
constexpr double kDim = 2.0;

double ComputePSKernelWeight(double distance, double variable_cutoff_radius) {
  if (distance <= kSmallEps || distance >= variable_cutoff_radius) {
    return 0.0;
  }
  return variable_cutoff_radius / distance - 1.0;
}

double ComputeReferenceNumberDensity(
    double smoothing_radius,
    double particle_spacing) {
  if (smoothing_radius <= kSmallEps || particle_spacing <= kSmallEps) {
    return 1.0;
  }

  // 基于均匀正方形虚拟粒子排布计算中心粒子粒子数密度:
  // 网格边长 = (2 * floor(r_e / l0) + 1) * l0
  const int span = static_cast<int>(std::floor(smoothing_radius / particle_spacing));
  double n0 = 0.0;
  for (int iy = -span; iy <= span; ++iy) {
    for (int ix = -span; ix <= span; ++ix) {
      if (ix == 0 && iy == 0) {
        continue;
      }
      const double dx = static_cast<double>(ix) * particle_spacing;
      const double dy = static_cast<double>(iy) * particle_spacing;
      const double dist = std::sqrt(dx * dx + dy * dy);
      n0 += WeightFunction(dist, smoothing_radius);
    }
  }

  return (n0 > kSmallEps) ? n0 : 1.0;
}

double ComputeVariableCutoffRadius(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double particle_spacing,
    double smoothing_radius) {
  const double2& pos_i = fluid_particles.position[particle_idx];
  double d_min = std::numeric_limits<double>::max();
  bool has_wall_neighbor = false;

  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    if (j < 0 || j >= static_cast<int>(solid_particles.position.size())) {
      continue;
    }
    has_wall_neighbor = true;
    const double dist = ComputeDistance(pos_i, solid_particles.position[j]);
    if (dist < d_min) {
      d_min = dist;
    }
  }

  if (!has_wall_neighbor) {
    return smoothing_radius;
  }

  return std::max(1.2 * particle_spacing, d_min);
}

double2 ComputeFreeSurfaceNormal(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double variable_cutoff_radius) {
  const double2& pos_i = fluid_particles.position[particle_idx];
  double2 sum_vec{0.0, 0.0};

  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    if (j < 0 || j >= fluid_particles.particle_num || j == particle_idx) {
      continue;
    }
    const double2& pos_j = fluid_particles.position[j];
    const double dx = pos_j.x - pos_i.x;
    const double dy = pos_j.y - pos_i.y;
    const double dist = std::sqrt(dx * dx + dy * dy);
    if (dist <= kSmallEps || dist >= variable_cutoff_radius) {
      continue;
    }
    const double weight = WeightFunction(dist, variable_cutoff_radius);
    sum_vec.x += (dx / dist) * weight;
    sum_vec.y += (dy / dist) * weight;
  }

  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    if (j < 0 || j >= static_cast<int>(solid_particles.position.size())) {
      continue;
    }
    const double2& pos_j = solid_particles.position[j];
    const double dx = pos_j.x - pos_i.x;
    const double dy = pos_j.y - pos_i.y;
    const double dist = std::sqrt(dx * dx + dy * dy);
    if (dist <= kSmallEps || dist >= variable_cutoff_radius) {
      continue;
    }
    const double weight = WeightFunction(dist, variable_cutoff_radius);
    sum_vec.x += (dx / dist) * weight;
    sum_vec.y += (dy / dist) * weight;
  }

  const double norm = std::sqrt(sum_vec.x * sum_vec.x + sum_vec.y * sum_vec.y);
  if (norm <= kSmallEps) {
    return {0.0, 0.0};
  }
  return {-sum_vec.x / norm, -sum_vec.y / norm};
}

double2 ComputeOPSDisplacement(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double particle_spacing,
    double reference_density,
    double variable_cutoff_radius) {
  const double2& pos_i = fluid_particles.position[particle_idx];
  double2 summation{0.0, 0.0};

  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    if (j < 0 || j >= fluid_particles.particle_num || j == particle_idx) {
      continue;
    }
    const double2& pos_j = fluid_particles.position[j];
    const double dx = pos_j.x - pos_i.x;
    const double dy = pos_j.y - pos_i.y;
    const double dist2 = dx * dx + dy * dy;
    const double dist = std::sqrt(dist2);
    if (dist2 <= kSmallEps || dist >= variable_cutoff_radius) {
      continue;
    }
    const double weight_ps = ComputePSKernelWeight(dist, variable_cutoff_radius);
    summation.x += weight_ps * dx / dist2;
    summation.y += weight_ps * dy / dist2;
  }

  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    if (j < 0 || j >= static_cast<int>(solid_particles.position.size())) {
      continue;
    }
    const double2& pos_j = solid_particles.position[j];
    const double dx = pos_j.x - pos_i.x;
    const double dy = pos_j.y - pos_i.y;
    const double dist2 = dx * dx + dy * dy;
    const double dist = std::sqrt(dist2);
    if (dist2 <= kSmallEps || dist >= variable_cutoff_radius) {
      continue;
    }
    const double weight_ps = ComputePSKernelWeight(dist, variable_cutoff_radius);
    summation.x += weight_ps * dx / dist2;
    summation.y += weight_ps * dy / dist2;
  }

  const double coeff =
      -kShiftLambda * kDim * particle_spacing * particle_spacing / reference_density;
  return {coeff * summation.x, coeff * summation.y};
}

double2 ComputeSlipDisplacement(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const std::vector<double2>& velocity_after_pressure,
    double variable_cutoff_radius,
    double rc) {
  const double2& pos_i = fluid_particles.position[particle_idx];
  const double2& vel_i = velocity_after_pressure[particle_idx];
  double2 slip{0.0, 0.0};

  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    if (j < 0 || j >= fluid_particles.particle_num || j == particle_idx ||
        j >= static_cast<int>(velocity_after_pressure.size())) {
      continue;
    }
    const double2& pos_j = fluid_particles.position[j];
    const double2& vel_j = velocity_after_pressure[j];
    const double dx = pos_j.x - pos_i.x;
    const double dy = pos_j.y - pos_i.y;
    const double dist2 = dx * dx + dy * dy;
    const double dist = std::sqrt(dist2);
    if (dist <= kSmallEps || dist >= variable_cutoff_radius || dist >= rc) {
      continue;
    }

    const double inv_dist = 1.0 / dist;
    const double ex = dx * inv_dist;
    const double ey = dy * inv_dist;
    const double2 du{vel_i.x - vel_j.x, vel_i.y - vel_j.y};
    const double radial_component = du.x * ex + du.y * ey;
    const double2 tangent{du.x - radial_component * ex, du.y - radial_component * ey};
    const double tangent_norm = std::sqrt(tangent.x * tangent.x + tangent.y * tangent.y);
    if (tangent_norm <= kSmallEps) {
      continue;
    }

    const double factor = 0.5 * std::sqrt(std::max(0.0, rc * rc - dist2));
    slip.x += tangent.x / tangent_norm * factor;
    slip.y += tangent.y / tangent_norm * factor;
  }

  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    if (j < 0 || j >= static_cast<int>(solid_particles.position.size())) {
      continue;
    }
    const double2& pos_j = solid_particles.position[j];
    const double2& vel_j = solid_particles.velocity[j];
    const double dx = pos_j.x - pos_i.x;
    const double dy = pos_j.y - pos_i.y;
    const double dist2 = dx * dx + dy * dy;
    const double dist = std::sqrt(dist2);
    if (dist <= kSmallEps || dist >= variable_cutoff_radius || dist >= rc) {
      continue;
    }

    const double inv_dist = 1.0 / dist;
    const double ex = dx * inv_dist;
    const double ey = dy * inv_dist;
    const double2 du{vel_i.x - vel_j.x, vel_i.y - vel_j.y};
    const double radial_component = du.x * ex + du.y * ey;
    const double2 tangent{du.x - radial_component * ex, du.y - radial_component * ey};
    const double tangent_norm = std::sqrt(tangent.x * tangent.x + tangent.y * tangent.y);
    if (tangent_norm <= kSmallEps) {
      continue;
    }

    const double factor = 0.5 * std::sqrt(std::max(0.0, rc * rc - dist2));
    slip.x += tangent.x / tangent_norm * factor;
    slip.y += tangent.y / tangent_norm * factor;
  }

  return slip;
}

}  // namespace

std::vector<double2> ParticleShifting::ComputeShiftingDisplacement(
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const std::vector<double2>& velocity_after_pressure,
    double smoothing_radius,
    double particle_spacing) {
  const int num_particles = fluid_particles.particle_num;
  std::vector<double2> displacement(num_particles, {0.0, 0.0});
  if (num_particles <= 0 || smoothing_radius <= kSmallEps ||
      particle_spacing <= kSmallEps) {
    return displacement;
  }

  const double n0 = ComputeReferenceNumberDensity(
      smoothing_radius, particle_spacing);
  const double rc = 0.7 * particle_spacing;
  const double max_displacement = 0.1 * particle_spacing;

  std::vector<double2> normals(num_particles, {0.0, 0.0});
  std::vector<double> variable_cutoff_radius(num_particles, smoothing_radius);
  for (int i = 0; i < num_particles; ++i) {
    variable_cutoff_radius[i] = ComputeVariableCutoffRadius(
        i, fluid_particles, solid_particles, particle_spacing, smoothing_radius);
    normals[i] = ComputeFreeSurfaceNormal(
        i, fluid_particles, solid_particles, variable_cutoff_radius[i]);
  }

  for (int i = 0; i < num_particles; ++i) {
    const SurfaceType surface_type =
        (i < static_cast<int>(fluid_particles.surface_type.size()))
            ? fluid_particles.surface_type[i]
            : SurfaceType::INNER;
    if (surface_type == SurfaceType::SPLASH) {
      displacement[i] = {0.0, 0.0};
      continue;
    }

    const double r_ve = variable_cutoff_radius[i];
    const double2 ops = ComputeOPSDisplacement(
        i, fluid_particles, solid_particles, particle_spacing, n0, r_ve);
    const double2 slip = ComputeSlipDisplacement(
        i, fluid_particles, solid_particles, velocity_after_pressure, r_ve, rc);
    double2 total{ops.x + slip.x, ops.y + slip.y};

    if (surface_type == SurfaceType::SURFACE) {
      const double2& n = normals[i];
      const double ndot = total.x * n.x + total.y * n.y;
      const double2 projected{
          total.x - ndot * n.x,
          total.y - ndot * n.y};
      const double factor = (kDim - 1.0) / kDim;
      total.x = factor * projected.x;
      total.y = factor * projected.y;
    }

    const double mag = std::sqrt(total.x * total.x + total.y * total.y);
    if (mag > max_displacement && mag > kSmallEps) {
      const double scale = max_displacement / mag;
      total.x *= scale;
      total.y *= scale;
    }
    displacement[i] = total;
  }

  return displacement;
}

}  // namespace mps2D

