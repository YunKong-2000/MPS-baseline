#include "ParticleShifting.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace mps2D {

namespace {
double WeightFunctionPS(double distance, double re_ps) {
  if (distance >= re_ps || re_ps <= 0.0) {
    return 0.0;
  }
  return re_ps / distance - 1.0;
}
}  // namespace

double ParticleShifting::ComputeReferenceNumberDensity(
    double particle_spacing,
    double re_ps) const {
  if (particle_spacing <= 0.0 || re_ps <= 0.0) {
    return 1.0;
  }
  // 9x9 规则均匀布点，取中心粒子的数密度作为 n0
  const int kGrid = 9;
  const int center = kGrid / 2;
  double n0 = 0.0;
  for (int iy = 0; iy < kGrid; ++iy) {
    for (int ix = 0; ix < kGrid; ++ix) {
      if (ix == center && iy == center) {
        continue;
      }
      const double dx = static_cast<double>(ix - center) * particle_spacing;
      const double dy = static_cast<double>(iy - center) * particle_spacing;
      const double dist = std::sqrt(dx * dx + dy * dy);
      if (dist < kEps) {
        continue;
      }
      n0 += WeightFunctionPS(dist, re_ps);
    }
  }
  return std::max(n0, 1e-8);
}

void ParticleShifting::ComputeMinWallDistanceAndRve(
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double smoothing_radius,
    double particle_spacing,
    std::vector<double>& rve) const {
  const int n = fluid_particles.particle_num;
  rve.assign(n, smoothing_radius);

  for (int i = 0; i < n; ++i) {
    const double2& pos_i = fluid_particles.position[i];
    double dmin_wall = std::numeric_limits<double>::max();

    for (int j : fluid_particles.solid_neighbour_list[i]) {
      const double dist = ComputeDistance(pos_i, solid_particles.position[j]);
      dmin_wall = std::min(dmin_wall, dist);
    }

    if (dmin_wall == std::numeric_limits<double>::max()) {
      rve[i] = smoothing_radius;
    } else {
      rve[i] = std::max(1.2 * particle_spacing, dmin_wall);
    }
  }
}

void ParticleShifting::ComputeFreeSurfaceNormals(
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const std::vector<double>& rve,
    std::vector<double2>& normals) const {
  const int n = fluid_particles.particle_num;
  normals.assign(n, double2{0.0, 0.0});

  for (int i = 0; i < n; ++i) {
    if (fluid_particles.surface_type[i] != SurfaceType::SURFACE) {
      continue;
    }
    const double2& pos_i = fluid_particles.position[i];
    double2 accum{0.0, 0.0};
    for (int j : fluid_particles.fluid_neighbour_list[i]) {
      const double2& pos_j = fluid_particles.position[j];
      const double2 rij{pos_j.x - pos_i.x, pos_j.y - pos_i.y};
      const double dist = ComputeDistance(pos_i, pos_j);
      if (dist < kEps || dist > rve[i]) {
        continue;
      }
      const double w = WeightFunction(dist, rve[i]);
      accum.x += (rij.x / dist) * w;
      accum.y += (rij.y / dist) * w;
    }
    for (int j : fluid_particles.solid_neighbour_list[i]) {
      const double2& pos_j = solid_particles.position[j];
      const double2 rij{pos_j.x - pos_i.x, pos_j.y - pos_i.y};
      const double dist = ComputeDistance(pos_i, pos_j);
      if (dist < kEps || dist > rve[i]) {
        continue;
      }
      const double w = WeightFunction(dist, rve[i]);
      accum.x += (rij.x / dist) * w;
      accum.y += (rij.y / dist) * w;
    }
    const double mag = ComputeVectorMagnitude(accum);
    if (mag > kEps) {
      normals[i] = double2{-accum.x / mag, -accum.y / mag};
    }
  }
}

double2 ParticleShifting::ComputeBasicShift(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double rve,
    double l0,
    double n0) const {
  const double2& pos_i = fluid_particles.position[particle_idx];
  double2 sum{0.0, 0.0};
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];
    const double2 rij{pos_j.x - pos_i.x, pos_j.y - pos_i.y};
    const double dist2 = rij.x * rij.x + rij.y * rij.y;
    if (dist2 < kEps) {
      continue;
    }
    const double dist = std::sqrt(dist2);
    const double w_ps = WeightFunctionPS(dist, rve);
    if (w_ps <= 0.0) {
      continue;
    }
    sum.x += w_ps * rij.x / dist2;
    sum.y += w_ps * rij.y / dist2;
  }
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    const double2& pos_j = solid_particles.position[j];
    const double2 rij{pos_j.x - pos_i.x, pos_j.y - pos_i.y};
    const double dist2 = rij.x * rij.x + rij.y * rij.y;
    if (dist2 < kEps) {
      continue;
    }
    const double dist = std::sqrt(dist2);
    const double w_ps = WeightFunctionPS(dist, rve);
    if (w_ps <= 0.0) {
      continue;
    }
    sum.x += w_ps * rij.x / dist2;
    sum.y += w_ps * rij.y / dist2;
  }

  const double coeff =
      -kLambdaShift * l0 * l0 * (static_cast<double>(kDimension) / n0);
  return {coeff * sum.x, coeff * sum.y};
}

double2 ParticleShifting::ComputeSlipShift(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double rve,
    double l0) const {
  const double2& pos_i = fluid_particles.position[particle_idx];
  const double2& u_i = fluid_particles.velocity[particle_idx];

  // 文档：r_c = 0.7 * l_0
  const double rc = 0.7 * l0;
  const double rc2 = rc * rc;

  const double rve2 = rve * rve;
  double2 sum{0.0, 0.0};

  auto AccumulateFromNeighbor = [&](const double2& pos_j,
                                     const double2& u_j) {
    const double2 rij{pos_j.x - pos_i.x, pos_j.y - pos_i.y};
    const double dist2 = rij.x * rij.x + rij.y * rij.y;
    if (dist2 < kEps || dist2 >= rve2 || dist2 >= rc2) {
      return;
    }

    const double dist = std::sqrt(dist2);
    const double2 e_ij{rij.x / dist, rij.y / dist};

    // 切向投影：T_ij(u_i - u_j) = (u_i - u_j) - e_ij * ((u_i - u_j) · e_ij)
    const double2 rel_u{u_i.x - u_j.x, u_i.y - u_j.y};
    const double dot_re = rel_u.x * e_ij.x + rel_u.y * e_ij.y;
    const double2 tangential_vec{rel_u.x - dot_re * e_ij.x,
                                  rel_u.y - dot_re * e_ij.y};

    const double tangential_mag = ComputeVectorMagnitude(tangential_vec);
    if (tangential_mag < kEps) {
      return;  // 保护：分母为 0 时 tau_ij = 0
    }

    const double2 tau_ij{tangential_vec.x / tangential_mag,
                          tangential_vec.y / tangential_mag};

    const double sqrt_term = std::sqrt(std::max(0.0, rc2 - dist2));
    sum.x += 0.5 * tau_ij.x * sqrt_term;
    sum.y += 0.5 * tau_ij.y * sqrt_term;
  };

  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    AccumulateFromNeighbor(fluid_particles.position[j], fluid_particles.velocity[j]);
  }
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    AccumulateFromNeighbor(solid_particles.position[j], solid_particles.velocity[j]);
  }

  return sum;
}

double2 ParticleShifting::ProjectForSurfaceConstraint(
    const double2& shift,
    const double2& normal) const {
  const double nn = normal.x * normal.x + normal.y * normal.y;
  if (nn < kEps) {
    return {0.0, 0.0};
  }

  const double dot_sn = shift.x * normal.x + shift.y * normal.y;
  const double2 tangential{
      shift.x - dot_sn * normal.x,
      shift.y - dot_sn * normal.y};
  const double factor = static_cast<double>(kDimension - 1) /
                        static_cast<double>(kDimension);
  return {factor * tangential.x, factor * tangential.y};
}

double2 ParticleShifting::ApplyMagnitudeLimiter(
    const double2& shift,
    double particle_spacing) const {
  const double mag = ComputeVectorMagnitude(shift);
  if (mag < kEps) {
    return shift;
  }
  const double max_mag = kMaxShiftFactor * particle_spacing;
  if (mag <= max_mag) {
    return shift;
  }
  const double scale = max_mag / mag;
  return {shift.x * scale, shift.y * scale};
}

double2 ParticleShifting::ApplyWallCollisionLimiter(
    int particle_idx,
    const double2& shift,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double particle_spacing) const {
  if (fluid_particles.solid_neighbour_list[particle_idx].empty()) {
    return shift;
  }

  const double2& pos_i = fluid_particles.position[particle_idx];
  const double2 predicted{pos_i.x + shift.x, pos_i.y + shift.y};
  const double d_safe = kSafeDistanceFactor * particle_spacing;

  bool close_to_wall = false;
  double min_dist = std::numeric_limits<double>::max();
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    const double dist_pred = ComputeDistance(predicted, solid_particles.position[j]);
    min_dist = std::min(min_dist, dist_pred);
  }
  if (min_dist < d_safe) {
    close_to_wall = true;
  }

  double2 limited = shift;
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    const double2 wall_pos = solid_particles.position[j];
    const double2 wall_n = NormalizeVector(solid_particles.normal_vector[j]);
    if (ComputeVectorMagnitude(wall_n) < kEps) {
      continue;
    }

    const double2 rp{predicted.x - wall_pos.x, predicted.y - wall_pos.y};
    const double signed_dist_pred = rp.x * wall_n.x + rp.y * wall_n.y;
    const double delta_n = limited.x * wall_n.x + limited.y * wall_n.y;

    if (close_to_wall || signed_dist_pred < d_safe) {
      if (delta_n < 0.0) {
        limited.x -= delta_n * wall_n.x;
        limited.y -= delta_n * wall_n.y;
      }
    }
  }
  return limited;
}

Eigen::Matrix2d ParticleShifting::ComputeVelocityGradient(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE,
                        CorrectiveMatrix::MATRIX_SIZE>& corrective_matrix,
    double smoothing_radius) const {
  const double2& pos_i = fluid_particles.position[particle_idx];
  const double2& vel_i = fluid_particles.velocity[particle_idx];

  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> M0 =
      corrective_matrix.row(0);
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> M1 =
      corrective_matrix.row(1);

  double dux_dx_sum = 0.0;
  double dux_dy_sum = 0.0;
  double duy_dx_sum = 0.0;
  double duy_dy_sum = 0.0;

  CorrectiveMatrix cm;
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];
    const double dist = ComputeDistance(pos_i, pos_j);
    if (dist < kEps || dist > smoothing_radius) {
      continue;
    }
    const double2& vel_j = fluid_particles.velocity[j];
    const double dvx = vel_j.x - vel_i.x;
    const double dvy = vel_j.y - vel_i.y;
    const double w = WeightFunction(dist, smoothing_radius);

    Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis =
        cm.ComputeBasisFunctions(pos_j.x - pos_i.x, pos_j.y - pos_i.y,
                                 smoothing_radius);
    const double M0P = (M0 * basis)(0, 0);
    const double M1P = (M1 * basis)(0, 0);

    dux_dx_sum += w * dvx * M0P;
    dux_dy_sum += w * dvx * M1P;
    duy_dx_sum += w * dvy * M0P;
    duy_dy_sum += w * dvy * M1P;
  }

  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    const double2& pos_j = solid_particles.position[j];
    const double dist = ComputeDistance(pos_i, pos_j);
    if (dist < kEps || dist > smoothing_radius) {
      continue;
    }
    const double2 vel_j = solid_particles.velocity[j];
    const double dvx = vel_j.x - vel_i.x;
    const double dvy = vel_j.y - vel_i.y;
    const double w = WeightFunction(dist, smoothing_radius);

    Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis =
        cm.ComputeBasisFunctions(pos_j.x - pos_i.x, pos_j.y - pos_i.y,
                                 smoothing_radius);
    const double M0P = (M0 * basis)(0, 0);
    const double M1P = (M1 * basis)(0, 0);

    dux_dx_sum += w * dvx * M0P;
    dux_dy_sum += w * dvx * M1P;
    duy_dx_sum += w * dvy * M0P;
    duy_dy_sum += w * dvy * M1P;
  }

  const double inv_rs = 1.0 / smoothing_radius;
  Eigen::Matrix2d grad_u;
  grad_u(0, 0) = inv_rs * dux_dx_sum;
  grad_u(0, 1) = inv_rs * dux_dy_sum;
  grad_u(1, 0) = inv_rs * duy_dx_sum;
  grad_u(1, 1) = inv_rs * duy_dy_sum;
  return grad_u;
}

void ParticleShifting::Apply(
    FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE,
                                    CorrectiveMatrix::MATRIX_SIZE>>&
        corrective_matrices_velocity,
    double smoothing_radius,
    double particle_spacing,
    std::vector<double2>* out_ps_displacement,
    std::vector<double>* out_rve,
    std::vector<double>* out_temp_number_density,
    std::vector<double2>* out_free_surface_normals,
    std::vector<double>* out_velocity_grad_xx,
    std::vector<double>* out_velocity_grad_xy,
    std::vector<double>* out_velocity_grad_yx,
    std::vector<double>* out_velocity_grad_yy) const {
  const int n = fluid_particles.particle_num;
  if (n <= 0) {
    return;
  }
  if (static_cast<int>(corrective_matrices_velocity.size()) != n) {
    return;
  }

  const double n0 = ComputeReferenceNumberDensity(particle_spacing, smoothing_radius);

  std::vector<double> rve;
  ComputeMinWallDistanceAndRve(
      fluid_particles, solid_particles, smoothing_radius, particle_spacing, rve);

  std::vector<double2> normals;
  ComputeFreeSurfaceNormals(fluid_particles, solid_particles, rve, normals);
  if (out_free_surface_normals != nullptr) {
    *out_free_surface_normals = normals;  // PS 前的法向（非自由面默认 0）
  }

  // PS 前的速度梯度输出（默认初始化为 0，便于 SPLASH 粒子等情况可视化）
  if (out_velocity_grad_xx != nullptr) {
    out_velocity_grad_xx->assign(n, 0.0);
  }
  if (out_velocity_grad_xy != nullptr) {
    out_velocity_grad_xy->assign(n, 0.0);
  }
  if (out_velocity_grad_yx != nullptr) {
    out_velocity_grad_yx->assign(n, 0.0);
  }
  if (out_velocity_grad_yy != nullptr) {
    out_velocity_grad_yy->assign(n, 0.0);
  }

  // 临时粒子数密度：使用每个粒子的可变截断半径 r_ve 定义邻域
  std::vector<double> temp_number_density(n, 0.0);
  for (int i = 0; i < n; ++i) {
    // 飞溅粒子不参与 PS 贡献计算
    if (static_cast<size_t>(i) < fluid_particles.surface_type.size() &&
        fluid_particles.surface_type[i] == SurfaceType::SPLASH) {
      continue;
    }
    const double2& pos_i = fluid_particles.position[i];
    const double rve_i = rve[i];
    double density_i = 0.0;
    for (int j : fluid_particles.fluid_neighbour_list[i]) {
      const double dist = ComputeDistance(pos_i, fluid_particles.position[j]);
      if (dist < kEps || dist >= rve_i) {
        continue;
      }
      density_i += WeightFunctionPS(dist, rve_i);
    }
    temp_number_density[i] = density_i;
  }

  std::vector<double2> displacements(n, double2{0.0, 0.0});
  std::vector<double2> new_position = fluid_particles.position;

  for (int i = 0; i < n; ++i) {
    SurfaceType type = SurfaceType::INNER;
    if (static_cast<size_t>(i) < fluid_particles.surface_type.size()) {
      type = fluid_particles.surface_type[i];
    }
    if (type == SurfaceType::SPLASH) {
      continue;
    }

    const double2 f_basic = ComputeBasicShift(i, fluid_particles, solid_particles,
                                              rve[i],
                                              particle_spacing, n0);

    // 文档：δr^PS = δr^OPS + δr^SL
    const double2 f_slip =
        ComputeSlipShift(i, fluid_particles, solid_particles, rve[i],
                         particle_spacing);
    double2 shift{f_basic.x + f_slip.x, f_basic.y + f_slip.y};

    if (type == SurfaceType::SURFACE) {
      shift = ProjectForSurfaceConstraint(shift, normals[i]);
    }

    shift = ApplyMagnitudeLimiter(shift, particle_spacing);
    shift = ApplyWallCollisionLimiter(
        i, shift, fluid_particles, solid_particles, particle_spacing);
    displacements[i] = shift;

    const Eigen::Matrix2d grad_u = ComputeVelocityGradient(
        i, fluid_particles, solid_particles, corrective_matrices_velocity[i],
        smoothing_radius);

    if (out_velocity_grad_xx != nullptr) {
      (*out_velocity_grad_xx)[i] = grad_u(0, 0);
    }
    if (out_velocity_grad_xy != nullptr) {
      (*out_velocity_grad_xy)[i] = grad_u(0, 1);
    }
    if (out_velocity_grad_yx != nullptr) {
      (*out_velocity_grad_yx)[i] = grad_u(1, 0);
    }
    if (out_velocity_grad_yy != nullptr) {
      (*out_velocity_grad_yy)[i] = grad_u(1, 1);
    }

    // 暂时关闭 PS 后的一阶泰勒速度修正：u ← u + ∇u·δr（位置仍按 δr 更新）
    new_position[i].x = fluid_particles.position[i].x + shift.x;
    new_position[i].y = fluid_particles.position[i].y + shift.y;
  }

  fluid_particles.position.swap(new_position);

  if (out_ps_displacement != nullptr) {
    *out_ps_displacement = displacements;
  }
  if (out_rve != nullptr) {
    *out_rve = rve;
  }
  if (out_temp_number_density != nullptr) {
    // PS 前的临时粒子数密度（仅速度/位置尚未更新时的值）
    *out_temp_number_density = temp_number_density;
  }
}

}  // namespace mps2D

