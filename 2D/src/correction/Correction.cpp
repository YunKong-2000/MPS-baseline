#include "Correction.hpp"
#include "../lsmps/CorrectiveMatrix.hpp"
#include <cmath>
#include <iostream>

namespace mps2D {

namespace {

// 原有的 LSMPS type-A 压力梯度离散（保留给自由面/飞溅粒子使用）
double2 ComputePressureGradientTypeAImpl(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>& moment_matrix_inverse,
    double smoothing_radius,
    double gravity_x,
    double gravity_y,
    double density) {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  double p_i = fluid_particles.pressure[particle_idx];
  
  // 提取压力corrective matrix的前两行（用于梯度计算，第二类边界条件）
  // 文档中对应 [M_{i,0}, M_{i,1}]
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> M0 = moment_matrix_inverse.row(0);  // x 方向
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> M1 = moment_matrix_inverse.row(1);  // y 方向
  
  double sum_x = 0.0;
  double sum_y = 0.0;
  
  CorrectiveMatrix corrective_matrix_calc;
  
  // 处理流体邻域粒子
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];
    double p_j = fluid_particles.pressure[j];
    
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    
    if (dist < 1e-10 || dist > smoothing_radius) {
      continue;
    }

    double dp = p_j - p_i;
    double weight = WeightFunction(dist, smoothing_radius);
    
    // 计算基函数（标准基函数 P_ij）
    Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis = 
        corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
    
    // 计算梯度贡献
    // 根据文档：∇φ_i = (1/r_s) Σ_j w_ij (φ_j - φ_i)
    //           [M_{i,0}; M_{i,1}] P_ij
    double M0P = (M0 * basis)(0, 0);
    double M1P = (M1 * basis)(0, 0);
    sum_x += weight * dp * M0P;
    sum_y += weight * dp * M1P;
  }
  
  // 处理壁面邻域粒子（第二类边界条件）
  // 壁面处压力梯度的法向分量：dp/dn = ρ * n · g
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    const double2& pos_j = solid_particles.position[j];
    const double2& normal = solid_particles.normal_vector[j];
    
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    
    if (dist < 1e-10 || dist > smoothing_radius) {
      continue;
    }
    
    double weight = WeightFunction(dist, smoothing_radius);
    
    // 使用壁面基函数（第二类边界条件）
    Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis_wall = 
        corrective_matrix_calc.ComputeBasisFunctionsForWall(
            dx, dy, normal.x, normal.y, smoothing_radius);
    
    // 壁面处压力梯度的法向分量：dp/dn = ρ * n · g
    // 归一化法向量
    double nn = std::sqrt(normal.x * normal.x + normal.y * normal.y);
    if (nn < 1e-10) {
      continue;  // 跳过无效的法向量
    }
    double n_x = normal.x / nn;
    double n_y = normal.y / nn;
    
    // 计算法向压力梯度：dp/dn = ρ * (n · g) = ρ * (n_x * g_x + n_y * g_y)
    double dp_dn = density * (n_x * gravity_x + n_y * gravity_y);
    
    // 计算梯度贡献
    // 文档中的壁面项：1/r_s Σ_{wall} w_ij r_s ρ g n [M_{i,0}; M_{i,1}] Q_ij
    double M0Q = (M0 * basis_wall)(0, 0);
    double M1Q = (M1 * basis_wall)(0, 0);
    sum_x += weight * smoothing_radius * dp_dn * M0Q;
    sum_y += weight * smoothing_radius * dp_dn * M1Q;
  }

  // 统一乘以 1 / r_s
  double inv_rs = 1.0 / smoothing_radius;
  double grad_x = inv_rs * sum_x;
  double grad_y = inv_rs * sum_y;

  return {grad_x, grad_y};
}

// Type-B LSMPS 压力梯度离散（用于内部粒子和近自由面粒子）
double2 ComputePressureGradientTypeBImpl(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>& moment_matrix_inverse,
    double smoothing_radius,
    double gravity_x,
    double gravity_y,
    double density) {

  const double2& pos_i = fluid_particles.position[particle_idx];

  // Type-B 使用 6 维基函数：p_ij, q_ij 参见 typeB.md
  constexpr int kBasisSizeB = 6;
  using MatrixB = Eigen::Matrix<double, kBasisSizeB, kBasisSizeB>;
  using VectorB = Eigen::Matrix<double, kBasisSizeB, 1>;
  using RowVectorB = Eigen::Matrix<double, 1, kBasisSizeB>;

  MatrixB moment = MatrixB::Zero();
  int neighbour_count = 0;

  if (smoothing_radius <= 0.0) {
    return ComputePressureGradientTypeAImpl(
        particle_idx, fluid_particles, solid_particles,
        moment_matrix_inverse, smoothing_radius,
        gravity_x, gravity_y, density);
  }

  const double inv_rs = 1.0 / smoothing_radius;
  const double inv_rs2 = inv_rs * inv_rs;

  // 构建流体粒子部分的矩矩阵 M_i^{fluid}
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];

    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    if (dist < 1e-10 || dist > smoothing_radius) {
      continue;
    }

    double weight = WeightFunction(dist, smoothing_radius);

    VectorB pij;
    pij[0] = 1.0;
    pij[1] = dx * inv_rs;
    pij[2] = dy * inv_rs;
    pij[3] = 0.5 * dx * dx * inv_rs2;
    pij[4] = 0.5 * dy * dy * inv_rs2;
    pij[5] = dx * dy * inv_rs2;

    moment += weight * pij * pij.transpose();
    ++neighbour_count;
  }

  // 构建壁面粒子部分的矩矩阵 M_i^{wall}
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    const double2& pos_j = solid_particles.position[j];
    const double2& normal = solid_particles.normal_vector[j];

    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    if (dist < 1e-10 || dist > smoothing_radius) {
      continue;
    }

    double nn = std::sqrt(normal.x * normal.x + normal.y * normal.y);
    if (nn < 1e-10) {
      continue;
    }
    double n_x = normal.x / nn;
    double n_y = normal.y / nn;

    double weight = WeightFunction(dist, smoothing_radius);

    VectorB qij;
    qij[0] = 0.0;
    qij[1] = n_x;
    qij[2] = n_y;
    qij[3] = dx * n_x * inv_rs;
    qij[4] = dy * n_y * inv_rs;
    qij[5] = (dx * n_y + dy * n_x) * inv_rs;

    moment += weight * qij * qij.transpose();
    ++neighbour_count;
  }

  // 邻域不足或矩阵不可逆时回退到 Type-A
  if (neighbour_count < kBasisSizeB) {
    return ComputePressureGradientTypeAImpl(
        particle_idx, fluid_particles, solid_particles,
        moment_matrix_inverse, smoothing_radius,
        gravity_x, gravity_y, density);
  }

  double det = moment.determinant();
  if (std::abs(det) < 1e-12) {
    return ComputePressureGradientTypeAImpl(
        particle_idx, fluid_particles, solid_particles,
        moment_matrix_inverse, smoothing_radius,
        gravity_x, gravity_y, density);
  }

  // LSMPS 矩阵 C_i = (M_i^{fluid} + M_i^{wall})^{-1}
  MatrixB C = moment.inverse();

  // 压力梯度对应的两行（φ_x, φ_y）
  RowVectorB Cx = C.row(1);
  RowVectorB Cy = C.row(2);

  double sum_fluid_x = 0.0;
  double sum_fluid_y = 0.0;
  double sum_wall_x = 0.0;
  double sum_wall_y = 0.0;

  // 流体邻域的压力梯度项：1/r_s Σ_j w_ij P_j [C_{i,2}; C_{i,3}] p_ij
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];
    double p_j = fluid_particles.pressure[j];

    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    if (dist < 1e-10 || dist > smoothing_radius) {
      continue;
    }

    double weight = WeightFunction(dist, smoothing_radius);

    VectorB pij;
    pij[0] = 1.0;
    pij[1] = dx * inv_rs;
    pij[2] = dy * inv_rs;
    pij[3] = 0.5 * dx * dx * inv_rs2;
    pij[4] = 0.5 * dy * dy * inv_rs2;
    pij[5] = dx * dy * inv_rs2;

    double coeff_x = Cx.dot(pij);
    double coeff_y = Cy.dot(pij);

    sum_fluid_x += weight * p_j * coeff_x;
    sum_fluid_y += weight * p_j * coeff_y;
  }

  // 壁面邻域的压力梯度项：Σ_j w_ij ρ (g · n_j) [C_{i,2}; C_{i,3}] q_ij
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    const double2& pos_j = solid_particles.position[j];
    const double2& normal = solid_particles.normal_vector[j];

    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    if (dist < 1e-10 || dist > smoothing_radius) {
      continue;
    }

    double nn = std::sqrt(normal.x * normal.x + normal.y * normal.y);
    if (nn < 1e-10) {
      continue;
    }
    double n_x = normal.x / nn;
    double n_y = normal.y / nn;

    double weight = WeightFunction(dist, smoothing_radius);

    double dp_dn = density * (n_x * gravity_x + n_y * gravity_y);

    VectorB qij;
    qij[0] = 0.0;
    qij[1] = n_x;
    qij[2] = n_y;
    qij[3] = dx * n_x * inv_rs;
    qij[4] = dy * n_y * inv_rs;
    qij[5] = (dx * n_y + dy * n_x) * inv_rs;

    double coeff_x = Cx.dot(qij);
    double coeff_y = Cy.dot(qij);

    sum_wall_x += weight * dp_dn * coeff_x;
    sum_wall_y += weight * dp_dn * coeff_y;
  }

  double grad_x = inv_rs * sum_fluid_x + sum_wall_x;
  double grad_y = inv_rs * sum_fluid_y + sum_wall_y;

  return {grad_x, grad_y};
}

}  // namespace

double2 Correction::ComputePressureGradient(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>& moment_matrix_inverse,
    double smoothing_radius,
    double gravity_x,
    double gravity_y,
    double density) {

  SurfaceType surface_type = SurfaceType::INNER;
  if (particle_idx >= 0 &&
      particle_idx < static_cast<int>(fluid_particles.surface_type.size())) {
    surface_type = fluid_particles.surface_type[particle_idx];
  }

  // 仅近自由面粒子使用 type-B，其余粒子使用原有 type-A
  if (surface_type == SurfaceType::NEAR_SURFACE) {
    return ComputePressureGradientTypeBImpl(
        particle_idx, fluid_particles, solid_particles,
        moment_matrix_inverse, smoothing_radius,
        gravity_x, gravity_y, density);
  }

  return ComputePressureGradientTypeAImpl(
      particle_idx, fluid_particles, solid_particles,
      moment_matrix_inverse, smoothing_radius,
      gravity_x, gravity_y, density);
}

double2 Correction::ComputeAcceleration(
    const double2& pressure_gradient,
    double density) {
  
  // 加速度 = -(压力梯度/密度)
  if (density < 1e-10) {
    std::cerr << "警告：密度过小，可能导致数值不稳定" << std::endl;
    return {0.0, 0.0};
  }
  
  return {-pressure_gradient.x / density, -pressure_gradient.y / density};
}

void Correction::UpdateVelocityAndPosition(
    int particle_idx,
    FluidParticle& fluid_particles,
    const double2& acceleration,
    double time_step) {
  
  // 更新速度：v_new = v_old + a * dt
  double2& velocity = fluid_particles.velocity[particle_idx];
  velocity.x += acceleration.x * time_step;
  velocity.y += acceleration.y * time_step;
  
  // 更新位置：x_new = x_old + v_new * dt
  double2& position = fluid_particles.position[particle_idx];
  position.x += velocity.x * time_step;
  position.y += velocity.y * time_step;
}

void Correction::ComputeAndUpdateAllParticles(
    FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>>& corrective_matrices,
    double smoothing_radius,
    double gravity_x,
    double gravity_y,
    double density,
    double time_step) {
  
  int num_particles = fluid_particles.particle_num;
  
  if (static_cast<int>(corrective_matrices.size()) != num_particles) {
    std::cerr << "错误：corrective_matrices数量与流体粒子数不匹配" << std::endl;
    return;
  }
  
  // 第一步：先计算所有粒子的压力梯度和加速度（使用原始位置）
  std::vector<double2> accelerations(num_particles);
  for (int i = 0; i < num_particles; ++i) {
    // 计算压力梯度（此时所有粒子的位置都还是原始值）
    double2 pressure_gradient = ComputePressureGradient(
        i,
        fluid_particles,
        solid_particles,
        corrective_matrices[i],
        smoothing_radius,
        gravity_x,
        gravity_y,
        density);
    
    // 计算加速度
    accelerations[i] = ComputeAcceleration(pressure_gradient, density);
  }
  
  // 第二步：统一更新所有粒子的速度和位置
  for (int i = 0; i < num_particles; ++i) {
    UpdateVelocityAndPosition(i, fluid_particles, accelerations[i], time_step);
  }
}

void Correction::ComputeAndUpdateAllParticles(
    FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>>& corrective_matrices,
    double smoothing_radius,
    double gravity_x,
    double gravity_y,
    double density,
    double time_step,
    std::vector<double2>& pressure_gradients) {
  
  int num_particles = fluid_particles.particle_num;
  
  if (static_cast<int>(corrective_matrices.size()) != num_particles) {
    std::cerr << "错误：corrective_matrices数量与流体粒子数不匹配" << std::endl;
    return;
  }
  
  // 确保输出向量大小正确
  pressure_gradients.resize(num_particles);
  
  // 第一步：先计算所有粒子的压力梯度和加速度（使用原始位置）
  std::vector<double2> accelerations(num_particles);
  for (int i = 0; i < num_particles; ++i) {
    // 计算压力梯度（此时所有粒子的位置都还是原始值）
    double2 pressure_gradient = ComputePressureGradient(
        i,
        fluid_particles,
        solid_particles,
        corrective_matrices[i],
        smoothing_radius,
        gravity_x,
        gravity_y,
        density);
    
    // 保存压力梯度（在位置更新之前）
    pressure_gradients[i] = pressure_gradient;
    
    // 计算加速度
    accelerations[i] = ComputeAcceleration(pressure_gradient, density);
  }
  
  // 第二步：统一更新所有粒子的速度和位置
  for (int i = 0; i < num_particles; ++i) {
    UpdateVelocityAndPosition(i, fluid_particles, accelerations[i], time_step);
  }
}

} // namespace mps2D

