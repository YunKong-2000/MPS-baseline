#include "CorrectiveMatrix.hpp"
#include <cmath>

namespace mps2D {

namespace {
constexpr double kSurfaceOperatorDiagonalRegularization = 1e-3;

bool NeedsSurfaceRegularization(
    int particle_idx,
    const FluidParticle& fluid_particles) {
  if (particle_idx < 0 ||
      particle_idx >= static_cast<int>(fluid_particles.surface_type.size())) {
    return false;
  }
  const SurfaceType surface_type = fluid_particles.surface_type[particle_idx];
  return surface_type == SurfaceType::SURFACE ||
         surface_type == SurfaceType::NEAR_SURFACE;
}
}  // namespace

Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>
CorrectiveMatrix::ComputeCorrectiveMatrix(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double smoothing_radius,
    bool border_condition) {
  // 兼容旧接口：直接复用moment matrix的计算
  return ComputeMomentMatrix(particle_idx, fluid_particles, solid_particles,
                             smoothing_radius, border_condition);
}

Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>
CorrectiveMatrix::ComputeCorrectiveMatrixFluidOnly(
    int particle_idx,
    const FluidParticle& fluid_particles,
    double smoothing_radius) {
  // 兼容旧接口：直接复用仅流体粒子的moment matrix计算
  return ComputeMomentMatrixFluidOnly(particle_idx, fluid_particles,
                                      smoothing_radius);
}

Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>
CorrectiveMatrix::ComputeMomentMatrix(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double smoothing_radius,
    bool border_condition) {
  
  // 统计邻域粒子总数
  int num_fluid_neighbors = fluid_particles.fluid_neighbour_list[particle_idx].size();
  int num_solid_neighbors = fluid_particles.solid_neighbour_list[particle_idx].size();
  int total_neighbors = num_fluid_neighbors + num_solid_neighbors;
  
  // 如果邻域粒子数不足，返回单位矩阵
  // 至少需要5个邻域粒子才能求解5x5系统
  if (total_neighbors < MATRIX_SIZE) {
    return Eigen::Matrix<double, MATRIX_SIZE, MATRIX_SIZE>::Identity();
  }
  
  // 构建moment矩矩阵M
  Eigen::MatrixXd M = BuildCoefficientMatrix(
      particle_idx, fluid_particles, solid_particles, smoothing_radius, border_condition);

  if (NeedsSurfaceRegularization(particle_idx, fluid_particles)) {
    // 自由面与近自由面粒子的算子矩阵在求逆前统一加对角正则，提升稳定性。
    M.diagonal().array() += kSurfaceOperatorDiagonalRegularization;
  }

  // 检查矩阵是否可逆
  if (!IsMatrixInvertible(M)) {
    return Eigen::Matrix<double, MATRIX_SIZE, MATRIX_SIZE>::Identity();
  }
  
  // 计算moment matrix的逆矩阵: M^(-1)
  Eigen::Matrix<double, MATRIX_SIZE, MATRIX_SIZE> moment_matrix = M.inverse();
  return moment_matrix;
}

Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>
CorrectiveMatrix::ComputeMomentMatrixFluidOnly(
    int particle_idx,
    const FluidParticle& fluid_particles,
    double smoothing_radius) {
  
  // 统计流体邻域粒子数
  int num_fluid_neighbors = fluid_particles.fluid_neighbour_list[particle_idx].size();
      
  // 如果流体邻域粒子数不足，返回单位矩阵
  // 至少需要5个邻域粒子才能求解5x5系统
  if (num_fluid_neighbors < MATRIX_SIZE) {
    return Eigen::Matrix<double, MATRIX_SIZE, MATRIX_SIZE>::Identity();
  }
  
  // 构建系数矩阵C（仅考虑流体粒子）
  const double2& pos_i = fluid_particles.position[particle_idx];
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(BASIS_SIZE, BASIS_SIZE);
  
  // 只处理流体邻域粒子
  // 对每个j流体邻域粒子，计算权重和基函数，然后累加到C中
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    if (dist < 1e-10 || dist > smoothing_radius) continue;
    double weight = WeightFunction(dist, smoothing_radius);
    Eigen::Vector<double, BASIS_SIZE> basis = ComputeBasisFunctions(dx, dy, smoothing_radius);
    C += weight * basis * basis.transpose();
  }

  if (NeedsSurfaceRegularization(particle_idx, fluid_particles)) {
    // 仅流体邻域版本与主流程保持一致，对表面相关粒子统一施加对角正则。
    C.diagonal().array() += kSurfaceOperatorDiagonalRegularization;
  }
  
  // 检查矩阵是否可逆
  if (!IsMatrixInvertible(C)) {
    return Eigen::Matrix<double, MATRIX_SIZE, MATRIX_SIZE>::Identity();
  }
  
  // 计算corrective matrix: C^(-1)
  Eigen::Matrix<double, MATRIX_SIZE, MATRIX_SIZE> corrective_matrix = C.inverse();
  return corrective_matrix;
}

Eigen::MatrixXd
CorrectiveMatrix::BuildCoefficientMatrix(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double smoothing_radius,
    bool border_condition) {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  
  // 创建系数矩阵C (5 x 5)
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(BASIS_SIZE, BASIS_SIZE);
  
  // 处理流体邻域粒子
  // 对每个j流体邻域粒子，计算权重和基函数，然后累加到C中
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    if (dist < 1e-10 || dist > smoothing_radius) continue;
    double weight = WeightFunction(dist, smoothing_radius);
    Eigen::Vector<double, BASIS_SIZE> basis = ComputeBasisFunctions(dx, dy, smoothing_radius);
    M += weight * basis * basis.transpose();
  }
  
  // 处理固体邻域粒子（壁面粒子）
  // 对每个j固体邻域粒子，使用壁面基函数计算方法
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    const double2& pos_j = solid_particles.position[j];
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    if (dist < 1e-10 || dist > smoothing_radius) continue;
    double weight = WeightFunction(dist, smoothing_radius);
    
    // 获取壁面法向量
    const double2& normal = solid_particles.normal_vector[j];
    Eigen::Vector<double, BASIS_SIZE> basis;
    if (border_condition) {
      basis = ComputeBasisFunctionsForWall(
          dx, dy, normal.x, normal.y, smoothing_radius);
    } else {
      basis = ComputeBasisFunctions(dx, dy, smoothing_radius);
    }
    M += weight * basis * basis.transpose();
  }
  
  return M;
}

Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE>
CorrectiveMatrix::ComputeBasisFunctions(
    double dx, double dy, double smoothing_radius) {
  // lsmps使用统一的r_s进行归一化，选取r_s=r_e
  // 基函数：x/r_e, y/r_e, (x/r_e)^2, (y/r_e)^2, (x/r_e)*(y/r_e)
  Eigen::Vector<double, BASIS_SIZE> basis;
  basis[0] = dx / smoothing_radius;
  basis[1] = dy / smoothing_radius;
  basis[2] = dx * dx / (smoothing_radius * smoothing_radius);
  basis[3] = dy * dy / (smoothing_radius * smoothing_radius);
  basis[4] = dx * dy / (smoothing_radius * smoothing_radius);
  return basis;
}

Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE>
CorrectiveMatrix::ComputeBasisFunctionsForWall(
    double dx, double dy,
    double normal_x, double normal_y,
    double smoothing_radius) {
  // 壁面基函数：n_x, n_y, 2*n_x*x/r_e, 2*n_y*y/r_e, (n_x*dy + n_y*dx)/r_e
  // 与 manual 测试保持一致
  Eigen::Vector<double, BASIS_SIZE> basis;
  basis[0] = normal_x;                                    // n_x
  basis[1] = normal_y;                                    // n_y
  basis[2] = 2.0 * normal_x * dx / smoothing_radius;      // 2*n_x*x/r_e
  basis[3] = 2.0 * normal_y * dy / smoothing_radius;      // 2*n_y*y/r_e
  basis[4] = (normal_x * dy + normal_y * dx) / smoothing_radius;  // (n_x*dy + n_y*dx)/r_e
  return basis;
}

bool
CorrectiveMatrix::IsMatrixInvertible(
    const Eigen::Matrix<double, MATRIX_SIZE, MATRIX_SIZE>& matrix,
    double tolerance) const {
  
  // 使用行列式判断矩阵是否可逆
  double determinant = matrix.determinant();
  return std::abs(determinant) > tolerance;
}

Eigen::MatrixXd
CorrectiveMatrix::BuildCoefficientMatrixForDiagnostics(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double smoothing_radius,
    bool border_condition) {
  // 直接调用私有的BuildCoefficientMatrix方法
  return BuildCoefficientMatrix(particle_idx, fluid_particles, solid_particles,
                                smoothing_radius, border_condition);
}

} // namespace mps2D

