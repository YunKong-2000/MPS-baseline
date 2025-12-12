#include "CorrectiveMatrix.hpp"
#include <cmath>

namespace mps2D {

Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>
CorrectiveMatrix::ComputeCorrectiveMatrix(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double smoothing_radius) {
  
  // 统计邻域粒子总数
  int num_fluid_neighbors = fluid_particles.fluid_neighbour_list[particle_idx].size();
  int num_solid_neighbors = fluid_particles.solid_neighbour_list[particle_idx].size();
  int total_neighbors = num_fluid_neighbors + num_solid_neighbors;
  
  // 如果邻域粒子数不足，返回单位矩阵
  // 至少需要5个邻域粒子才能求解5x5系统
  if (total_neighbors < MATRIX_SIZE) {
    return Eigen::Matrix<double, MATRIX_SIZE, MATRIX_SIZE>::Identity();
  }
  
  // 构建系数矩阵A和权重矩阵W
  Eigen::MatrixXd C = BuildCoefficientMatrix(
      particle_idx, fluid_particles, solid_particles, smoothing_radius);
    
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
    double smoothing_radius) {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  
  // 创建系数矩阵C (5 x 5)
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(BASIS_SIZE, BASIS_SIZE);
  
  // 处理流体邻域粒子
  // 对每个j流体邻域粒子，计算权重和基函数，然后累加到C中
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    double weight = WeightFunction(dist, smoothing_radius);
    Eigen::Vector<double, BASIS_SIZE> basis = ComputeBasisFunctions(dx, dy, smoothing_radius);
    C += weight * basis * basis.transpose();
  }
  
  // 处理固体邻域粒子（壁面粒子）
  // 对每个j固体邻域粒子，使用壁面基函数计算方法
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    const double2& pos_j = solid_particles.position[j];
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    double weight = WeightFunction(dist, smoothing_radius);
    
    // 获取壁面法向量
    const double2& normal = solid_particles.normal_vector[j];
    Eigen::Vector<double, BASIS_SIZE> basis = ComputeBasisFunctionsForWall(
        dx, dy, normal.x, normal.y, smoothing_radius);
    C += weight * basis * basis.transpose();
  }
  
  return C;
}

Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE>
CorrectiveMatrix::ComputeBasisFunctions(
    double dx, double dy, double smoothing_radius) {
  
  // 归一化相对位置
  double dist = std::sqrt(dx * dx + dy * dy);
  
  // 基函数：x/r_e, y/r_e, (x/r_e)^2, (y/r_e)^2, (x/r_e)*(y/r_e)
  Eigen::Vector<double, BASIS_SIZE> basis;
  basis[0] = dx / dist;
  basis[1] = dy / dist;
  basis[2] = dx * dx / (dist * smoothing_radius);
  basis[3] = dy * dy / (dist * smoothing_radius);
  basis[4] = dx * dy / (dist * smoothing_radius);
  
  return basis;
}

Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE>
CorrectiveMatrix::ComputeBasisFunctionsForWall(
    double dx, double dy,
    double normal_x, double normal_y,
    double smoothing_radius) {
  
  // 归一化相对位置
  double x_norm = dx / smoothing_radius;
  double y_norm = dy / smoothing_radius;
  
  // 壁面基函数：n_x, n_y, 2*n_x*x/r_e, 2*n_y*y/r_e, (n_x*x + n_y*y)/r_e
  Eigen::Vector<double, BASIS_SIZE> basis;
  basis[0] = normal_x;                                    // n_x
  basis[1] = normal_y;                                    // n_y
  basis[2] = 2.0 * normal_x * x_norm;                     // 2*n_x*x/r_e
  basis[3] = 2.0 * normal_y * y_norm;                     // 2*n_y*y/r_e
  basis[4] = (normal_x * x_norm + normal_y * y_norm);     // (n_x*x + n_y*y)/r_e
  
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

} // namespace mps2D

