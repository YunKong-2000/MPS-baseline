#include "PPEMatrixBuilder.hpp"
#include "../lsmps/CorrectiveMatrix.hpp"
#include <cmath>
#include <iostream>

namespace mps2D {

bool PPEMatrixBuilder::BuildPPEMatrix(
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>>& corrective_matrices,
    double smoothing_radius,
    double density,
    double time_step,
      double gravity_x,
      double gravity_y,
      Eigen::SparseMatrix<double, Eigen::ColMajor>& A,
      Eigen::VectorXd& b) {
  
  int num_fluid_particles = fluid_particles.particle_num;
  
  if (num_fluid_particles == 0) {
    std::cerr << "错误：流体粒子数为0" << std::endl;
    return false;
  }
  
  if (static_cast<int>(corrective_matrices.size()) != num_fluid_particles) {
    std::cerr << "错误：corrective matrix数量与流体粒子数不匹配" << std::endl;
    return false;
  }
  
  // 统计非零元素个数并预分配内存
  int nnz = CountNonZeros(fluid_particles);
  
  // 初始化稀疏矩阵（使用三元组形式构建）
  // 使用三元组格式可以高效地构建稀疏矩阵，然后转换为压缩格式（CSC）
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(nnz);
  
  // 初始化右边项向量
  b = Eigen::VectorXd::Zero(num_fluid_particles);
  
  // 遍历所有流体粒子，构建系数矩阵和右边项
  for (int i = 0; i < num_fluid_particles; ++i) {
    BuildParticleRow(
        i, fluid_particles, solid_particles, corrective_matrices[i],
        smoothing_radius, density, time_step, gravity_x, gravity_y,
        triplets, b);
  }
  
  // 构建稀疏矩阵
  // 使用ColMajor（列主序）格式，这是CSC（Compressed Sparse Column）格式
  // 对于迭代求解器，CSC格式通常更高效
  A.resize(num_fluid_particles, num_fluid_particles);
  
  // 对三元组按列优先排序，可以提高压缩效率
  // 排序规则：先按列索引排序，再按行索引排序
  std::sort(triplets.begin(), triplets.end(), 
            [](const Eigen::Triplet<double>& a, const Eigen::Triplet<double>& b) {
              if (a.col() != b.col()) {
                return a.col() < b.col();
              }
              return a.row() < b.row();
            });
  
  // 从三元组构建稀疏矩阵
  // setFromTriplets会自动处理重复元素（累加）
  A.setFromTriplets(triplets.begin(), triplets.end());
  
  // 转换为压缩格式（CSC）
  // 压缩格式可以显著减少内存占用和提高矩阵-向量乘法性能
  A.makeCompressed();
  
  return true;
}

int PPEMatrixBuilder::CountNonZeros(const FluidParticle& fluid_particles) const {
  int nnz = 0;
  int num_fluid_particles = fluid_particles.particle_num;
  
  // 每个粒子至少有一个对角线元素
  nnz += num_fluid_particles;
  
  // 加上所有邻域粒子对应的非对角线元素
  for (int i = 0; i < num_fluid_particles; ++i) {
    nnz += fluid_particles.fluid_neighbour_list[i].size();
  }
  
  return nnz;
}

void PPEMatrixBuilder::BuildParticleRow(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>& corrective_matrix,
    double smoothing_radius,
    double density,
    double time_step,
    double gravity_x,
    double gravity_y,
    std::vector<Eigen::Triplet<double>>& triplets,
    Eigen::VectorXd& b) const {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  const double2& vel_i = fluid_particles.velocity[particle_idx];
  
  // 创建CorrectiveMatrix实例用于计算基函数
  CorrectiveMatrix corrective_matrix_calc;
  
  // 提取corrective matrix的行向量
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C1 = corrective_matrix.row(0);
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C2 = corrective_matrix.row(1);
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C3 = corrective_matrix.row(2);
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C4 = corrective_matrix.row(3);
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C3_plus_C4 = C3 + C4;
  
  // 系数矩阵的系数因子
  double coeff_factor = 2.0 / (smoothing_radius * density);
  
  // 初始化对角线系数累加和右边项
  double diag_sum = 0.0;
  double divergence = 0.0;
  double wall_pressure_term = 0.0;  // 壁面压力边界条件项
  
  // 遍历流体邻域粒子，同时计算对角线系数、非对角线系数和右边项
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];
    const double2& vel_j = fluid_particles.velocity[j];
    
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    
    if (dist < 1e-10 || dist > smoothing_radius) {
      continue;
    }
    
    double weight = WeightFunction(dist, smoothing_radius);
    
    // 使用CorrectiveMatrix的方法计算基函数（第一类边界条件，用于流体粒子）
    Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis = 
        corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
    
    // 计算公共项：(w_ij / r_ij) * [C_3 + C_4] * P_ij（用于压力拉普拉斯算子）
    // 注意：C3_plus_C4 * basis 返回1x1矩阵，使用(0, 0)访问标量值
    double common_coeff = (weight / dist) * (C3_plus_C4 * basis)(0, 0);
    
    // 累加对角线系数（负号）
    diag_sum += common_coeff;
    
    // 添加非对角线系数（正号）
    double off_diag_coeff = coeff_factor * common_coeff;
    triplets.push_back(Eigen::Triplet<double>(particle_idx, j, off_diag_coeff));
    
    // 计算速度散度项（用于右边项，使用第一类边界条件的基函数）
    double dux_dr = (vel_j.x - vel_i.x) / dist;
    double duy_dr = (vel_j.y - vel_i.y) / dist;
    double C1P = (C1 * basis)(0, 0);
    double C2P = (C2 * basis)(0, 0);
    divergence += weight * (C1P * dux_dr + C2P * duy_dr);
  }
  
  // 处理壁面邻域粒子（近壁面内部粒子）
  // 注意：压力拉普拉斯算子使用第二类边界条件，速度散度使用第一类边界条件
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    const double2& pos_j = solid_particles.position[j];
    const double2& vel_wall = solid_particles.velocity[j];  // 壁面速度（通常为零）
    const double2& normal = solid_particles.normal_vector[j];  // 壁面法向量
    
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    
    if (dist < 1e-10 || dist > smoothing_radius) {
      continue;
    }
    
    double weight = WeightFunction(dist, smoothing_radius);
    
    // 对于速度散度（右边项）：使用第一类边界条件的基函数（标准基函数）
    Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis_velocity = 
        corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
    
    // 对于壁面压力边界条件项（右边项）：使用第二类边界条件的基函数（壁面基函数）
    Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis_pressure = 
        corrective_matrix_calc.ComputeBasisFunctionsForWall(
            dx, dy, normal.x, normal.y, smoothing_radius);
    
    // 计算速度散度项（壁面速度对右边项的贡献，使用第一类边界条件）
    // 根据文档公式：Σ(w_ij * (u_wall - u_i)/r_ij * [C_1; C_2] * P_ij)
    double dux_dr = (vel_wall.x - vel_i.x) / dist;
    double duy_dr = (vel_wall.y - vel_i.y) / dist;
    double C1P = (C1 * basis_velocity)(0, 0);
    double C2P = (C2 * basis_velocity)(0, 0);
    divergence += weight * (C1P * dux_dr + C2P * duy_dr);
    
    // 计算壁面压力边界条件项（使用第二类边界条件的基函数）
    double n_dot_g = normal.x * gravity_x + normal.y * gravity_y;
    double wall_pressure_coeff = -weight * (density * n_dot_g) * (C3_plus_C4 * basis_pressure)(0, 0);
    wall_pressure_term += wall_pressure_coeff;
  }
  
  // 设置对角线系数
  double diag_coeff = -coeff_factor * diag_sum;
  triplets.push_back(Eigen::Triplet<double>(particle_idx, particle_idx, diag_coeff));
  
  // 设置右边项
  b(particle_idx) = (1.0 / time_step) * divergence + coeff_factor * wall_pressure_term;
}

bool PPEMatrixBuilder::BuildPPEMatrixWithDebug(
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>>& corrective_matrices,
    double smoothing_radius,
    double density,
    double time_step,
    double gravity_x,
    double gravity_y,
    Eigen::SparseMatrix<double, Eigen::ColMajor>& A,
    Eigen::VectorXd& b,
    DebugInfo& debug_info) {
  
  int num_fluid_particles = fluid_particles.particle_num;
  
  if (num_fluid_particles == 0) {
    std::cerr << "错误：流体粒子数为0" << std::endl;
    return false;
  }
  
  if (static_cast<int>(corrective_matrices.size()) != num_fluid_particles) {
    std::cerr << "错误：corrective matrix数量与流体粒子数不匹配" << std::endl;
    return false;
  }
  
  // 统计非零元素个数并预分配内存
  int nnz = CountNonZeros(fluid_particles);
  
  // 初始化稀疏矩阵（使用三元组形式构建）
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(nnz);
  
  // 初始化右边项向量
  b = Eigen::VectorXd::Zero(num_fluid_particles);
  
  // 初始化调试信息
  debug_info.diagonal_coefficients = Eigen::VectorXd::Zero(num_fluid_particles);
  debug_info.off_diagonal_row_sums = Eigen::VectorXd::Zero(num_fluid_particles);
  debug_info.divergence_terms = Eigen::VectorXd::Zero(num_fluid_particles);
  debug_info.wall_pressure_terms = Eigen::VectorXd::Zero(num_fluid_particles);
  debug_info.right_hand_side = Eigen::VectorXd::Zero(num_fluid_particles);
  
  // 遍历所有流体粒子，构建系数矩阵和右边项
  for (int i = 0; i < num_fluid_particles; ++i) {
    BuildParticleRowWithDebug(
        i, fluid_particles, solid_particles, corrective_matrices[i],
        smoothing_radius, density, time_step, gravity_x, gravity_y,
        triplets, b, debug_info);
  }
  
  // 构建稀疏矩阵
  A.resize(num_fluid_particles, num_fluid_particles);
  
  // 对三元组按列优先排序
  std::sort(triplets.begin(), triplets.end(), 
            [](const Eigen::Triplet<double>& a, const Eigen::Triplet<double>& b) {
              if (a.col() != b.col()) {
                return a.col() < b.col();
              }
              return a.row() < b.row();
            });
  
  A.setFromTriplets(triplets.begin(), triplets.end());
  A.makeCompressed();
  
  // 保存完整的右边项
  debug_info.right_hand_side = b;
  
  return true;
}

void PPEMatrixBuilder::BuildParticleRowWithDebug(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>& corrective_matrix,
    double smoothing_radius,
    double density,
    double time_step,
    double gravity_x,
    double gravity_y,
    std::vector<Eigen::Triplet<double>>& triplets,
    Eigen::VectorXd& b,
    DebugInfo& debug_info) const {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  const double2& vel_i = fluid_particles.velocity[particle_idx];
  
  // 创建CorrectiveMatrix实例用于计算基函数
  CorrectiveMatrix corrective_matrix_calc;
  
  // 提取corrective matrix的行向量
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C1 = corrective_matrix.row(0);
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C2 = corrective_matrix.row(1);
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C3 = corrective_matrix.row(2);
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C4 = corrective_matrix.row(3);
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C3_plus_C4 = C3 + C4;
  
  // 系数矩阵的系数因子
  double coeff_factor = 2.0 / (smoothing_radius * density);
  
  // 初始化对角线系数累加和右边项
  double diag_sum = 0.0;
  double divergence = 0.0;
  double wall_pressure_term = 0.0;  // 壁面压力边界条件项
  double off_diag_sum = 0.0;  // 非对角线系数之和（用于调试）
  
  // 遍历流体邻域粒子，同时计算对角线系数、非对角线系数和右边项
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];
    const double2& vel_j = fluid_particles.velocity[j];
    
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    
    if (dist < 1e-10 || dist > smoothing_radius) {
      continue;
    }
    
    double weight = WeightFunction(dist, smoothing_radius);
    
    // 使用CorrectiveMatrix的方法计算基函数（第一类边界条件，用于流体粒子）
    Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis = 
        corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
    
    // 计算公共项：(w_ij / r_ij) * [C_3 + C_4] * P_ij（用于压力拉普拉斯算子）
    double common_coeff = (weight / dist) * (C3_plus_C4 * basis)(0, 0);
    
    // 累加对角线系数（负号）
    diag_sum += common_coeff;
    
    // 添加非对角线系数（正号）
    double off_diag_coeff = coeff_factor * common_coeff;
    triplets.push_back(Eigen::Triplet<double>(particle_idx, j, off_diag_coeff));
    off_diag_sum += off_diag_coeff;  // 累加非对角线系数
    
    // 计算速度散度项（用于右边项，使用第一类边界条件的基函数）
    double dux_dr = (vel_j.x - vel_i.x) / dist;
    double duy_dr = (vel_j.y - vel_i.y) / dist;
    double C1P = (C1 * basis)(0, 0);
    double C2P = (C2 * basis)(0, 0);
    divergence += weight * (C1P * dux_dr + C2P * duy_dr);
  }
  
  // 处理壁面邻域粒子（近壁面内部粒子）
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    const double2& pos_j = solid_particles.position[j];
    const double2& vel_wall = solid_particles.velocity[j];
    const double2& normal = solid_particles.normal_vector[j];
    
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    
    if (dist < 1e-10 || dist > smoothing_radius) {
      continue;
    }
    
    double weight = WeightFunction(dist, smoothing_radius);
    
    // 对于速度散度（右边项）：使用第一类边界条件的基函数（标准基函数）
    Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis_velocity = 
        corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
    
    // 对于壁面压力边界条件项（右边项）：使用第二类边界条件的基函数（壁面基函数）
    Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis_pressure = 
        corrective_matrix_calc.ComputeBasisFunctionsForWall(
            dx, dy, normal.x, normal.y, smoothing_radius);
    
    // 计算速度散度项（壁面速度对右边项的贡献，使用第一类边界条件）
    double dux_dr = (vel_wall.x - vel_i.x) / dist;
    double duy_dr = (vel_wall.y - vel_i.y) / dist;
    double C1P = (C1 * basis_velocity)(0, 0);
    double C2P = (C2 * basis_velocity)(0, 0);
    divergence += weight * (C1P * dux_dr + C2P * duy_dr);
    
    // 计算壁面压力边界条件项（使用第二类边界条件的基函数）
    double n_dot_g = normal.x * gravity_x + normal.y * gravity_y;
    double wall_pressure_coeff = -weight * (density * n_dot_g) * (C3_plus_C4 * basis_pressure)(0, 0);
    wall_pressure_term += wall_pressure_coeff;
  }
  
  // 设置对角线系数
  double diag_coeff = -coeff_factor * diag_sum;
  triplets.push_back(Eigen::Triplet<double>(particle_idx, particle_idx, diag_coeff));
  
  // 保存调试信息
  debug_info.diagonal_coefficients(particle_idx) = diag_coeff;
  debug_info.off_diagonal_row_sums(particle_idx) = off_diag_sum;
  debug_info.divergence_terms(particle_idx) = divergence;  // 未除以时间步长
  debug_info.wall_pressure_terms(particle_idx) = wall_pressure_term;  // 未乘以系数因子
  
  // 设置右边项
  b(particle_idx) = (1.0 / time_step) * divergence + coeff_factor * wall_pressure_term;
}

} // namespace mps2D
