#include "PPEMatrixBuilder.hpp"
#include "../lsmps/CorrectiveMatrix.hpp"
#include "../core/FileOperator.hpp"
#include "../mps/OriginalMPS.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <filesystem>

// PETSc头文件（用于直接构建PETSc矩阵）
#include <petsc.h>
#include <petscmat.h>
#include <petscvec.h>

namespace mps2D {

// 静态标志，确保PETSc只初始化一次
static bool petsc_initialized_in_builder = false;

// 确保PETSc已初始化（用于BuildPPEMatrixPetsc）
static void EnsurePetscInitialized() {
  if (!petsc_initialized_in_builder) {
    // 检查PETSc是否已经初始化
    PetscBool initialized = PETSC_FALSE;
    PetscInitialized(&initialized);
    
    if (!initialized) {
      int argc = 0;
      char** argv = nullptr;
      PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr);
      if (ierr) {
        std::cerr << "错误：PETSc初始化失败" << std::endl;
        return;
      }
      
      // 设置PETSc选项：不显示版权信息
      PetscOptionsSetValue(nullptr, "-options_left", "false");
    }
    petsc_initialized_in_builder = true;
  }
}

bool PPEMatrixBuilder::BuildPPEMatrixPetsc(
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>>& corrective_matrices_velocity,
    const std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>>& corrective_matrices_pressure,
    double smoothing_radius,
    double density,
    double time_step,
    double particle_spacing,
    double gravity_x,
    double gravity_y,
    Mat& A_petsc,
    Vec& b_petsc) {
  
  // 确保PETSc已初始化
  EnsurePetscInitialized();
  
  int num_fluid_particles = fluid_particles.particle_num;
  
  if (num_fluid_particles == 0) {
    std::cerr << "错误：流体粒子数为0" << std::endl;
    return false;
  }
  
  if (static_cast<int>(corrective_matrices_velocity.size()) != num_fluid_particles) {
    std::cerr << "错误：速度corrective matrix数量与流体粒子数不匹配" << std::endl;
    return false;
  }
  
  if (static_cast<int>(corrective_matrices_pressure.size()) != num_fluid_particles) {
    std::cerr << "错误：压力corrective matrix数量与流体粒子数不匹配" << std::endl;
    return false;
  }
  
  // 统计每行的非零元素数
  std::vector<int> nnz_per_row(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
    nnz_per_row[i] = 1 + static_cast<int>(fluid_particles.fluid_neighbour_list[i].size());
  }
  
  // 初始化PETSc矩阵和向量
  InitializePetscMatrixAndVector(num_fluid_particles, nnz_per_row, A_petsc, b_petsc);
  
  // 创建CorrectiveMatrix实例用于计算基函数
  CorrectiveMatrix corrective_matrix_calc;
  
  // 系数矩阵的系数因子（统一用于所有粒子）
  double coeff_factor = 2.0 / (smoothing_radius * density);
  
  // 遍历所有流体粒子，构建系数矩阵和右边项
  // 统一使用LSMPS离散方法（所有粒子，包括自由面粒子）
  for (int particle_idx = 0; particle_idx < num_fluid_particles; ++particle_idx) {
    // 所有粒子都使用相同的LSMPS离散方法
    // 速度散度使用第一类边界条件的corrective matrix
    // 压力拉普拉斯算子使用第二类边界条件的corrective matrix
      BuildInnerParticleRow(
        particle_idx, fluid_particles, solid_particles, 
        corrective_matrices_velocity[particle_idx],
        corrective_matrices_pressure[particle_idx],
          corrective_matrix_calc, smoothing_radius, density, time_step,
          gravity_x, gravity_y, coeff_factor, A_petsc, b_petsc);
  }
  
  // 组装矩阵和向量
  MatAssemblyBegin(A_petsc, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A_petsc, MAT_FINAL_ASSEMBLY);
  
  VecAssemblyBegin(b_petsc);
  VecAssemblyEnd(b_petsc);
  
  return true;
}

bool PPEMatrixBuilder::BuildPenaltySystem(
    Mat A_petsc,
    Vec b_petsc,
    const FluidParticle& fluid_particles,
    double penalty_parameter,
    Mat& K_petsc,
    Vec& f_petsc) const {
  
  if (A_petsc == NULL || b_petsc == NULL) {
    std::cerr << "错误：输入矩阵或向量为NULL" << std::endl;
    return false;
  }
  
  int num_fluid_particles = fluid_particles.particle_num;
  if (num_fluid_particles == 0) {
    std::cerr << "错误：流体粒子数为0" << std::endl;
    return false;
  }
  
  // 确保PETSc已初始化
  EnsurePetscInitialized();
  
  // 步骤1：计算 A^T
  Mat A_transpose = NULL;
  MatTranspose(A_petsc, MAT_INITIAL_MATRIX, &A_transpose);
  
  // 步骤2：计算 K = A^T A
  MatMatMult(A_transpose, A_petsc, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &K_petsc);
  
  // 步骤3：构建对角惩罚矩阵 D 并添加到 K
  // D_{ii} = μ (如果 i 是自由面粒子)，否则 0
  Vec diagonal_penalty = NULL;
  VecCreate(PETSC_COMM_WORLD, &diagonal_penalty);
  VecSetSizes(diagonal_penalty, PETSC_DECIDE, num_fluid_particles);
  VecSetType(diagonal_penalty, VECSEQ);
  VecSetFromOptions(diagonal_penalty);
  VecSet(diagonal_penalty, 0.0);
  
  // 设置自由面粒子的惩罚项
  std::vector<PetscInt> surface_indices;
  for (int i = 0; i < num_fluid_particles; ++i) {
    if (fluid_particles.surface_type[i] == SurfaceType::SURFACE) {
      surface_indices.push_back(static_cast<PetscInt>(i));
    }
  }
  
  if (!surface_indices.empty()) {
    std::vector<PetscScalar> penalty_values(surface_indices.size(), penalty_parameter);
    VecSetValues(diagonal_penalty, surface_indices.size(), surface_indices.data(), 
                 penalty_values.data(), INSERT_VALUES);
    VecAssemblyBegin(diagonal_penalty);
    VecAssemblyEnd(diagonal_penalty);
  }
  
  // 将惩罚矩阵 D 添加到 K 的对角线上
  // K = K + D（其中 D 是对角矩阵）
  MatDiagonalSet(K_petsc, diagonal_penalty, ADD_VALUES);
  
  // 组装矩阵
  MatAssemblyBegin(K_petsc, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K_petsc, MAT_FINAL_ASSEMBLY);
  
  // 标记矩阵为对称正定（K = A^T A + D 是对称正定矩阵）
  MatSetOption(K_petsc, MAT_SYMMETRIC, PETSC_TRUE);
  MatSetOption(K_petsc, MAT_SPD, PETSC_TRUE);  // 对称正定
  
  // 步骤4：计算 f = A^T b
  VecCreate(PETSC_COMM_WORLD, &f_petsc);
  VecSetSizes(f_petsc, PETSC_DECIDE, num_fluid_particles);
  VecSetType(f_petsc, VECSEQ);
  VecSetFromOptions(f_petsc);
  
  MatMultTranspose(A_petsc, b_petsc, f_petsc);
  
  // 清理临时对象
  MatDestroy(&A_transpose);
  VecDestroy(&diagonal_penalty);
  
  return true;
}

bool PPEMatrixBuilder::WriteDebugInfoToVTK(
    const FluidParticle& fluid_particles,
    Mat A_petsc,
    Vec b_petsc,
    const std::string& filename) const {
  
  if (A_petsc == NULL || b_petsc == NULL) {
    std::cerr << "错误：矩阵或向量为NULL" << std::endl;
    return false;
  }
  
  int num_fluid_particles = fluid_particles.particle_num;
  if (num_fluid_particles == 0) {
    std::cerr << "错误：流体粒子数为0" << std::endl;
    return false;
  }
  
  // 创建输出目录（如果不存在）
  std::filesystem::path file_path(filename);
  std::filesystem::path dir_path = file_path.parent_path();
  if (!dir_path.empty() && !std::filesystem::exists(dir_path)) {
    try {
      std::filesystem::create_directories(dir_path);
    } catch (const std::filesystem::filesystem_error& e) {
      std::cerr << "错误：无法创建输出目录 " << dir_path << ": " << e.what() << std::endl;
      return false;
    }
  }
  
  // 检查文件是否已存在，如果存在则只追加数据，否则创建新文件
  FileOperator file_op;
  bool file_exists = std::filesystem::exists(filename);
  
  if (!file_exists) {
    // 文件不存在，创建基础VTK文件（包含粒子位置和速度）
  if (!file_op.writeVTKBase(filename, fluid_particles)) {
    std::cerr << "错误：无法创建VTK基础文件" << std::endl;
    return false;
  }
  }
  // 如果文件已存在，直接追加数据（假设基础结构已经存在）
  
  // 从PETSc矩阵中提取对角线元素
  std::vector<double> diagonal_elements(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
    PetscInt row = static_cast<PetscInt>(i);
    PetscInt col = static_cast<PetscInt>(i);
    PetscScalar value;
    PetscErrorCode ierr = MatGetValues(A_petsc, 1, &row, 1, &col, &value);
    if (ierr != 0) {
      std::cerr << "警告：无法获取矩阵对角线元素 " << i << std::endl;
      diagonal_elements[i] = 0.0;
    } else {
      diagonal_elements[i] = value;
    }
  }
  
  // 从PETSc向量中提取右边项（批量获取，更高效）
  std::vector<double> right_hand_side(num_fluid_particles);
  std::vector<PetscInt> indices(num_fluid_particles);
  std::vector<PetscScalar> values(num_fluid_particles);
  for (int i = 0; i < num_fluid_particles; ++i) {
    indices[i] = static_cast<PetscInt>(i);
  }
  PetscErrorCode ierr = VecGetValues(b_petsc, num_fluid_particles, indices.data(), values.data());
  if (ierr != 0) {
    std::cerr << "错误：无法从向量中获取值" << std::endl;
    return false;
  }
  for (int i = 0; i < num_fluid_particles; ++i) {
    right_hand_side[i] = values[i];
  }
  
  // 追加对角线元素到VTK文件
  if (!file_op.appendVTKScalar(filename, "diagonal_coefficients", diagonal_elements)) {
    std::cerr << "警告：无法追加对角线元素到VTK文件" << std::endl;
    return false;
  }
  
  // 追加右边项到VTK文件
  if (!file_op.appendVTKScalar(filename, "right_hand_side", right_hand_side)) {
    std::cerr << "警告：无法追加右边项到VTK文件" << std::endl;
    return false;
  }
  
  return true;
}

void PPEMatrixBuilder::InitializePetscMatrixAndVector(
    int num_particles,
    const std::vector<int>& nnz_per_row,
    Mat& A_petsc,
    Vec& b_petsc) const {
  
  PetscInt m = static_cast<PetscInt>(num_particles);
  PetscInt n = static_cast<PetscInt>(num_particles);
  
  // 转换为PetscInt向量
  std::vector<PetscInt> nnz_per_row_petsc(nnz_per_row.begin(), nnz_per_row.end());
  
  // 创建或重置PETSc矩阵
  if (A_petsc == NULL) {
    MatCreate(PETSC_COMM_WORLD, &A_petsc);
    MatSetSizes(A_petsc, PETSC_DECIDE, PETSC_DECIDE, m, n);
    MatSetType(A_petsc, MATSEQAIJ);
    MatSeqAIJSetPreallocation(A_petsc, 0, nnz_per_row_petsc.data());
    MatSetUp(A_petsc);
  } else {
    MatZeroEntries(A_petsc);
  }
  
  // 创建或重置PETSc向量
  if (b_petsc == NULL) {
    VecCreate(PETSC_COMM_WORLD, &b_petsc);
    VecSetSizes(b_petsc, PETSC_DECIDE, m);
    VecSetType(b_petsc, VECSEQ);
    VecSetFromOptions(b_petsc);
    VecSet(b_petsc, 0.0);
  } else {
    VecSet(b_petsc, 0.0);
  }
}

void PPEMatrixBuilder::BuildSurfaceParticleRow(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    double smoothing_radius,
    double density,
    double time_step,
    double reference_density,
    double lambda,
    Mat& A_petsc,
    Vec& b_petsc) const {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  const double2& vel_i = fluid_particles.velocity[particle_idx];
  PetscInt row = static_cast<PetscInt>(particle_idx);
  
  // 计算系数：4/(n0*lambda*rho)
  double surface_coeff_factor = 4.0 / (reference_density * lambda * density);
  
  // 初始化累加变量
  // sum_weight 累加所有邻域粒子的权重，实际上就是 n_i*（粒子数密度，考虑壁面粒子）
  double sum_weight = 0.0;  // sum(w_ij) for all j ≠ i (包括壁面粒子) = n_i*
  double divergence = 0.0;
  
  // 遍历流体邻域粒子 (j ∈ fluid, j ≠ i)
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
    sum_weight += weight;
    
    // 非对角线系数：4*w_ij/(n0*lambda*rho)
    PetscInt col = static_cast<PetscInt>(j);
    double off_diag_coeff = surface_coeff_factor * weight;
    MatSetValue(A_petsc, row, col, off_diag_coeff, ADD_VALUES);
    
    // 计算速度散度项：2/(n0*dt) * sum[(u_j - u_i)/r_ij · r_ij/r_ij * w_ij]
    // 根据文档：速度散度计算中 j ≠ i 包括所有邻域粒子
    double2 vel_diff = {vel_j.x - vel_i.x, vel_j.y - vel_i.y};
    double r_ij_mag = dist;
    double2 r_ij_unit = {dx / r_ij_mag, dy / r_ij_mag};
    double dot_product = vel_diff.x * r_ij_unit.x + vel_diff.y * r_ij_unit.y;
    divergence += (dot_product / r_ij_mag) * weight;
  }
  
  // 遍历壁面邻域粒子 (j ∈ wall, j ≠ i)
  // 注意：壁面粒子不参与压力矩阵的非对角线项，但参与速度散度计算和权重累加
  for (int j : fluid_particles.solid_neighbour_list[particle_idx]) {
    const double2& pos_j = solid_particles.position[j];
    const double2& vel_wall = solid_particles.velocity[j];
    
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    
    if (dist < 1e-10 || dist > smoothing_radius) {
      continue;
    }
    
    double weight = WeightFunction(dist, smoothing_radius);
    sum_weight += weight;  // 壁面粒子也参与权重累加
    
    // 计算速度散度项（壁面粒子参与速度散度计算）
    double2 vel_diff = {vel_wall.x - vel_i.x, vel_wall.y - vel_i.y};
    double r_ij_mag = dist;
    double2 r_ij_unit = {dx / r_ij_mag, dy / r_ij_mag};
    double dot_product = vel_diff.x * r_ij_unit.x + vel_diff.y * r_ij_unit.y;
    divergence += (dot_product / r_ij_mag) * weight;
  }
  
  // 循环结束后，sum_weight 就是 n_i*（粒子数密度，考虑壁面粒子）
  double n_i_star = sum_weight;
  
  // 根据文档：n_i' = min(n_i*, n_0)
  double n_i_prime = std::min(n_i_star, reference_density);
  
  // 根据文档第36行：对角线系数 = -4/(n0*lambda*rho) * (sum(w_ij) + (n_0 - n_i'))
  // 其中 sum(w_ij) 包括所有 j ≠ i 的权重（流体+壁面）
  double diag_coeff = -surface_coeff_factor * (sum_weight + (reference_density - n_i_prime));
  // double diag_coeff = -surface_coeff_factor * (sum_weight);
  
  MatSetValue(A_petsc, row, row, diag_coeff, ADD_VALUES);
  
  // 右边项：2/(n0*dt) * sum[...]
  // 根据文档：速度散度计算中 j ≠ i 包括所有邻域粒子（流体+壁面）
  double b_value = (2.0 / (reference_density * time_step)) * divergence;
  VecSetValue(b_petsc, row, b_value, INSERT_VALUES);
}

void PPEMatrixBuilder::BuildInnerParticleRow(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>& corrective_matrix_velocity,
    const Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>& corrective_matrix_pressure,
    CorrectiveMatrix& corrective_matrix_calc,
    double smoothing_radius,
    double density,
    double time_step,
    double gravity_x,
    double gravity_y,
    double coeff_factor,
    Mat& A_petsc,
    Vec& b_petsc) const {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  const double2& vel_i = fluid_particles.velocity[particle_idx];
  PetscInt row = static_cast<PetscInt>(particle_idx);
  
  // 提取速度corrective matrix的行向量（用于速度散度计算，第一类边界条件）
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C1_velocity = corrective_matrix_velocity.row(0);
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C2_velocity = corrective_matrix_velocity.row(1);
  
  // 提取压力corrective matrix的行向量（用于压力拉普拉斯算子计算，第二类边界条件）
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C3_pressure = corrective_matrix_pressure.row(2);
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C4_pressure = corrective_matrix_pressure.row(3);
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C3_plus_C4_pressure = C3_pressure + C4_pressure;
  
  // 初始化累加变量
  double diag_sum = 0.0;
  double divergence = 0.0;
  double wall_pressure_term = 0.0;
  
  // 遍历流体邻域粒子
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
    
    // 计算基函数
    Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis = 
        corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
    
    // 计算公共项（使用压力corrective matrix，第二类边界条件）
    double common_coeff = (weight / dist) * (C3_plus_C4_pressure * basis)(0, 0);
    
    // 累加对角线系数（负号）
    diag_sum += common_coeff;
    
    // 添加非对角线系数（正号）
    PetscInt col = static_cast<PetscInt>(j);
    double off_diag_coeff = coeff_factor * common_coeff;
    MatSetValue(A_petsc, row, col, off_diag_coeff, ADD_VALUES);
    
    // 计算速度散度项（使用速度corrective matrix，第一类边界条件）
    double dux_dr = (vel_j.x - vel_i.x) / dist;
    double duy_dr = (vel_j.y - vel_i.y) / dist;
    double C1P = (C1_velocity * basis)(0, 0);
    double C2P = (C2_velocity * basis)(0, 0);
    divergence += weight * (C1P * dux_dr + C2P * duy_dr);
  }
  
  // 处理壁面邻域粒子
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
    
    // 对于速度散度：使用第一类边界条件的基函数（标准基函数）
    Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis_velocity = 
        corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
    
    // 对于壁面压力边界条件项：使用第二类边界条件的基函数（壁面基函数）
    Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis_pressure = 
        corrective_matrix_calc.ComputeBasisFunctionsForWall(
            dx, dy, normal.x, normal.y, smoothing_radius);
    
    // 计算速度散度项（使用速度corrective matrix，第一类边界条件）
    double dux_dr = (vel_wall.x - vel_i.x) / dist;
    double duy_dr = (vel_wall.y - vel_i.y) / dist;
    double C1P = (C1_velocity * basis_velocity)(0, 0);
    double C2P = (C2_velocity * basis_velocity)(0, 0);
    divergence += weight * (C1P * dux_dr + C2P * duy_dr);
    
    // 计算壁面压力边界条件项（使用压力corrective matrix，第二类边界条件）
    double n_dot_g = normal.x * gravity_x + normal.y * gravity_y;
    double wall_pressure_coeff = -weight * (density * n_dot_g) * (C3_plus_C4_pressure * basis_pressure)(0, 0);
    wall_pressure_term += wall_pressure_coeff;
  }
  
  // 设置对角线系数
  double diag_coeff = -coeff_factor * diag_sum;
  MatSetValue(A_petsc, row, row, diag_coeff, ADD_VALUES);
  
  // 设置右边项
  double b_value = (1.0 / time_step) * divergence + coeff_factor * wall_pressure_term;
  VecSetValue(b_petsc, row, b_value, INSERT_VALUES);
}

} // namespace mps2D
