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
    const std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>>& corrective_matrices,
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
  
  if (static_cast<int>(corrective_matrices.size()) != num_fluid_particles) {
    std::cerr << "错误：corrective matrix数量与流体粒子数不匹配" << std::endl;
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
  
  // 计算参考粒子数密度n0和lambda（用于自由面粒子）
  OriginalMPS mps_calculator;
  double reference_density = mps_calculator.ComputeReferenceDensity(
      fluid_particles, solid_particles, smoothing_radius);
  // lambda计算公式：λ = (1/5) * r_e^2
  double lambda = (1.0 / 5.0) * smoothing_radius * smoothing_radius;
  
  // 计算w_l0（在l0距离处的权重函数值，l0为粒子初始间距）
  double w_l0 = WeightFunction(particle_spacing, smoothing_radius);
  
  // 系数矩阵的系数因子（用于内部粒子）
  double coeff_factor = 2.0 / (smoothing_radius * density);
  
  // 遍历所有流体粒子，构建系数矩阵和右边项
  for (int particle_idx = 0; particle_idx < num_fluid_particles; ++particle_idx) {
    bool is_surface_particle = (fluid_particles.surface_type[particle_idx] == SurfaceType::SURFACE);
    
    if (is_surface_particle) {
      // 自由面粒子：使用文档中的离散方法
      BuildSurfaceParticleRow(
          particle_idx, fluid_particles, smoothing_radius, density, time_step,
          reference_density, lambda, w_l0, A_petsc, b_petsc);
    } else {
      // 内部粒子：使用原有的LSMPS方法
      BuildInnerParticleRow(
          particle_idx, fluid_particles, solid_particles, corrective_matrices[particle_idx],
          corrective_matrix_calc, smoothing_radius, density, time_step,
          gravity_x, gravity_y, coeff_factor, A_petsc, b_petsc);
    }
  }
  
  // 组装矩阵和向量
  MatAssemblyBegin(A_petsc, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A_petsc, MAT_FINAL_ASSEMBLY);
  
  VecAssemblyBegin(b_petsc);
  VecAssemblyEnd(b_petsc);
  
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
  
  // 使用FileOperator创建基础VTK文件（包含粒子位置和速度）
  FileOperator file_op;
  if (!file_op.writeVTKBase(filename, fluid_particles)) {
    std::cerr << "错误：无法创建VTK基础文件" << std::endl;
    return false;
  }
  
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
    double smoothing_radius,
    double density,
    double time_step,
    double reference_density,
    double lambda,
    double w_l0,
    Mat& A_petsc,
    Vec& b_petsc) const {
  
  const double2& pos_i = fluid_particles.position[particle_idx];
  const double2& vel_i = fluid_particles.velocity[particle_idx];
  PetscInt row = static_cast<PetscInt>(particle_idx);
  
  // 计算n_i*（粒子数密度，只计算流体粒子）
  double n_i_star = 0.0;
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    double dist = ComputeDistance(pos_i, fluid_particles.position[j]);
    if (dist > 1e-10 && dist <= smoothing_radius) {
      n_i_star += WeightFunction(dist, smoothing_radius);
    }
  }
  
  // 计算tilde_n_i = n_0 + 0.1 * sum(w_ij - w_l0)
  double sum_w_diff = 0.0;
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    double dist = ComputeDistance(pos_i, fluid_particles.position[j]);
    if (dist > 1e-10 && dist <= smoothing_radius) {
      double w_ij = WeightFunction(dist, smoothing_radius);
      sum_w_diff += (w_ij - w_l0);
    }
  }
  double tilde_n_i = reference_density + 0.1 * sum_w_diff;
  
  // n_i' = max(n_i*, tilde_n_i)
  double n_i_prime = std::max(n_i_star, tilde_n_i);
  
  // 计算系数：4/(n0*lambda*rho)
  double surface_coeff_factor = 4.0 / (reference_density * lambda * density);
  
  // 初始化累加变量
  double sum_weight = 0.0;
  double divergence = 0.0;
  
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
    sum_weight += weight;
    
    // 非对角线系数：4*w_ij/(n0*lambda*rho)
    PetscInt col = static_cast<PetscInt>(j);
    double off_diag_coeff = surface_coeff_factor * weight;
    MatSetValue(A_petsc, row, col, off_diag_coeff, ADD_VALUES);
    
    // 计算速度散度项：2/(n0*dt) * sum[(u_j - u_i)/r_ij · r_ij/r_ij * w_ij]
    double2 vel_diff = {vel_j.x - vel_i.x, vel_j.y - vel_i.y};
    double r_ij_mag = dist;
    double2 r_ij_unit = {dx / r_ij_mag, dy / r_ij_mag};
    double dot_product = vel_diff.x * r_ij_unit.x + vel_diff.y * r_ij_unit.y;
    divergence += (dot_product / r_ij_mag) * weight;
  }
  
  // 对角线系数：-4/(n0*lambda*rho) * (sum(w_ij) - (n_i' - n_i*))
  double diag_coeff = -surface_coeff_factor * (sum_weight + (n_i_prime - n_i_star));
  
  // 调试输出：针对2775和2997号粒子
  if (particle_idx == 2775 || particle_idx == 2997) {
    std::cout << "\n=== 粒子 " << particle_idx << " 的调试信息 ===" << std::endl;
    std::cout << "  n_i* (实际粒子数密度): " << n_i_star << std::endl;
    std::cout << "  tilde_n_i: " << tilde_n_i << std::endl;
    std::cout << "  n_i' (max(n_i*, tilde_n_i)): " << n_i_prime << std::endl;
    std::cout << "  n_i' - n_i*: " << (n_i_prime - n_i_star) << std::endl;
    std::cout << "  sum_weight (sum(w_ij)): " << sum_weight << std::endl;
    std::cout << "  sum_weight - (n_i' - n_i*): " << (sum_weight - (n_i_prime - n_i_star)) << std::endl;
    std::cout << "  surface_coeff_factor (4/(n0*lambda*rho)): " << surface_coeff_factor << std::endl;
    std::cout << "  对角线系数 (diag_coeff): " << diag_coeff << std::endl;
    std::cout << "  邻域粒子数量: " << fluid_particles.fluid_neighbour_list[particle_idx].size() << std::endl;
    std::cout << "  reference_density (n_0): " << reference_density << std::endl;
    std::cout << "  lambda: " << lambda << std::endl;
    std::cout << "  w_l0: " << w_l0 << std::endl;
    std::cout << "  sum_w_diff: " << sum_w_diff << std::endl;
    std::cout << "========================================\n" << std::endl;
  }
  
  MatSetValue(A_petsc, row, row, diag_coeff, ADD_VALUES);
  
  // 右边项：2/(n0*dt) * sum[...]
  double b_value = (2.0 / (reference_density * time_step)) * divergence;
  VecSetValue(b_petsc, row, b_value, INSERT_VALUES);
}

void PPEMatrixBuilder::BuildInnerParticleRow(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles,
    const Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>& corrective_matrix,
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
  
  // 提取corrective matrix的行向量
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C1 = corrective_matrix.row(0);
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C2 = corrective_matrix.row(1);
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C3 = corrective_matrix.row(2);
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C4 = corrective_matrix.row(3);
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C3_plus_C4 = C3 + C4;
  
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
    
    // 计算公共项
    double common_coeff = (weight / dist) * (C3_plus_C4 * basis)(0, 0);
    
    // 累加对角线系数（负号）
    diag_sum += common_coeff;
    
    // 添加非对角线系数（正号）
    PetscInt col = static_cast<PetscInt>(j);
    double off_diag_coeff = coeff_factor * common_coeff;
    MatSetValue(A_petsc, row, col, off_diag_coeff, ADD_VALUES);
    
    // 计算速度散度项
    double dux_dr = (vel_j.x - vel_i.x) / dist;
    double duy_dr = (vel_j.y - vel_i.y) / dist;
    double C1P = (C1 * basis)(0, 0);
    double C2P = (C2 * basis)(0, 0);
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
    
    // 对于速度散度：使用第一类边界条件的基函数
    Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis_velocity = 
        corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
    
    // 对于壁面压力边界条件项：使用第二类边界条件的基函数
    Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis_pressure = 
        corrective_matrix_calc.ComputeBasisFunctionsForWall(
            dx, dy, normal.x, normal.y, smoothing_radius);
    
    // 计算速度散度项
    double dux_dr = (vel_wall.x - vel_i.x) / dist;
    double duy_dr = (vel_wall.y - vel_i.y) / dist;
    double C1P = (C1 * basis_velocity)(0, 0);
    double C2P = (C2 * basis_velocity)(0, 0);
    divergence += weight * (C1P * dux_dr + C2P * duy_dr);
    
    // 计算壁面压力边界条件项
    double n_dot_g = normal.x * gravity_x + normal.y * gravity_y;
    double wall_pressure_coeff = -weight * (density * n_dot_g) * (C3_plus_C4 * basis_pressure)(0, 0);
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
