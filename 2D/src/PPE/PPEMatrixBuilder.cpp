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
  // 根据文档：拉普拉斯算子前的系数为 2 / (r_s^2)
  // 这里同时除以密度 density
  double coeff_factor = 2.0 / (smoothing_radius * smoothing_radius * density);
  
  // 遍历所有流体粒子，构建系数矩阵和右边项
  for (int particle_idx = 0; particle_idx < num_fluid_particles; ++particle_idx) {
    // 自由面粒子：使用行修改法，直接施加 p_i = 0
    if (fluid_particles.surface_type[particle_idx] == SurfaceType::SURFACE) {
      BuildSurfaceParticleRowAdjusted(
          particle_idx, particle_spacing, density, A_petsc, b_petsc);
      continue;
    }

    // 非自由面粒子：使用LSMPS离散方法
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

void PPEMatrixBuilder::BuildSurfaceParticleRowAdjusted(
    int particle_idx,
    double particle_spacing,
    double density,
    Mat& A_petsc,
    Vec& b_petsc) const {
  PetscInt row = static_cast<PetscInt>(particle_idx);

  // 参考 `PPEadjust.md`：自由面粒子施加 p_i = 0
  // 取 c 接近内部粒子对角线尺度：c = 1 / (Δx^2 ρ)
  const double dx = (particle_spacing > 0.0) ? particle_spacing : 1e-12;
  const double c = 1.0 / (dx * dx * density);

  MatSetValue(A_petsc, row, row, c, INSERT_VALUES);
  VecSetValue(b_petsc, row, 0.0, INSERT_VALUES);
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
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> M2_pressure = corrective_matrix_pressure.row(2);
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> M3_pressure = corrective_matrix_pressure.row(3);
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> M2_plus_M3_pressure = M2_pressure + M3_pressure;
  
  // 初始化累加变量
  double diag_sum = 0.0;
  double divergence_sum = 0.0;
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
    
    // 计算基函数 P_ij
    Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis = 
        corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
    
    // 压力拉普拉斯公共项（仅使用基函数和 moment 矩阵逆）
    // 对应文档中的 [M_{i,2} + M_{i,3}] P_ij
    double common_coeff = weight * (M2_plus_M3_pressure * basis)(0, 0);
    
    // 累加对角线系数（负号）
    diag_sum += common_coeff;
    
    // 添加非对角线系数（正号）
    PetscInt col = static_cast<PetscInt>(j);
    double off_diag_coeff = coeff_factor * common_coeff;
    MatSetValue(A_petsc, row, col, off_diag_coeff, ADD_VALUES);
    
    // 计算速度散度项（使用速度corrective matrix，第一类边界条件）
    double dvx = vel_j.x - vel_i.x;
    double dvy = vel_j.y - vel_i.y;
    double C1P = (C1_velocity * basis)(0, 0);
    double C2P = (C2_velocity * basis)(0, 0);
    // 根据文档：∇·u 使用 1/r_s 系数，而差分项为 (u_j - u_i)
    divergence_sum += weight * (C1P * dvx + C2P * dvy);
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
    double dvx = vel_wall.x - vel_i.x;
    double dvy = vel_wall.y - vel_i.y;
    double C1P = (C1_velocity * basis_velocity)(0, 0);
    double C2P = (C2_velocity * basis_velocity)(0, 0);
    divergence_sum += weight * (C1P * dvx + C2P * dvy);
    
    // 计算壁面压力边界条件项（使用压力corrective matrix，第二类边界条件）
    double n_dot_g = normal.x * gravity_x + normal.y * gravity_y;
    // 根据文档：壁面项中包含 r_s ρ g n
    double wall_pressure_coeff = -weight * smoothing_radius * density * n_dot_g *
        (M2_plus_M3_pressure * basis_pressure)(0, 0);
    wall_pressure_term += wall_pressure_coeff;
  }
  
  // 设置对角线系数
  double diag_coeff = -coeff_factor * diag_sum;
  MatSetValue(A_petsc, row, row, diag_coeff, ADD_VALUES);
  
  // 设置右边项
  // 速度散度项前的系数 1 / (r_s * Δt)
  double b_value = (1.0 / (smoothing_radius * time_step)) * divergence_sum +
      coeff_factor * wall_pressure_term;
  VecSetValue(b_petsc, row, b_value, INSERT_VALUES);
}

} // namespace mps2D
