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

namespace {

// 与 BuildInnerParticleRow 中流体邻域项一致，仅返回拉普拉斯离散对角系数（不含壁面项）
double ComputeLaplacianDiagonalCoeffOnly(
    int particle_idx,
    const FluidParticle& fluid_particles,
    const Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>& corrective_matrix_pressure,
    CorrectiveMatrix& corrective_matrix_calc,
    double smoothing_radius,
    double coeff_factor) {
  const double2& pos_i = fluid_particles.position[particle_idx];
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> M2_pressure =
      corrective_matrix_pressure.row(2);
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> M3_pressure =
      corrective_matrix_pressure.row(3);
  Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> M2_plus_M3_pressure =
      M2_pressure + M3_pressure;

  double diag_sum = 0.0;
  for (int j : fluid_particles.fluid_neighbour_list[particle_idx]) {
    const double2& pos_j = fluid_particles.position[j];
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = ComputeDistance(pos_i, pos_j);
    if (dist < 1e-10 || dist > smoothing_radius) {
      continue;
    }
    double weight = WeightFunction(dist, smoothing_radius);
    Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis =
        corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
    double common_coeff = weight * (M2_plus_M3_pressure * basis)(0, 0);
    diag_sum += common_coeff;
  }
  return -coeff_factor * diag_sum;
}

// 取与内部粒子同阶的对角元尺度，用于飞溅粒子 PPE 行
double ChooseSplashReferenceDiagonal(
    const FluidParticle& fluid_particles,
    const std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>>& corrective_matrices_pressure,
    double smoothing_radius,
    double coeff_factor) {
  CorrectiveMatrix corrective_matrix_calc;
  auto diag_for = [&](int idx) -> double {
    return ComputeLaplacianDiagonalCoeffOnly(
        idx, fluid_particles, corrective_matrices_pressure[idx],
        corrective_matrix_calc, smoothing_radius, coeff_factor);
  };
  const int n = fluid_particles.particle_num;
  for (int i = 0; i < n; ++i) {
    if (static_cast<size_t>(i) >= fluid_particles.surface_type.size()) {
      break;
    }
    if (fluid_particles.surface_type[i] == SurfaceType::INNER) {
      double v = diag_for(i);
      if (std::abs(v) > 1e-30) {
        return v;
      }
    }
  }
  for (int i = 0; i < n; ++i) {
    if (static_cast<size_t>(i) >= fluid_particles.surface_type.size()) {
      break;
    }
    if (fluid_particles.surface_type[i] != SurfaceType::SPLASH) {
      double v = diag_for(i);
      if (std::abs(v) > 1e-30) {
        return v;
      }
    }
  }
  return -std::max(std::abs(coeff_factor), 1e-12);
}

// 罚函数正规方程 Kp=f 上强制 p_i=0（齐次 Dirichlet），保证飞溅自由度严格为零
void ApplySplashDirichletOnPenaltySystem(Mat K_petsc, Vec f_petsc,
                                         const FluidParticle& fluid_particles) {
  const int n = fluid_particles.particle_num;
  std::vector<PetscInt> rows;
  for (int i = 0; i < n; ++i) {
    if (static_cast<size_t>(i) < fluid_particles.surface_type.size() &&
        fluid_particles.surface_type[i] == SurfaceType::SPLASH) {
      rows.push_back(static_cast<PetscInt>(i));
    }
  }
  if (rows.empty()) {
    return;
  }
  PetscScalar diag_bc = 1.0;
  for (int i = 0; i < n; ++i) {
    if (static_cast<size_t>(i) < fluid_particles.surface_type.size() &&
        fluid_particles.surface_type[i] != SurfaceType::SPLASH) {
      PetscInt r = static_cast<PetscInt>(i);
      PetscScalar v = 0.0;
      PetscErrorCode ierr_g = MatGetValue(K_petsc, r, r, &v);
      if (ierr_g != 0) {
        continue;
      }
      double vd = std::abs(static_cast<double>(v));
      if (vd > 1e-30) {
        diag_bc = static_cast<PetscScalar>(vd);
        break;
      }
    }
  }
  PetscErrorCode ierr = MatZeroRowsColumns(
      K_petsc, static_cast<PetscInt>(rows.size()), rows.data(), diag_bc,
      nullptr, f_petsc);
  if (ierr != 0) {
    std::cerr << "警告：MatZeroRowsColumns（飞溅粒子 Dirichlet）失败， ierr=" << ierr
              << std::endl;
    return;
  }
  MatAssemblyBegin(K_petsc, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K_petsc, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(f_petsc);
  VecAssemblyEnd(f_petsc);
}

}  // namespace

// 确保PETSc已在程序入口处初始化（用于BuildPPEMatrixPetsc）
static bool EnsurePetscInitialized() {
  PetscBool initialized = PETSC_FALSE;
  PetscInitialized(&initialized);
  if (!initialized) {
    std::cerr << "错误：PETSc尚未初始化，请在程序入口调用PetscInitialize。" << std::endl;
    return false;
  }
  return true;
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
    Vec& b_petsc,
    std::vector<double>* velocity_divergence_out,
    std::vector<double>* diagonal_abs_out,
    std::vector<double>* rhs_abs_out) {
  
  // 确保PETSc已在程序入口处初始化
  if (!EnsurePetscInitialized()) {
    return false;
  }
  
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
  
  // 如果需要调试输出速度散度，则初始化输出向量
  if (velocity_divergence_out != nullptr) {
    velocity_divergence_out->assign(num_fluid_particles, 0.0);
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

  const double splash_ref_diag = ChooseSplashReferenceDiagonal(
      fluid_particles, corrective_matrices_pressure, smoothing_radius,
      coeff_factor);
  
  // 遍历所有流体粒子，构建系数矩阵和右边项
  for (int particle_idx = 0; particle_idx < num_fluid_particles; ++particle_idx) {
    // 飞溅粒子：仅对角元（尺度与内部离散同阶）、非对角为零、右端项为零，使 p_i = 0
    if (static_cast<size_t>(particle_idx) < fluid_particles.surface_type.size() &&
        fluid_particles.surface_type[particle_idx] == SurfaceType::SPLASH) {
      PetscInt row = static_cast<PetscInt>(particle_idx);
      MatSetValue(A_petsc, row, row, splash_ref_diag, ADD_VALUES);
      VecSetValue(b_petsc, row, 0.0, INSERT_VALUES);
      if (velocity_divergence_out != nullptr) {
        (*velocity_divergence_out)[particle_idx] = 0.0;
      }
      continue;
    }
    // 速度散度使用第一类边界条件的corrective matrix
    // 压力拉普拉斯算子使用第二类边界条件的corrective matrix
    double divergence_value = 0.0;
    BuildInnerParticleRow(
        particle_idx, fluid_particles, solid_particles,
        corrective_matrices_velocity[particle_idx],
        corrective_matrices_pressure[particle_idx],
        corrective_matrix_calc, smoothing_radius, density, time_step,
        gravity_x, gravity_y, coeff_factor, A_petsc, b_petsc,
        (velocity_divergence_out != nullptr) ? &divergence_value : nullptr);
    if (velocity_divergence_out != nullptr) {
      (*velocity_divergence_out)[particle_idx] = divergence_value;
    }
  }
  
  // 组装矩阵和向量
  MatAssemblyBegin(A_petsc, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A_petsc, MAT_FINAL_ASSEMBLY);
  
  VecAssemblyBegin(b_petsc);
  VecAssemblyEnd(b_petsc);
  
  // 可选调试输出：提取 PPE 系数矩阵 A 的主对角线绝对值，以及源项向量 b 的绝对值
  if (diagonal_abs_out != nullptr || rhs_abs_out != nullptr) {
    if (diagonal_abs_out != nullptr) {
      diagonal_abs_out->assign(num_fluid_particles, 0.0);
    }
    if (rhs_abs_out != nullptr) {
      rhs_abs_out->assign(num_fluid_particles, 0.0);
    }

    // 提取对角线元素 |A_ii|
    if (diagonal_abs_out != nullptr) {
      for (int i = 0; i < num_fluid_particles; ++i) {
        PetscInt row = static_cast<PetscInt>(i);
        PetscInt col = static_cast<PetscInt>(i);
        PetscScalar value = 0.0;
        PetscErrorCode ierr = MatGetValues(A_petsc, 1, &row, 1, &col, &value);
        if (ierr == 0) {
          (*diagonal_abs_out)[i] = std::abs(static_cast<double>(value));
        } else {
          (*diagonal_abs_out)[i] = 0.0;
        }
      }
    }

    // 提取源项向量元素 |b_i|
    if (rhs_abs_out != nullptr) {
      std::vector<PetscInt> indices(num_fluid_particles);
      std::vector<PetscScalar> values(num_fluid_particles);
      for (int i = 0; i < num_fluid_particles; ++i) {
        indices[i] = static_cast<PetscInt>(i);
      }
      PetscErrorCode ierr =
          VecGetValues(b_petsc, num_fluid_particles, indices.data(), values.data());
      if (ierr == 0) {
        for (int i = 0; i < num_fluid_particles; ++i) {
          (*rhs_abs_out)[i] = std::abs(static_cast<double>(values[i]));
        }
      } else {
        // 失败则保持默认 0.0
      }
    }
  }

  return true;
}

bool PPEMatrixBuilder::BuildPPEPenaltyNormalEquationPetsc(
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
    double penalty_mu,
    Mat& K_petsc,
    Vec& f_petsc,
    std::vector<double>* velocity_divergence_out,
    std::vector<double>* diagonal_abs_out,
    std::vector<double>* rhs_abs_out) {
  // 1) 先构建原始系统 A p = b（不对自由面行做行修改）
  Mat A_petsc = NULL;
  Vec b_petsc = NULL;
  if (!BuildPPEMatrixPetsc(
          fluid_particles, solid_particles,
          corrective_matrices_velocity, corrective_matrices_pressure,
          smoothing_radius, density, time_step, particle_spacing,
          gravity_x, gravity_y,
          A_petsc, b_petsc,
          velocity_divergence_out, diagonal_abs_out, rhs_abs_out)) {
    if (A_petsc != NULL) MatDestroy(&A_petsc);
    if (b_petsc != NULL) VecDestroy(&b_petsc);
    return false;
  }

  // 2) 计算 K = A^T A
  if (K_petsc != NULL) {
    MatDestroy(&K_petsc);
    K_petsc = NULL;
  }
  PetscErrorCode ierr = MatTransposeMatMult(
      A_petsc, A_petsc, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &K_petsc);
  if (ierr != 0 || K_petsc == NULL) {
    std::cerr << "错误：构建 K = A^T A 失败" << std::endl;
    if (A_petsc != NULL) MatDestroy(&A_petsc);
    if (b_petsc != NULL) VecDestroy(&b_petsc);
    if (K_petsc != NULL) MatDestroy(&K_petsc);
    K_petsc = NULL;
    return false;
  }

  // 3) 计算 f = A^T b
  if (f_petsc == NULL) {
    ierr = VecDuplicate(b_petsc, &f_petsc);
    if (ierr != 0 || f_petsc == NULL) {
      std::cerr << "错误：创建向量 f 失败" << std::endl;
      MatDestroy(&A_petsc);
      VecDestroy(&b_petsc);
      MatDestroy(&K_petsc);
      K_petsc = NULL;
      return false;
    }
  } else {
    VecSet(f_petsc, 0.0);
  }
  ierr = MatMultTranspose(A_petsc, b_petsc, f_petsc);
  if (ierr != 0) {
    std::cerr << "错误：构建 f = A^T b 失败" << std::endl;
    MatDestroy(&A_petsc);
    VecDestroy(&b_petsc);
    MatDestroy(&K_petsc);
    K_petsc = NULL;
    VecDestroy(&f_petsc);
    f_petsc = NULL;
    return false;
  }

  // 4) 加入罚项 D（自由面粒子对角线加 penalty_mu）
  if (penalty_mu > 0.0) {
    const int n = fluid_particles.particle_num;
    for (int i = 0; i < n; ++i) {
      if (fluid_particles.surface_type[i] == SurfaceType::SURFACE) {
        PetscInt row = static_cast<PetscInt>(i);
        MatSetValue(K_petsc, row, row, penalty_mu, ADD_VALUES);
      }
    }
  }

  MatAssemblyBegin(K_petsc, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K_petsc, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(f_petsc);
  VecAssemblyEnd(f_petsc);

  // 飞溅粒子：在正规方程上施加 p=0（A 中已隔离行不足以单独保证 K=A^T A 下 p_i=0）
  ApplySplashDirichletOnPenaltySystem(K_petsc, f_petsc, fluid_particles);

  // 5) 清理临时 A、b
  MatDestroy(&A_petsc);
  VecDestroy(&b_petsc);
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
    Vec& b_petsc,
    double* velocity_divergence_out) const {
  
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
    // 壁面贡献使用有效速度 vel_wall + Δt*g，不修改原始 solid 数组
    double vel_wall_eff_x = vel_wall.x + time_step * gravity_x;
    double vel_wall_eff_y = vel_wall.y + time_step * gravity_y;
    double dvx = vel_wall_eff_x - vel_i.x;
    double dvy = vel_wall_eff_y - vel_i.y;
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
  double divergence_term = (1.0 / (smoothing_radius * time_step)) * divergence_sum;
  double b_value = divergence_term + coeff_factor * wall_pressure_term;
  VecSetValue(b_petsc, row, b_value, INSERT_VALUES);
  
  // 如果需要调试输出，则返回当前粒子的临时速度散度（仅速度项）
  if (velocity_divergence_out != nullptr) {
    *velocity_divergence_out = divergence_term;
  }
}

} // namespace mps2D
