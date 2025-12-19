#include "PPESolver.hpp"
#include <iostream>
#include <cstring>  // for strcmp

// PETSc头文件
#include <petsc.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>

namespace mps2D {

// 静态标志，确保PETSc只初始化一次
static bool petsc_initialized = false;

PPESolver::PPESolver(const SolverConfig& config) : config_(config) {
  // 初始化PETSc（只初始化一次）
  if (!petsc_initialized) {
    int argc = 0;
    char** argv = nullptr;
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    petsc_initialized = true;
    
    // 设置PETSc选项：不显示版权信息
    PetscOptionsSetValue(nullptr, "-options_left", "false");
  }
}

PPESolver::~PPESolver() {
  // PETSc的清理在程序结束时统一处理，不需要在每个对象析构时处理
  // 如果需要，可以添加PetscFinalize()的调用，但通常不建议
}

bool PPESolver::Solve(
    Mat A_petsc,
    Vec b_petsc,
    Vec& p_petsc) {
  
  if (A_petsc == NULL || b_petsc == NULL) {
    std::cerr << "错误：PETSc矩阵或向量为NULL" << std::endl;
    return false;
  }
  
  // 创建解向量（如果尚未创建）
  if (p_petsc == NULL) {
    VecDuplicate(b_petsc, &p_petsc);
    VecSet(p_petsc, 0.0);  // 初始化为零向量
  } else {
    VecSet(p_petsc, 0.0);  // 重置为零向量
  }
  
  // 创建KSP求解器
  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  
  // 设置矩阵
  KSPSetOperators(ksp, A_petsc, A_petsc);
  
  // 根据配置选择求解器类型
  if (config_.solver_type == SolverType::BICGSTAB) {
    KSPSetType(ksp, KSPBCGS);
  } else if (config_.solver_type == SolverType::GMRES) {
    KSPSetType(ksp, KSPGMRES);
    KSPGMRESSetRestart(ksp, config_.restart);
  } else {
    std::cerr << "错误：未知的求解器类型" << std::endl;
    KSPDestroy(&ksp);
    return false;
  }
  
  // 设置预处理类型
  PC pc;
  KSPGetPC(ksp, &pc);
  
  // 检查矩阵规模，决定使用直接求解器还是迭代求解器
  PetscInt m_global, n_global;
  MatGetSize(A_petsc, &m_global, &n_global);
  
  // 对于中等规模矩阵（<10000），使用直接求解器（LU分解）可以达到机器精度
  // 直接求解器的精度：理论上可以达到机器精度（双精度约1e-15到1e-16）
  if (m_global < 10000) {
    // 使用直接求解器（LU分解）
    // 注意：对于直接求解器，KSP实际上只需要1次迭代
    KSPSetType(ksp, KSPPREONLY);  // 只应用预处理器，不迭代
    PCSetType(pc, PCLU);  // 使用LU分解作为"预处理器"（实际上是直接求解）
    
    // 尝试使用PETSc内置的LU求解器
    // 如果可用，使用更高效的求解器（如MUMPS, SuperLU等）
    PCFactorSetMatSolverType(pc, MATSOLVERPETSC);
    
    std::cout << "    使用直接求解器（LU分解），目标精度：机器精度（~1e-15）" << std::endl;
  } else {
    // 对于大矩阵，使用增强的ILU预处理
    PCSetType(pc, PCILU);
    PCFactorSetLevels(pc, 5);  // 使用5级fill-in
    PCFactorSetDropTolerance(pc, 1e-14, PETSC_DEFAULT, PETSC_DEFAULT);
    PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO);
    PCFactorSetShiftAmount(pc, 1e-10);
    PCFactorSetMatOrderingType(pc, MATORDERINGND);
    std::cout << "    使用迭代求解器（ILU预处理）" << std::endl;
  }
  
  // 计算右边项的范数，用于相对残差计算
  PetscReal b_norm;
  VecNorm(b_petsc, NORM_2, &b_norm);
  
  // 设置求解器参数
  // 使用相对残差：rtol = tolerance（相对容差），atol = tolerance * b_norm（绝对容差）
  // 这样PETSc会使用相对残差进行收敛判断
  KSPSetTolerances(ksp, config_.tolerance, config_.tolerance * b_norm, 
                   PETSC_DEFAULT, config_.max_iterations);
  
  // 设置使用相对残差进行收敛判断
  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
  
  // 设置从命令行选项读取参数（可选）
  KSPSetFromOptions(ksp);
  
  // 求解
  KSPSolve(ksp, b_petsc, p_petsc);
  
  // 获取求解信息
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);
  PetscInt iter_num;
  KSPGetIterationNumber(ksp, &iter_num);
  last_iterations_ = static_cast<int>(iter_num);
  
  // 检查是否使用了直接求解器
  KSPType ksp_type;
  KSPGetType(ksp, &ksp_type);
  bool is_direct_solver = (strcmp(ksp_type, KSPPREONLY) == 0);
  
  // 获取KSP计算的残差（这是PETSc内部使用的残差，可能已经预处理）
  PetscReal ksp_residual_norm = 0.0;
  KSPGetResidualNorm(ksp, &ksp_residual_norm);
  
  // 计算最终残差（使用相对残差，这是真实的未预处理残差）
  Vec residual;
  VecDuplicate(b_petsc, &residual);
  MatMult(A_petsc, p_petsc, residual);
  VecAXPY(residual, -1.0, b_petsc);
  PetscReal residual_norm;
  VecNorm(residual, NORM_2, &residual_norm);
  
  // 计算相对残差
  PetscReal b_norm_final;
  VecNorm(b_petsc, NORM_2, &b_norm_final);
  if (b_norm_final > 1e-15) {
    last_residual_ = residual_norm / b_norm_final;  // 相对残差
  } else {
    last_residual_ = residual_norm;  // 如果b的范数太小，使用绝对残差
  }
  
  VecDestroy(&residual);
  
  // 判断收敛性
  bool ksp_converged = (reason > 0);
  
  // 对于直接求解器，理论上应该达到机器精度
  // 机器精度约为1e-15（双精度），但由于矩阵条件数的影响，实际精度可能略低
  // 对于条件数κ的矩阵，直接求解器的精度约为：机器精度 × κ
  const double machine_epsilon = 1e-15;
  const double direct_solver_tolerance = machine_epsilon * 1e6;  // 允许条件数影响（约1e-9）
  
  if (is_direct_solver) {
    // 直接求解器：如果残差接近机器精度级别，认为成功
    // 注意：即使使用直接求解器，如果矩阵条件数很大，残差也可能较大
    if (last_residual_ < direct_solver_tolerance) {
      last_converged_ = true;
      std::cout << "    直接求解器达到高精度（残差: " << last_residual_ 
                << "，接近机器精度级别）" << std::endl;
    } else {
      last_converged_ = false;
      std::cout << "    警告：直接求解器残差较大（" << last_residual_ 
                << "），可能由于矩阵条件数较大" << std::endl;
      std::cout << "    理论机器精度: ~" << machine_epsilon << std::endl;
      std::cout << "    实际残差: " << last_residual_ << std::endl;
      std::cout << "    估计矩阵条件数: ~" << (last_residual_ / machine_epsilon) << std::endl;
    }
  } else {
    // 迭代求解器：检查是否满足容忍度
    bool residual_converged = (last_residual_ <= config_.tolerance);
    last_converged_ = ksp_converged && residual_converged;
    
    // 输出残差信息
    if (last_converged_) {
      std::cout << "    求解成功！" << std::endl;
    } else {
      if (!ksp_converged) {
        std::cout << "    求解未收敛（KSP收敛原因: " << reason << "，负数表示未收敛）" << std::endl;
      }
      if (!residual_converged) {
        std::cout << "    相对残差未达到目标（当前: " << last_residual_ 
                  << "，目标: " << config_.tolerance << "）" << std::endl;
      }
    }
  }
  
  // 对于直接求解器，输出额外信息
  if (is_direct_solver && !last_converged_) {
    // 已经在上面输出了
  } else if (is_direct_solver && last_converged_) {
    // 已经在上面输出了
  } else if (!is_direct_solver) {
    // 迭代求解器的输出已经在上面处理了
  }
  std::cout << "    迭代次数: " << last_iterations_ << std::endl;
  std::cout << "    绝对残差: " << residual_norm << std::endl;
  std::cout << "    相对残差: " << last_residual_ << " (目标: " << config_.tolerance << ")" << std::endl;
  if (ksp_residual_norm > 0) {
    std::cout << "    KSP内部残差: " << ksp_residual_norm << std::endl;
  }
  
  // 清理KSP对象
  KSPDestroy(&ksp);
  
  if (!last_converged_) {
    std::cerr << "警告：求解未收敛，迭代次数：" << last_iterations_ 
              << "，残差：" << last_residual_ << std::endl;
  }
  
  return last_converged_;
}

} // namespace mps2D
