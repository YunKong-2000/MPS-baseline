#include "PPESolver.hpp"
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>

namespace mps2D {

PPESolver::PPESolver(const SolverConfig& config) : config_(config) {
}

bool PPESolver::Solve(
    const Eigen::SparseMatrix<double, Eigen::ColMajor>& A,
    const Eigen::VectorXd& b,
    Eigen::VectorXd& p) {
  
  if (A.rows() != A.cols() || A.rows() != b.size()) {
    std::cerr << "错误：矩阵和向量维度不匹配" << std::endl;
    return false;
  }
  
  // 确保矩阵是压缩格式
  if (!A.isCompressed()) {
    std::cerr << "警告：矩阵未压缩，正在压缩..." << std::endl;
    const_cast<Eigen::SparseMatrix<double, Eigen::ColMajor>&>(A).makeCompressed();
  }
  
  // 初始化解向量为零向量
  p = Eigen::VectorXd::Zero(b.size());
  
  // 根据配置选择求解器
  bool success = false;
  if (config_.solver_type == SolverType::BICGSTAB) {
    success = SolveBiCGSTAB(A, b, p);
  } else if (config_.solver_type == SolverType::GMRES) {
    success = SolveGMRES(A, b, p);
  } else {
    std::cerr << "错误：未知的求解器类型" << std::endl;
    return false;
  }
  
  return success;
}

bool PPESolver::SolveBiCGSTAB(
    const Eigen::SparseMatrix<double, Eigen::ColMajor>& A,
    const Eigen::VectorXd& b,
    Eigen::VectorXd& p) {
  
  // 使用Eigen的BiCGSTAB求解器
  // 使用压缩格式的稀疏矩阵（CSC格式）可以提高矩阵-向量乘法性能
  // 预处理器在模板参数中指定，通过compute方法设置
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::DiagonalPreconditioner<double>> solver;
  
  // 设置求解器参数
  solver.setMaxIterations(config_.max_iterations);
  solver.setTolerance(config_.tolerance);
  
  // 计算对角预处理器（通过compute方法设置）
  solver.compute(A);
  
  // 求解
  p = solver.solveWithGuess(b, p);
  
  // 记录求解信息
  last_iterations_ = solver.iterations();
  last_residual_ = solver.error();
  last_converged_ = (solver.info() == Eigen::Success);
  
  if (!last_converged_) {
    std::cerr << "警告：BiCGSTAB求解未收敛，迭代次数：" << last_iterations_ 
              << "，残差：" << last_residual_ << std::endl;
  }
  
  return last_converged_;
}

bool PPESolver::SolveGMRES(
    const Eigen::SparseMatrix<double, Eigen::ColMajor>& A,
    const Eigen::VectorXd& b,
    Eigen::VectorXd& p) {
  
  // 注意：Eigen库可能不包含GMRES求解器
  // 如果Eigen版本不支持GMRES，可以使用BiCGSTAB作为替代
  // 这里暂时使用BiCGSTAB作为GMRES的替代实现
  std::cerr << "警告：Eigen库可能不支持GMRES求解器，使用BiCGSTAB替代" << std::endl;
  
  // 使用BiCGSTAB作为替代
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::DiagonalPreconditioner<double>> solver;
  
  // 设置求解器参数
  solver.setMaxIterations(config_.max_iterations);
  solver.setTolerance(config_.tolerance);
  
  // 计算对角预处理器
  solver.compute(A);
  
  // 求解
  p = solver.solveWithGuess(b, p);
  
  // 记录求解信息
  last_iterations_ = solver.iterations();
  last_residual_ = solver.error();
  last_converged_ = (solver.info() == Eigen::Success);
  
  if (!last_converged_) {
    std::cerr << "警告：求解未收敛，迭代次数：" << last_iterations_ 
              << "，残差：" << last_residual_ << std::endl;
  }
  
  return last_converged_;
}

} // namespace mps2D
