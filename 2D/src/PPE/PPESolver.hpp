#pragma once
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <string>

namespace mps2D {

// PPE（压力泊松方程）迭代求解器
// 用于求解大规模非对称稀疏线性方程组 Ap = b
class PPESolver {
public:
  // 求解器类型
  enum class SolverType {
    BICGSTAB,  // BiCGSTAB迭代法（推荐用于非对称矩阵）
    GMRES      // GMRES迭代法
  };

  // 求解器配置参数
  struct SolverConfig {
    SolverType solver_type = SolverType::BICGSTAB;
    int max_iterations = 1000;           // 最大迭代次数
    double tolerance = 1e-6;              // 收敛容差
    int restart = 30;                    // GMRES重启参数（仅用于GMRES）
  };

  PPESolver() = default;
  explicit PPESolver(const SolverConfig& config);
  ~PPESolver() = default;

  // 求解PPE方程 Ap = b
  // 参数：
  //   A: 系数矩阵（稀疏矩阵，压缩格式CSC）
  //   b: 右边项向量
  //   p: 输出的压力解向量
  // 返回：是否成功求解
  bool Solve(
      const Eigen::SparseMatrix<double, Eigen::ColMajor>& A,
      const Eigen::VectorXd& b,
      Eigen::VectorXd& p);

  // 设置求解器配置
  void SetConfig(const SolverConfig& config) { config_ = config; }

  // 获取求解器配置
  const SolverConfig& GetConfig() const { return config_; }

  // 获取最后一次求解的迭代次数
  int GetLastIterations() const { return last_iterations_; }

  // 获取最后一次求解的残差
  double GetLastResidual() const { return last_residual_; }

  // 获取最后一次求解是否收敛
  bool GetLastConverged() const { return last_converged_; }

private:
  SolverConfig config_;
  int last_iterations_ = 0;
  double last_residual_ = 0.0;
  bool last_converged_ = false;

  // 使用BiCGSTAB方法求解
  bool SolveBiCGSTAB(
      const Eigen::SparseMatrix<double, Eigen::ColMajor>& A,
      const Eigen::VectorXd& b,
      Eigen::VectorXd& p);

  // 使用GMRES方法求解
  bool SolveGMRES(
      const Eigen::SparseMatrix<double, Eigen::ColMajor>& A,
      const Eigen::VectorXd& b,
      Eigen::VectorXd& p);
};

} // namespace mps2D
