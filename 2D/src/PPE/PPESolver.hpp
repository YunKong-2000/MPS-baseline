#pragma once

// PETSc前向声明
struct _p_Mat;
struct _p_Vec;
typedef struct _p_Mat* Mat;
typedef struct _p_Vec* Vec;

namespace mps2D {

// PPE（压力泊松方程）迭代求解器
// 使用PETSc库求解大规模非对称稀疏线性方程组 Ap = b
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
  ~PPESolver();

  // 求解PPE方程 Ap = b（使用PETSc矩阵和向量）
  // 参数：
  //   A_petsc: PETSc系数矩阵
  //   b_petsc: PETSc右边项向量
  //   p_petsc: 输出的PETSc压力解向量（如果为NULL，函数会创建；如果非NULL，会使用现有向量）
  // 返回：是否成功求解
  // 注意：此方法不会销毁输入的A_petsc和b_petsc，但会创建p_petsc（如果为NULL）
  //       调用者负责管理这些PETSc对象的内存
  bool Solve(
      Mat A_petsc,
      Vec b_petsc,
      Vec& p_petsc);

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
};

} // namespace mps2D
