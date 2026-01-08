#include "config/MPSConfig2D.hpp"
#include "core/Particle.hpp"
#include "core/FileOperator.hpp"
#include "core/TimeStepManager.h"
#include "neighbour_list/NeighborListSearcher.hpp"
#include "lsmps/CorrectiveMatrix.hpp"
#include "explicit_force/ExplicitForce.hpp"
#include "surface_detection/SurfaceDetector.hpp"
#include "PPE/PPEMatrixBuilder.hpp"
#include "PPE/PPESolver.hpp"
#include "correction/Correction.hpp"
#include "core/ErrorHandling.h"
#include "../../third_party/ini/SimpleIni.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <cstring>
#include <csignal>
#include <cstdlib>

// PETSc头文件
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>

using namespace mps2D;

// 全局变量用于信号处理
static volatile bool g_interrupted = false;
static int g_current_iteration = 0;

// 信号处理函数
void signal_handler(int signal) {
  g_interrupted = true;
  std::cerr << "\n\n收到中断信号 (" << signal << ")，正在安全退出..." << std::endl;
  std::cerr << "当前迭代次数: " << g_current_iteration << std::endl;
}

// 强制刷新输出流
inline void flush_log() {
  std::cout.flush();
  std::cerr.flush();
}

// 计算进度百分比
inline double calculate_progress(double current_time, double total_time, int iteration, int max_iterations) {
  double time_progress = (total_time > 0) ? (current_time / total_time * 100.0) : 0.0;
  double iter_progress = (max_iterations > 0) ? (static_cast<double>(iteration) / max_iterations * 100.0) : 0.0;
  return std::max(time_progress, iter_progress);
}

int main(int argc, char* argv[]) {
  // 注册信号处理函数
  signal(SIGINT, signal_handler);   // Ctrl+C
  signal(SIGTERM, signal_handler); // 终止信号
  
  // 在PETSc初始化之前保存命令行参数
  // 因为PetscInitialize会修改argc和argv，移除它识别的命令行参数
  std::string program_name = (argc > 0) ? argv[0] : "MPSBaseline2D";
  std::string config_file = "config.ini";
  if (argc > 1) {
    config_file = argv[1];
  }
  
  std::cout << "========== MPS Baseline 2D 模拟程序 ==========" << std::endl;
  std::cout << "程序启动时间: " << std::chrono::duration_cast<std::chrono::seconds>(
      std::chrono::system_clock::now().time_since_epoch()).count() << std::endl;
  flush_log();
  
  // 初始化PETSc（PPE求解需要）
  PetscErrorCode petsc_err = PetscInitialize(&argc, &argv, NULL, NULL);
  if (petsc_err != 0) {
    std::cerr << "错误：PETSc初始化失败" << std::endl;
    return 1;
  }

  try {
    // ========== 步骤1：读取模拟参数 ==========
    std::cout << "========== 步骤1：读取模拟参数 ==========" << std::endl;
    
    // 验证配置文件路径
    if (config_file.empty() || config_file == "-") {
      std::cerr << "错误：无效的配置文件路径" << std::endl;
      std::cerr << "用法: " << program_name << " [配置文件路径]" << std::endl;
      std::cerr << "默认使用: config.ini" << std::endl;
      PetscFinalize();
      return 1;
    }
    
    // 检查配置文件是否存在
    struct stat file_stat;
    if (stat(config_file.c_str(), &file_stat) != 0) {
      std::cerr << "错误：配置文件不存在: " << config_file << std::endl;
      char cwd[1024];
      if (getcwd(cwd, sizeof(cwd)) != NULL) {
        std::cerr << "当前工作目录: " << cwd << std::endl;
      }
      std::cerr << "请检查文件路径是否正确" << std::endl;
      PetscFinalize();
      return 1;
    }
    
    // 加载配置文件
    SimpleIni ini;
    if (!ini.LoadFile(config_file)) {
      std::cerr << "错误：无法加载配置文件 " << config_file << std::endl;
      std::cerr << "请检查文件格式是否正确" << std::endl;
      PetscFinalize();
      return 1;
    }
    std::cout << "成功加载配置文件: " << config_file << std::endl;
    
    // 创建配置对象并加载参数
    MPSConfig2D config;
    if (!config.LoadFromConfig(ini)) {
      std::cerr << "错误：配置文件参数验证失败" << std::endl;
      PetscFinalize();
      return 1;
    }
    
    // 获取配置参数
    const auto& file_config = config.GetFileConfig();
    const auto& sim_config = config.GetSimulationConfig();
    const auto& particle_config = config.GetParticleConfig();
    
    std::cout << "配置参数加载完成" << std::endl;
    std::cout << "  输入目录: " << file_config.input_dir << std::endl;
    std::cout << "  输出目录: " << file_config.output_dir << std::endl;
    
    // 创建输出目录（如果不存在）
    struct stat info;
    if (stat(file_config.output_dir.c_str(), &info) != 0) {
      // 目录不存在，创建它
      #ifdef _WIN32
        std::string cmd = "mkdir " + file_config.output_dir;
      #else
        std::string cmd = "mkdir -p " + file_config.output_dir;
      #endif
      int result = system(cmd.c_str());
      if (result != 0) {
        std::cerr << "警告：无法创建输出目录 " << file_config.output_dir << std::endl;
      } else {
        std::cout << "已创建输出目录: " << file_config.output_dir << std::endl;
      }
    }
    std::cout << "  时间步长: " << sim_config.time_step << " s" << std::endl;
    std::cout << "  总仿真时间: " << sim_config.total_time << " s" << std::endl;
    std::cout << "  密度: " << sim_config.density << " kg/m³" << std::endl;
    std::cout << "  平滑半径: " << particle_config.smoothing_radius << " m" << std::endl;
    
    // ========== 步骤2：读取前处理文件 ==========
    std::cout << "\n========== 步骤2：读取前处理文件 ==========" << std::endl;
    
    // 创建粒子对象
    FluidParticle fluid_particles("Fluid");
    SolidParticle solid_particles("Solid");
    
    // 读取流体粒子文件
    FileOperator file_operator;
    std::string fluid_file = file_config.input_dir + "/" + file_config.fluid_particle_file;
    int num_fluid = file_operator.getParticleFromFile(fluid_file, fluid_particles);
    if (num_fluid <= 0) {
      std::cerr << "错误：无法读取流体粒子文件 " << fluid_file << std::endl;
      PetscFinalize();
      return 1;
    }
    std::cout << "成功读取流体粒子: " << num_fluid << " 个" << std::endl;
    
    // 初始化流体粒子的其他属性
    fluid_particles.density.resize(num_fluid, sim_config.density);
    fluid_particles.pressure.resize(num_fluid, 0.0);
    fluid_particles.surface_type.resize(num_fluid, SurfaceType::INNER);
    fluid_particles.fluid_neighbour_list.resize(num_fluid);
    fluid_particles.solid_neighbour_list.resize(num_fluid);
    
    // 读取固体粒子文件
    std::string solid_file = file_config.input_dir + "/" + file_config.solid_particle_file;
    int num_solid = file_operator.getParticleFromFile(solid_file, solid_particles);
    if (num_solid < 0) {
      std::cerr << "错误：无法读取固体粒子文件 " << solid_file << std::endl;
      PetscFinalize();
      return 1;
    }
    std::cout << "成功读取固体粒子: " << num_solid << " 个" << std::endl;
    
    // 验证固体粒子数据
    std::cout << "\n========== 验证固体粒子数据 ==========" << std::endl;
    if (num_solid == 0) {
      std::cerr << "警告：没有读取到任何固体粒子！" << std::endl;
    } else {
      // 检查数据一致性
      bool data_valid = true;
      if (static_cast<int>(solid_particles.position.size()) != num_solid) {
        std::cerr << "错误：位置数据数量不匹配！期望 " << num_solid 
                  << "，实际 " << solid_particles.position.size() << std::endl;
        data_valid = false;
      }
      if (static_cast<int>(solid_particles.velocity.size()) != num_solid) {
        std::cerr << "错误：速度数据数量不匹配！期望 " << num_solid 
                  << "，实际 " << solid_particles.velocity.size() << std::endl;
        data_valid = false;
      }
      if (static_cast<int>(solid_particles.normal_vector.size()) != num_solid) {
        std::cerr << "错误：法向量数据数量不匹配！期望 " << num_solid 
                  << "，实际 " << solid_particles.normal_vector.size() << std::endl;
        data_valid = false;
      }
      
      if (data_valid) {
        std::cout << "✓ 数据一致性检查通过" << std::endl;
        
        // 计算位置范围
        double min_x = solid_particles.position[0].x, max_x = solid_particles.position[0].x;
        double min_y = solid_particles.position[0].y, max_y = solid_particles.position[0].y;
        for (int i = 1; i < num_solid; ++i) {
          min_x = std::min(min_x, solid_particles.position[i].x);
          max_x = std::max(max_x, solid_particles.position[i].x);
          min_y = std::min(min_y, solid_particles.position[i].y);
          max_y = std::max(max_y, solid_particles.position[i].y);
        }
        std::cout << "  位置范围: x=[" << std::fixed << std::setprecision(6) 
                  << min_x << ", " << max_x << "], y=[" << min_y << ", " << max_y << "]" << std::endl;
        
        // 检查法向量
        int zero_normal_count = 0;
        double max_normal_magnitude = 0.0;
        for (int i = 0; i < num_solid; ++i) {
          double nx = solid_particles.normal_vector[i].x;
          double ny = solid_particles.normal_vector[i].y;
          double magnitude = std::sqrt(nx * nx + ny * ny);
          max_normal_magnitude = std::max(max_normal_magnitude, magnitude);
          if (magnitude < 1e-10) {
            zero_normal_count++;
          }
        }
        std::cout << "  法向量统计: 最大模长=" << max_normal_magnitude 
                  << ", 零法向量数量=" << zero_normal_count << std::endl;
        if (zero_normal_count > 0) {
          std::cerr << "  警告：有 " << zero_normal_count << " 个粒子的法向量为零或接近零！" << std::endl;
        }
        
        // 检查速度
        int non_zero_velocity_count = 0;
        for (int i = 0; i < num_solid; ++i) {
          double vx = solid_particles.velocity[i].x;
          double vy = solid_particles.velocity[i].y;
          if (std::abs(vx) > 1e-10 || std::abs(vy) > 1e-10) {
            non_zero_velocity_count++;
          }
        }
        std::cout << "  速度统计: 非零速度粒子数量=" << non_zero_velocity_count 
                  << " / " << num_solid << std::endl;
        
        // 输出前5个粒子的详细信息
        int sample_count = std::min(5, num_solid);
        std::cout << "\n  前 " << sample_count << " 个固体粒子的详细信息:" << std::endl;
        for (int i = 0; i < sample_count; ++i) {
          std::cout << "    粒子 " << i << ": pos=(" << std::fixed << std::setprecision(6)
                    << solid_particles.position[i].x << ", " << solid_particles.position[i].y << ")"
                    << ", vel=(" << solid_particles.velocity[i].x << ", " 
                    << solid_particles.velocity[i].y << ")"
                    << ", normal=(" << solid_particles.normal_vector[i].x << ", " 
                    << solid_particles.normal_vector[i].y << ")" << std::endl;
        }
        
        // 输出最后5个粒子的信息
        if (num_solid > 5) {
          std::cout << "\n  最后 " << std::min(5, num_solid - 5) << " 个固体粒子的详细信息:" << std::endl;
          int start_idx = std::max(5, num_solid - 5);
          for (int i = start_idx; i < num_solid; ++i) {
            std::cout << "    粒子 " << i << ": pos=(" << std::fixed << std::setprecision(6)
                      << solid_particles.position[i].x << ", " << solid_particles.position[i].y << ")"
                      << ", vel=(" << solid_particles.velocity[i].x << ", " 
                      << solid_particles.velocity[i].y << ")"
                      << ", normal=(" << solid_particles.normal_vector[i].x << ", " 
                      << solid_particles.normal_vector[i].y << ")" << std::endl;
          }
        }
      } else {
        std::cerr << "错误：固体粒子数据验证失败！" << std::endl;
        PetscFinalize();
        return 1;
      }
    }
    flush_log();
    
    // ========== 输出壁面粒子（仅在程序开始时输出一次）==========
    if (num_solid > 0) {
      std::cout << "\n========== 输出壁面粒子（初始状态）==========" << std::endl;
      std::string solid_output_file = file_config.output_dir + "/solid_initial.vtk";
      
      bool solid_write_success = file_operator.writeVTKBase(solid_output_file, solid_particles);
      if (solid_write_success) {
        // 追加法向量（vector格式）
        file_operator.appendVTKVector(solid_output_file, "normal_vector", solid_particles.normal_vector);
        
        std::cout << "壁面粒子已输出到: " << solid_output_file << std::endl;
        flush_log();
      } else {
        std::cerr << "警告：输出壁面粒子文件失败: " << solid_output_file << std::endl;
        flush_log();
      }
    }
    
    // ========== 初始化时间步管理器 ==========
    std::cout << "\n========== 初始化时间步管理器 ==========" << std::endl;
    TimeStepManager time_manager(config, particle_config.particle_spacing);
    std::cout << "时间步管理器初始化完成" << std::endl;
    
    // ========== 创建模块实例 ==========
    NeighborListSearcher neighbor_searcher;
    CorrectiveMatrix corrective_matrix_calc;
    ExplicitForce explicit_force;
    SurfaceDetector surface_detector;
    PPEMatrixBuilder ppe_matrix_builder;
    PPESolver ppe_solver;
    Correction correction;
    
    // PPE求解器配置
    // 注意：使用罚函数方法将问题转化为对称正定系统 K·p = f，然后使用CG方法求解
    // 与测试程序保持一致的方法
    PPESolver::SolverConfig solver_config;
    solver_config.solver_type = PPESolver::SolverType::CG;  // 使用CG方法
    solver_config.max_iterations = 10000;  // 最大迭代次数（2万粒子可能需要更多迭代）
    solver_config.tolerance = 1e-6;
    solver_config.force_iterative = true;  // 强制使用迭代求解器
    solver_config.is_symmetric_positive_definite = true;  // 使用罚函数方法后的系统是对称正定的
    ppe_solver.SetConfig(solver_config);
    
    // 罚函数参数
    double penalty_parameter = 1e3;  // 罚函数参数μ（与测试程序保持一致）
    
    // ========== 模拟循环 ==========
    std::cout << "\n========== 开始模拟循环 ==========" << std::endl;
    std::cout << "最大迭代次数: " << sim_config.max_iterations << std::endl;
    std::cout << "总仿真时间: " << sim_config.total_time << " s" << std::endl;
    std::cout << "输出间隔: " << sim_config.output_interval << " s" << std::endl;
    flush_log();
    
    int iteration = 0;
    int output_count = 0;
    auto start_time = std::chrono::steady_clock::now();
    auto last_status_time = start_time;
    
    // 用于存储corrective matrix的向量（在循环外声明以避免重复分配）
    std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>>
        corrective_matrices_explicit(num_fluid);
    
    while (!time_manager.IsSimulationFinished() && iteration < sim_config.max_iterations && !g_interrupted) {
      try {
        ++iteration;
        g_current_iteration = iteration;
        double current_time = time_manager.GetCurrentTime();
        double time_step = time_manager.GetTimeStep();
        
        // 计算进度
        double progress = calculate_progress(current_time, sim_config.total_time, iteration, sim_config.max_iterations);
        
        // 计算已用时间
        auto current_wall_time = std::chrono::steady_clock::now();
        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(current_wall_time - start_time).count();
        auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(current_wall_time - start_time).count();
        
        // 每10个时间步或每5秒输出一次详细进度
        bool should_output_detail = (iteration % 10 == 0 || iteration == 1);
        auto time_since_last_status = std::chrono::duration_cast<std::chrono::seconds>(current_wall_time - last_status_time).count();
        bool should_output_status = (time_since_last_status >= 5);
        
        if (should_output_detail || should_output_status) {
          if (should_output_status) {
            last_status_time = current_wall_time;
          }
          std::cout << "\n--- 时间步 " << iteration 
                    << " | 时间: " << std::fixed << std::setprecision(6) << current_time 
                    << " s / " << sim_config.total_time << " s"
                    << " | 进度: " << std::setprecision(2) << progress << "%"
                    << " | 时间步长: " << std::setprecision(6) << time_step << " s"
                    << " | 已用时间: " << elapsed_seconds << " s ---" << std::endl;
          flush_log();
        }
      
        // 构建邻域列表（在每次循环开始时执行）
        auto step_start = std::chrono::steady_clock::now();
        if (should_output_detail) {
          std::cout << "  [步骤1/8] 构建邻域列表..." << std::flush;
        }
        try {
          neighbor_searcher.BuildNeighborList(
              fluid_particles, solid_particles,
              particle_config.particle_radius,
              particle_config.smoothing_radius,
              particle_config.cell_size);
          if (should_output_detail) {
            auto step_end = std::chrono::steady_clock::now();
            auto step_duration = std::chrono::duration_cast<std::chrono::milliseconds>(step_end - step_start).count();
            std::cout << " 完成 (" << step_duration << " ms)" << std::endl;
            flush_log();
          }
        } catch (const std::exception& e) {
          std::cerr << "\n错误：构建邻域列表失败 - " << e.what() << std::endl;
          flush_log();
          throw;
        }
        
        // 表面检测（在每次循环开始时执行）
        step_start = std::chrono::steady_clock::now();
        if (should_output_detail) {
          std::cout << "  [步骤2/8] 表面检测..." << std::flush;
        }
        try {
          surface_detector.DetectSurfaceParticles(
              fluid_particles, solid_particles,
              particle_config.smoothing_radius,
              particle_config.particle_spacing);
          if (should_output_detail) {
            auto step_end = std::chrono::steady_clock::now();
            auto step_duration = std::chrono::duration_cast<std::chrono::milliseconds>(step_end - step_start).count();
            std::cout << " 完成 (" << step_duration << " ms)" << std::endl;
            flush_log();
          }
        } catch (const std::exception& e) {
          std::cerr << "\n错误：表面检测失败 - " << e.what() << std::endl;
          flush_log();
          throw;
        }
        
        // 计算corrective matrix（用于显式力，第一类边界条件）
        step_start = std::chrono::steady_clock::now();
        if (should_output_detail) {
          std::cout << "  [步骤3/8] 计算corrective matrix（显式力用）..." << std::flush;
        }
        try {
          for (int i = 0; i < num_fluid; ++i) {
            corrective_matrices_explicit[i] = corrective_matrix_calc.ComputeCorrectiveMatrix(
                i, fluid_particles, solid_particles,
                particle_config.smoothing_radius, false);  // 第一类边界条件
          }
          if (should_output_detail) {
            auto step_end = std::chrono::steady_clock::now();
            auto step_duration = std::chrono::duration_cast<std::chrono::milliseconds>(step_end - step_start).count();
            std::cout << " 完成 (" << step_duration << " ms)" << std::endl;
            flush_log();
          }
        } catch (const std::exception& e) {
          std::cerr << "\n错误：计算corrective matrix失败 - " << e.what() << std::endl;
          flush_log();
          throw;
        }
        
        // 步骤4：显式更新模块（只更新速度作为临时速度，不更新位置）
        step_start = std::chrono::steady_clock::now();
        if (should_output_detail) {
          std::cout << "  [步骤4/8] 显式更新模块：计算粘性力和重力，更新临时速度..." << std::flush;
        }
        // 用于保存显式更新后的速度（临时速度）
        std::vector<double2> velocity_explicit(num_fluid);
        // 用于保存粘性力加速度（从explicit_force模块中获取）
        std::vector<double2> viscous_acceleration(num_fluid);
        try {
          // 使用重载方法，直接获取explicit_force模块计算出的粘性力加速度
          explicit_force.ComputeAndUpdateVelocity(
              fluid_particles, solid_particles,
              corrective_matrices_explicit,
              particle_config.smoothing_radius,
              sim_config.kinematic_viscosity,
              sim_config.gravity_x, sim_config.gravity_y,
              time_step,
              viscous_acceleration);
          
          // 保存显式更新后的速度
          for (int i = 0; i < num_fluid; ++i) {
            velocity_explicit[i] = fluid_particles.velocity[i];
          }
          if (should_output_detail) {
            auto step_end = std::chrono::steady_clock::now();
            auto step_duration = std::chrono::duration_cast<std::chrono::milliseconds>(step_end - step_start).count();
            std::cout << " 完成 (" << step_duration << " ms)" << std::endl;
            flush_log();
          }
        } catch (const std::exception& e) {
          std::cerr << "\n错误：显式更新失败 - " << e.what() << std::endl;
          flush_log();
          throw;
        }
      
        // 步骤5：计算corrective matrix（用于PPE）
        // 注意：显式更新只更新速度，不更新位置，因此用于速度散度的corrective matrix
        // 与显式力计算时使用的相同，可以直接复用corrective_matrices_explicit
        // 只需要计算用于压力梯度的corrective matrix（第二类边界条件）
        step_start = std::chrono::steady_clock::now();
        if (should_output_detail) {
          std::cout << "  [步骤5/8] 计算corrective matrix（PPE用，压力梯度）..." << std::flush;
        }
        std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>>
            corrective_matrices_ppe_pressure(num_fluid);  // 用于压力梯度（第二类边界条件）
        
        try {
          for (int i = 0; i < num_fluid; ++i) {
            // 压力梯度使用第二类边界条件
            corrective_matrices_ppe_pressure[i] = corrective_matrix_calc.ComputeCorrectiveMatrix(
                i, fluid_particles, solid_particles,
                particle_config.smoothing_radius, true);
          }
          if (should_output_detail) {
            auto step_end = std::chrono::steady_clock::now();
            auto step_duration = std::chrono::duration_cast<std::chrono::milliseconds>(step_end - step_start).count();
            std::cout << " 完成 (" << step_duration << " ms)" << std::endl;
            flush_log();
          }
        } catch (const std::exception& e) {
          std::cerr << "\n错误：计算PPE用corrective matrix失败 - " << e.what() << std::endl;
          flush_log();
          throw;
        }
        
        // 复用显式力计算时的corrective matrix作为速度散度用（第一类边界条件）
        // 因为显式更新不改变粒子位置，所以可以直接复用
        const auto& corrective_matrices_ppe_velocity = corrective_matrices_explicit;
      
        // 步骤6：构建PPE
        step_start = std::chrono::steady_clock::now();
        if (should_output_detail) {
          std::cout << "  [步骤6/8] 构建PPE系数矩阵和右边项..." << std::flush;
        }
        Mat A_petsc = NULL;
        Vec b_petsc = NULL;
        
        try {
          bool build_success = ppe_matrix_builder.BuildPPEMatrixPetsc(
              fluid_particles, solid_particles,
              corrective_matrices_ppe_velocity,  // 用于速度散度
              corrective_matrices_ppe_pressure,  // 用于压力边界条件
              particle_config.smoothing_radius,
              sim_config.density,
              time_step,
              particle_config.particle_spacing,
              sim_config.gravity_x, sim_config.gravity_y,
              A_petsc, b_petsc);
          
          if (!build_success) {
            std::cerr << "\n错误：PPE矩阵构建失败（时间步 " << iteration << "）" << std::endl;
            flush_log();
            if (A_petsc != NULL) MatDestroy(&A_petsc);
            if (b_petsc != NULL) VecDestroy(&b_petsc);
            throw std::runtime_error("PPE矩阵构建失败");
          }
          if (should_output_detail) {
            auto step_end = std::chrono::steady_clock::now();
            auto step_duration = std::chrono::duration_cast<std::chrono::milliseconds>(step_end - step_start).count();
            std::cout << " 完成 (" << step_duration << " ms)" << std::endl;
            flush_log();
          }
        } catch (const std::exception& e) {
          std::cerr << "\n错误：构建PPE矩阵时发生异常 - " << e.what() << std::endl;
          flush_log();
          if (A_petsc != NULL) MatDestroy(&A_petsc);
          if (b_petsc != NULL) VecDestroy(&b_petsc);
          throw;
        }
        
        // 使用罚函数方法构建对称正定系统 K·p = f
        // K = A^T A + D，f = A^T b（与测试程序保持一致）
        step_start = std::chrono::steady_clock::now();
        if (should_output_detail) {
          std::cout << "  [步骤7/8] 构建罚函数系统（对称正定）..." << std::flush;
        }
        Mat K_petsc = NULL;
        Vec f_petsc = NULL;
        
        try {
          bool penalty_success = ppe_matrix_builder.BuildPenaltySystem(
              A_petsc, b_petsc, fluid_particles, penalty_parameter, K_petsc, f_petsc);
          
          if (!penalty_success) {
            std::cerr << "\n错误：构建罚函数系统失败（时间步 " << iteration << "）" << std::endl;
            flush_log();
            if (A_petsc != NULL) MatDestroy(&A_petsc);
            if (b_petsc != NULL) VecDestroy(&b_petsc);
            if (K_petsc != NULL) MatDestroy(&K_petsc);
            if (f_petsc != NULL) VecDestroy(&f_petsc);
            throw std::runtime_error("构建罚函数系统失败");
          }
          if (should_output_detail) {
            auto step_end = std::chrono::steady_clock::now();
            auto step_duration = std::chrono::duration_cast<std::chrono::milliseconds>(step_end - step_start).count();
            std::cout << " 完成 (" << step_duration << " ms)" << std::endl;
            flush_log();
          }
        } catch (const std::exception& e) {
          std::cerr << "\n错误：构建罚函数系统时发生异常 - " << e.what() << std::endl;
          flush_log();
          if (A_petsc != NULL) MatDestroy(&A_petsc);
          if (b_petsc != NULL) VecDestroy(&b_petsc);
          if (K_petsc != NULL) MatDestroy(&K_petsc);
          if (f_petsc != NULL) VecDestroy(&f_petsc);
          throw;
        }
        
        // 步骤7：求解PPE（使用对称正定系统 K·p = f）
        step_start = std::chrono::steady_clock::now();
        if (should_output_detail) {
          std::cout << "  [步骤8/8] 求解PPE（CG方法，对称正定系统）..." << std::flush;
        }
        Vec p_petsc = NULL;
        bool solve_success = false;
        
        try {
          solve_success = ppe_solver.Solve(K_petsc, f_petsc, p_petsc);
          
          if (!solve_success) {
            std::cerr << "\n警告：PPE求解未收敛（时间步 " << iteration 
                      << ", 迭代次数: " << ppe_solver.GetLastIterations() 
                      << ", 残差: " << ppe_solver.GetLastResidual() << "）" << std::endl;
            flush_log();
            // 即使求解失败，p_petsc可能已经被创建，需要清理
            if (p_petsc != NULL) {
              VecDestroy(&p_petsc);
              p_petsc = NULL;
            }
            // 清理其他PETSc对象
            if (A_petsc != NULL) MatDestroy(&A_petsc);
            if (b_petsc != NULL) VecDestroy(&b_petsc);
            if (K_petsc != NULL) MatDestroy(&K_petsc);
            if (f_petsc != NULL) VecDestroy(&f_petsc);
            // 继续执行，不退出（允许程序继续运行）
          } else {
            if (should_output_detail) {
              auto step_end = std::chrono::steady_clock::now();
              auto step_duration = std::chrono::duration_cast<std::chrono::milliseconds>(step_end - step_start).count();
              std::cout << " 完成 (" << step_duration << " ms, 迭代: " 
                        << ppe_solver.GetLastIterations() 
                        << ", 残差: " << std::scientific << std::setprecision(3)
                        << ppe_solver.GetLastResidual() << std::fixed << ")" << std::endl;
              flush_log();
            }
          }
        } catch (const std::exception& e) {
          std::cerr << "\n错误：求解PPE时发生异常 - " << e.what() << std::endl;
          flush_log();
          if (p_petsc != NULL) VecDestroy(&p_petsc);
          if (A_petsc != NULL) MatDestroy(&A_petsc);
          if (b_petsc != NULL) VecDestroy(&b_petsc);
          if (K_petsc != NULL) MatDestroy(&K_petsc);
          if (f_petsc != NULL) VecDestroy(&f_petsc);
          throw;
        }
        
        // 从PETSc向量提取压力值（只有在求解成功且p_petsc不为NULL时）
        if (solve_success && p_petsc != NULL) {
          PetscScalar* p_array;
          VecGetArray(p_petsc, &p_array);
          for (int i = 0; i < num_fluid; ++i) {
            fluid_particles.pressure[i] = p_array[i];
          }
          VecRestoreArray(p_petsc, &p_array);
        } else if (!solve_success) {
          // 如果求解失败，保持压力为0或使用上一次的值
          // 这里可以选择跳过压力修正步骤
        }
        
        // 步骤8：Correction模块
        step_start = std::chrono::steady_clock::now();
        if (should_output_detail) {
          std::cout << "  [步骤9/9] 压力修正：更新速度和位置..." << std::flush;
        }
        
        // 保存压力梯度（在位置更新之前计算）
        std::vector<double2> pressure_gradients_for_output;
        try {
          correction.ComputeAndUpdateAllParticles(
              fluid_particles, solid_particles,
              corrective_matrices_ppe_pressure,  // 使用第二类边界条件的corrective matrix
              particle_config.smoothing_radius,
              sim_config.gravity_x, sim_config.gravity_y,
              sim_config.density,
              time_step,
              pressure_gradients_for_output);  // 保存压力梯度
          
          if (should_output_detail) {
            auto step_end = std::chrono::steady_clock::now();
            auto step_duration = std::chrono::duration_cast<std::chrono::milliseconds>(step_end - step_start).count();
            std::cout << " 完成 (" << step_duration << " ms)" << std::endl;
            flush_log();
          }
        } catch (const std::exception& e) {
          std::cerr << "\n错误：压力修正失败 - " << e.what() << std::endl;
          flush_log();
          throw;
        }
        
        // 注意：corrective_matrices_explicit在下次循环开始时会重新计算，不需要复制
        
        // 清理PETSc对象（确保所有对象都被正确清理）
        // 注意：按照创建顺序的逆序销毁，避免依赖问题
        // 如果求解失败，这些对象可能已经被清理，需要检查
        if (p_petsc != NULL) {
          VecDestroy(&p_petsc);
          p_petsc = NULL;
        }
        if (f_petsc != NULL) {
          VecDestroy(&f_petsc);
          f_petsc = NULL;
        }
        if (K_petsc != NULL) {
          MatDestroy(&K_petsc);
          K_petsc = NULL;
        }
        if (b_petsc != NULL) {
          VecDestroy(&b_petsc);
          b_petsc = NULL;
        }
        if (A_petsc != NULL) {
          MatDestroy(&A_petsc);
          A_petsc = NULL;
        }
        
        // 动态调整时间步（时间步会自动更新到time_manager中）
        try {
          time_manager.AdjustTimeStep(fluid_particles);
        } catch (const std::exception& e) {
          std::cerr << "\n警告：调整时间步失败 - " << e.what() << "，继续使用当前时间步" << std::endl;
          flush_log();
        }
        
        // 步骤9：判断是否需要输出计算结果
        if (time_manager.ShouldOutput()) {
          ++output_count;
          std::cout << "\n  [输出] 输出计算结果（第 " << output_count << " 次）..." << std::flush;
          flush_log();
          
          // 构建输出文件名（不包含时间）
          std::ostringstream oss;
          oss << file_config.output_dir << "/result_" 
              << std::setfill('0') << std::setw(6) << output_count << ".vtk";
          std::string output_file = oss.str();
          
          try {
            // 输出流体粒子VTK文件
            bool write_success = file_operator.writeVTKBase(output_file, fluid_particles);
            if (write_success) {
              // 追加压力和密度数据（标量）
              file_operator.appendVTKScalar(output_file, "pressure", fluid_particles.pressure);
              file_operator.appendVTKScalar(output_file, "density", fluid_particles.density);
              
              // 追加显式更新后的速度向量（vector格式）
              file_operator.appendVTKVector(output_file, "velocity_explicit", velocity_explicit);
              
              // 追加correction后的速度向量（vector格式）
              file_operator.appendVTKVector(output_file, "velocity_corrected", fluid_particles.velocity);
              
              // 追加压力梯度向量（vector格式）- 使用correction模块中计算的压力梯度
              file_operator.appendVTKVector(output_file, "pressure_gradient", pressure_gradients_for_output);
              
              // 追加表面类型（作为标量）
              std::vector<int> surface_type_int(num_fluid);
              for (int i = 0; i < num_fluid; ++i) {
                surface_type_int[i] = static_cast<int>(fluid_particles.surface_type[i]);
              }
              file_operator.appendVTKScalar(output_file, "surface_type", surface_type_int);
              
              // 追加粘性力向量（vector格式）- 使用explicit_force模块计算出的粘性力加速度
              file_operator.appendVTKVector(output_file, "viscous_force", viscous_acceleration);
              
              std::cout << " 完成: " << output_file << std::endl;
              flush_log();
            } else {
              std::cerr << "\n警告：输出流体粒子文件失败: " << output_file << std::endl;
              flush_log();
            }
          } catch (const std::exception& e) {
            std::cerr << "\n错误：输出文件时发生异常 - " << e.what() << std::endl;
            flush_log();
            // 不抛出异常，继续执行
          }
          
          time_manager.MarkOutputDone();
        }
        
        // 更新时间
        try {
          time_manager.AdvanceTime();
        } catch (const std::exception& e) {
          std::cerr << "\n错误：推进时间失败 - " << e.what() << std::endl;
          flush_log();
          throw;
        }
      
        // 判断是否结束
        if (g_interrupted) {
          std::cout << "\n\n========== 模拟被用户中断 ==========" << std::endl;
          std::cout << "中断时间步: " << iteration << std::endl;
          std::cout << "中断时间: " << std::fixed << std::setprecision(6) << current_time << " s" << std::endl;
          flush_log();
          break;
        }
        
        if (time_manager.IsSimulationFinished()) {
          std::cout << "\n\n========== 模拟已完成 ==========" << std::endl;
          std::cout << "完成原因: 达到总时间 " << sim_config.total_time << " s" << std::endl;
          std::cout << "最终时间: " << std::fixed << std::setprecision(6) << current_time << " s" << std::endl;
          flush_log();
          break;
        }
        
        if (iteration >= sim_config.max_iterations) {
          std::cout << "\n\n========== 模拟已结束 ==========" << std::endl;
          std::cout << "完成原因: 达到最大迭代次数 " << sim_config.max_iterations << std::endl;
          std::cout << "最终时间: " << std::fixed << std::setprecision(6) << current_time << " s" << std::endl;
          flush_log();
          break;
        }
        
      } catch (const std::exception& e) {
        std::cerr << "\n\n========== 模拟循环中发生异常 ==========" << std::endl;
        std::cerr << "异常时间步: " << iteration << std::endl;
        std::cerr << "异常时间: " << std::fixed << std::setprecision(6) 
                  << time_manager.GetCurrentTime() << " s" << std::endl;
        std::cerr << "异常信息: " << e.what() << std::endl;
        flush_log();
        throw;  // 重新抛出异常，由外层catch处理
      }
    }
    
    // ========== 模拟结束 ==========
    auto end_time = std::chrono::steady_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    double final_time = time_manager.GetCurrentTime();
    
    std::cout << "\n========== 模拟结束统计 ==========" << std::endl;
    std::cout << "总迭代次数: " << iteration << " / " << sim_config.max_iterations << std::endl;
    std::cout << "最终模拟时间: " << std::fixed << std::setprecision(6) << final_time 
              << " s / " << sim_config.total_time << " s" << std::endl;
    std::cout << "总计算时间: " << elapsed_seconds << " 秒 (" 
              << (elapsed_milliseconds / 1000.0) << " 秒)" << std::endl;
    if (iteration > 0) {
      std::cout << "平均每个时间步: " << std::fixed << std::setprecision(3) 
                << (elapsed_milliseconds / static_cast<double>(iteration)) << " 毫秒" << std::endl;
      std::cout << "模拟速度: " << std::fixed << std::setprecision(2) 
                << (final_time / (elapsed_milliseconds / 1000.0)) << " 倍实时速度" << std::endl;
    }
    std::cout << "输出文件数: " << output_count << std::endl;
    std::cout << "完成状态: " << (g_interrupted ? "用户中断" : 
                                  (time_manager.IsSimulationFinished() ? "正常完成" : 
                                   (iteration >= sim_config.max_iterations ? "达到最大迭代次数" : "未知"))) << std::endl;
    flush_log();
    
  } catch (const MPSException& e) {
    std::cerr << "错误: " << e.what() << std::endl;
    PetscFinalize();
    return 1;
  } catch (const std::exception& e) {
    std::cerr << "标准异常: " << e.what() << std::endl;
    PetscFinalize();
    return 1;
  } catch (...) {
    std::cerr << "未知异常" << std::endl;
    PetscFinalize();
    return 1;
  }
  
  // 清理PETSc
  PetscFinalize();
  return 0;
}

