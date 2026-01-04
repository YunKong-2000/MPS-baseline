#include "../src/PPE/PPEMatrixBuilder.hpp"
#include "../src/PPE/PPESolver.hpp"
#include "../src/lsmps/CorrectiveMatrix.hpp"
#include "../src/neighbour_list/NeighborListSearcher.hpp"
#include "../src/surface_detection/SurfaceDetector.hpp"
#include "../src/core/Particle.hpp"
#include "../src/core/FileOperator.hpp"
#include "core/Types.h"
#include "core/MPSUtils.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <map>
#include <chrono>
#include <fstream>
#include <Eigen/Dense>

// PETSc头文件
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>

using namespace mps2D;

// 生成不规则几何壁面的静水压力测试场景
// 参数：
//   container_width: 容器宽度
//   container_height: 容器高度
//   water_height: 水位高度
//   particle_spacing: 粒子间距
//   fluid_particles: 输出流体粒子
//   solid_particles: 输出固体粒子
//   wall_type: 壁面类型 ("step": 阶梯形, "wave": 波浪形, "slope": 倾斜形, "polygon": 多边形)
void GenerateIrregularHydrostaticTest(
    double container_width,
    double container_height,
    double water_height,
    double particle_spacing,
    const std::string& wall_type,
    FluidParticle& fluid_particles,
    SolidParticle& solid_particles) {
  
  // 计算流体粒子区域
  double center_x = container_width / 2.0;
  double min_x = 0.0;
  double max_x = container_width;
  double min_y = 0.0;
  double max_y = water_height;
  
  // 计算x方向的粒子数
  int nx_fluid = static_cast<int>(std::round(container_width / particle_spacing)) + 1;
  int ny_fluid = static_cast<int>(std::round(water_height / particle_spacing)) + 1;
  
  // 生成流体粒子
  std::vector<mps2D::double2> fluid_positions;
  for (int j = 0; j < ny_fluid; ++j) {
    double y = j * particle_spacing;
    if (y > max_y) break;
    for (int i = 0; i < nx_fluid; ++i) {
      double x = i * particle_spacing;
      if (x > max_x) break;
      
      // 检查是否在有效区域内（避开不规则壁面）
      bool inside = true;
      
      // 根据壁面类型检查
      if (wall_type == "step") {
        // 阶梯形：左侧有阶梯
        double step_x1 = container_width * 0.2;
        double step_x2 = container_width * 0.4;
        double step_y1 = water_height * 0.3;
        if (x < step_x1 && y < step_y1) {
          inside = false;  // 在阶梯内部
        }
        if (x >= step_x1 && x < step_x2 && y < water_height * 0.5) {
          inside = false;  // 在阶梯内部
        }
      } else if (wall_type == "wave") {
        // 波浪形：右侧壁面是波浪形
        // 流体粒子应该在波浪壁面内侧，距离壁面一个粒子间距
        double wave_amplitude = container_width * 0.05;
        double wave_frequency = 2.0 * M_PI / (container_height * 0.5);
        double wave_x = container_width - wave_amplitude * (1.0 + std::sin(wave_frequency * y)) - particle_spacing;
        if (x > wave_x) {
          inside = false;  // 在波浪壁面外侧
        }
      } else if (wall_type == "slope") {
        // 倾斜形：底部倾斜
        double slope_angle = M_PI / 6.0;  // 30度
        double slope_y = x * std::tan(slope_angle);
        if (y < slope_y) {
          inside = false;  // 在倾斜底部下方
        }
      } else if (wall_type == "polygon") {
        // 多边形：左侧是五边形
        double poly_center_x = container_width * 0.15;
        double poly_center_y = water_height * 0.4;
        double poly_radius = container_width * 0.12;
        int num_sides = 5;
        
        // 检查是否在多边形内部
        double dx = x - poly_center_x;
        double dy = y - poly_center_y;
        double dist = std::sqrt(dx * dx + dy * dy);
        if (dist < poly_radius) {
          // 进一步检查角度
          double angle = std::atan2(dy, dx);
          if (angle < 0) angle += 2.0 * M_PI;
          double sector_angle = 2.0 * M_PI / num_sides;
          int sector = static_cast<int>(angle / sector_angle);
          double sector_center = sector * sector_angle + sector_angle / 2.0;
          double dist_to_edge = dist * std::cos(angle - sector_center);
          if (dist_to_edge < poly_radius * 0.9) {
            inside = false;
          }
        }
      }
      
      if (inside) {
        fluid_positions.push_back({x, y});
      }
    }
  }
  
  int num_fluid = static_cast<int>(fluid_positions.size());
  fluid_particles.particle_num = num_fluid;
  fluid_particles.position.resize(num_fluid);
  fluid_particles.velocity.resize(num_fluid);
  fluid_particles.pressure.resize(num_fluid);
  fluid_particles.density.resize(num_fluid);
  fluid_particles.surface_type.resize(num_fluid);
  fluid_particles.fluid_neighbour_list.resize(num_fluid);
  fluid_particles.solid_neighbour_list.resize(num_fluid);
  
  for (int i = 0; i < num_fluid; ++i) {
    fluid_particles.position[i] = fluid_positions[i];
    fluid_particles.velocity[i] = {0.0, 0.0};
    fluid_particles.pressure[i] = 0.0;
    fluid_particles.density[i] = 1000.0;
    fluid_particles.surface_type[i] = SurfaceType::INNER;
  }
  
  // 生成不规则壁面粒子
  std::vector<mps2D::double2> solid_positions;
  std::vector<mps2D::double2> solid_normals;
  
  // 底部壁面（根据类型可能倾斜）
  if (wall_type == "slope") {
    // 倾斜底部
    double slope_angle = M_PI / 6.0;
    for (int i = 0; i <= nx_fluid; ++i) {
      double x = i * particle_spacing;
      if (x > max_x) break;
      double y = x * std::tan(slope_angle) - particle_spacing;
      solid_positions.push_back({x, y});
      // 法向量垂直于倾斜面
      double nx = std::sin(slope_angle);
      double ny = std::cos(slope_angle);
      solid_normals.push_back({nx, ny});
    }
  } else {
    // 平底
    for (int i = 0; i <= nx_fluid; ++i) {
      double x = i * particle_spacing;
      if (x > max_x) break;
      solid_positions.push_back({x, -particle_spacing});
      solid_normals.push_back({0.0, 1.0});
    }
  }
  
  // 左侧壁面（根据类型）
  if (wall_type == "step") {
    // 阶梯形左侧壁面
    double step_x1 = container_width * 0.2;
    double step_x2 = container_width * 0.4;
    double step_y1 = water_height * 0.3;
    double step_y2 = water_height * 0.5;
    
    // 第一段（底部到第一个阶梯）
    for (int j = 0; j * particle_spacing < step_y1; ++j) {
      double y = j * particle_spacing;
      solid_positions.push_back({-particle_spacing, y});
      solid_normals.push_back({1.0, 0.0});
    }
    // 第二段（第一个阶梯平台）
    for (int i = 0; i * particle_spacing < step_x1; ++i) {
      double x = i * particle_spacing;
      solid_positions.push_back({x, step_y1 - particle_spacing});
      solid_normals.push_back({0.0, 1.0});
    }
    // 第三段（第一个阶梯到第二个阶梯）
    for (int j = 0; step_y1 + j * particle_spacing < step_y2; ++j) {
      double y = step_y1 + j * particle_spacing;
      solid_positions.push_back({step_x1 - particle_spacing, y});
      solid_normals.push_back({1.0, 0.0});
    }
    // 第四段（第二个阶梯平台）
    for (int i = 0; step_x1 + i * particle_spacing < step_x2; ++i) {
      double x = step_x1 + i * particle_spacing;
      solid_positions.push_back({x, step_y2 - particle_spacing});
      solid_normals.push_back({0.0, 1.0});
    }
    // 第五段（第二个阶梯到顶部）
    for (int j = 0; step_y2 + j * particle_spacing < water_height + particle_spacing; ++j) {
      double y = step_y2 + j * particle_spacing;
      solid_positions.push_back({step_x2 - particle_spacing, y});
      solid_normals.push_back({1.0, 0.0});
    }
  } else if (wall_type == "polygon") {
    // 多边形左侧：在五边形周围生成壁面粒子
    double poly_center_x = container_width * 0.15;
    double poly_center_y = water_height * 0.4;
    double poly_radius = container_width * 0.12;
    int num_sides = 5;
    
    // 生成五边形周围的壁面粒子
    for (int side = 0; side < num_sides; ++side) {
      double angle1 = 2.0 * M_PI * side / num_sides;
      double angle2 = 2.0 * M_PI * (side + 1) / num_sides;
      int num_particles = static_cast<int>(poly_radius * std::abs(angle2 - angle1) / particle_spacing) + 1;
      for (int k = 0; k < num_particles; ++k) {
        double t = static_cast<double>(k) / num_particles;
        double angle = angle1 + t * (angle2 - angle1);
        double x = poly_center_x + (poly_radius + particle_spacing) * std::cos(angle);
        double y = poly_center_y + (poly_radius + particle_spacing) * std::sin(angle);
        if (x < 0 || y < 0 || y > water_height) continue;
        solid_positions.push_back({x, y});
        // 法向量指向多边形中心
        double nx = -std::cos(angle);
        double ny = -std::sin(angle);
        solid_normals.push_back({nx, ny});
      }
    }
    
    // 左侧其余部分
    for (int j = 0; j * particle_spacing < water_height + particle_spacing; ++j) {
      double y = j * particle_spacing;
      double x = -particle_spacing;
      // 检查是否与多边形重叠
      bool overlap = false;
      for (const auto& pos : solid_positions) {
        if (std::abs(pos.x - x) < particle_spacing * 0.5 && 
            std::abs(pos.y - y) < particle_spacing * 0.5) {
          overlap = true;
          break;
        }
      }
      if (!overlap) {
        solid_positions.push_back({x, y});
        solid_normals.push_back({1.0, 0.0});
      }
    }
  } else {
    // 普通左侧壁面
    for (int j = 0; j * particle_spacing < water_height + particle_spacing; ++j) {
      double y = j * particle_spacing;
      solid_positions.push_back({-particle_spacing, y});
      solid_normals.push_back({1.0, 0.0});
    }
  }
  
  // 右侧壁面（根据类型）
  if (wall_type == "wave") {
    // 波浪形右侧壁面
    // 壁面粒子应该在波浪表面，与流体粒子距离一个粒子间距
    double wave_amplitude = container_width * 0.05;
    double wave_frequency = 2.0 * M_PI / (container_height * 0.5);
    for (int j = 0; j * particle_spacing < water_height + particle_spacing; ++j) {
      double y = j * particle_spacing;
      // 波浪表面的x坐标
      double wave_x = container_width - wave_amplitude * (1.0 + std::sin(wave_frequency * y));
      // 壁面粒子位置：在波浪表面外侧一个粒子间距
      solid_positions.push_back({wave_x + particle_spacing, y});
      // 法向量垂直于波浪表面（指向流体）
      double dx = -wave_amplitude * wave_frequency * std::cos(wave_frequency * y);
      double norm = std::sqrt(1.0 + dx * dx);
      solid_normals.push_back({-1.0 / norm, dx / norm});
    }
  } else {
    // 普通右侧壁面：与流体粒子距离一个粒子间距
    for (int j = 0; j * particle_spacing < water_height + particle_spacing; ++j) {
      double y = j * particle_spacing;
      solid_positions.push_back({container_width + particle_spacing, y});
      solid_normals.push_back({-1.0, 0.0});
    }
  }
  
  int num_solid = static_cast<int>(solid_positions.size());
  solid_particles.particle_num = num_solid;
  solid_particles.position.resize(num_solid);
  solid_particles.velocity.resize(num_solid);
  solid_particles.normal_vector.resize(num_solid);
  
  for (int i = 0; i < num_solid; ++i) {
    solid_particles.position[i] = solid_positions[i];
    solid_particles.velocity[i] = {0.0, 0.0};
    solid_particles.normal_vector[i] = solid_normals[i];
  }
  
  std::cout << "生成不规则几何壁面静水测试场景：" << std::endl;
  std::cout << "  壁面类型: " << wall_type << std::endl;
  std::cout << "  流体粒子数: " << num_fluid << std::endl;
  std::cout << "  壁面粒子数: " << num_solid << std::endl;
}

int main(int argc, char** argv) {
  // 初始化PETSc
  PetscInitialize(&argc, &argv, NULL, NULL);
  
  // 测试参数
  const double particle_spacing = 0.02;  // 粒子间距 (m)
  const double container_width = 1.0;   // 容器宽度 (m)
  const double container_height = 1.0;   // 容器高度 (m)
  const double water_height = 0.6;      // 水位高度 (m)
  const double g = 9.81;                // 重力加速度 (m/s²)
  const double rho0 = 1000.0;           // 参考密度 (kg/m³)
  
  // 选择壁面类型（可以通过命令行参数指定）
  std::string wall_type = "step";  // 默认阶梯形
  
  // 处理命令行参数
  if (argc > 1) {
    std::string arg = argv[1];
    if (arg == "-h" || arg == "--help") {
      std::cout << "用法: " << argv[0] << " [壁面类型]" << std::endl;
      std::cout << std::endl;
      std::cout << "壁面类型选项:" << std::endl;
      std::cout << "  step    - 阶梯形壁面（左侧有阶梯，默认）" << std::endl;
      std::cout << "  wave    - 波浪形壁面（右侧是波浪形）" << std::endl;
      std::cout << "  slope   - 倾斜形壁面（底部倾斜）" << std::endl;
      std::cout << "  polygon - 多边形壁面（左侧有五边形障碍物）" << std::endl;
      std::cout << std::endl;
      std::cout << "示例:" << std::endl;
      std::cout << "  " << argv[0] << " step" << std::endl;
      std::cout << "  " << argv[0] << " wave" << std::endl;
      PetscFinalize();
      return 0;
    } else {
      wall_type = arg;
    }
  }
  
  // 验证壁面类型
  std::vector<std::string> valid_types = {"step", "wave", "slope", "polygon"};
  bool valid_type = false;
  for (const auto& type : valid_types) {
    if (wall_type == type) {
      valid_type = true;
      break;
    }
  }
  
  if (!valid_type) {
    std::cerr << "错误：无效的壁面类型 '" << wall_type << "'" << std::endl;
    std::cerr << "支持的壁面类型: step, wave, slope, polygon" << std::endl;
    std::cerr << "使用 -h 或 --help 查看帮助信息" << std::endl;
    PetscFinalize();
    return 1;
  }
  
  std::cout << "========================================" << std::endl;
  std::cout << "不规则几何壁面静水压力测试" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "壁面类型: " << wall_type << std::endl;
  std::cout << "容器尺寸: " << container_width << " x " << container_height << " m" << std::endl;
  std::cout << "水位高度: " << water_height << " m" << std::endl;
  std::cout << "粒子间距: " << particle_spacing << " m" << std::endl;
  
  // 生成测试场景
  FluidParticle fluid_particles("fluid");
  SolidParticle solid_particles("solid");
  GenerateIrregularHydrostaticTest(
      container_width, container_height, water_height, 
      particle_spacing, wall_type,
      fluid_particles, solid_particles);
  
  // 计算参数
  double particle_radius = particle_spacing * 0.5;
  double smoothing_radius = particle_spacing * 2.1;
  double cell_size = particle_spacing * 2.5;
  double time_step = 0.001;  // 时间步长
  double gravity_x = 0.0;
  double gravity_y = -g;
  
  // 构建邻居列表
  std::cout << "\n构建邻居列表..." << std::endl;
  NeighborListSearcher neighbor_searcher;
  neighbor_searcher.BuildNeighborList(
      fluid_particles, solid_particles,
      particle_radius, smoothing_radius, cell_size);
  
  // 检测自由面
  std::cout << "检测自由面..." << std::endl;
  SurfaceDetector surface_detector;
  surface_detector.DetectSurfaceParticles(
      fluid_particles, solid_particles, smoothing_radius, particle_spacing);
  
  int num_surface = 0;
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    if (fluid_particles.surface_type[i] == SurfaceType::SURFACE) {
      ++num_surface;
    }
  }
  std::cout << "  自由面粒子数: " << num_surface << std::endl;
  
  // 计算修正矩阵
  std::cout << "\n计算修正矩阵..." << std::endl;
  CorrectiveMatrix corrective_matrix_calc;
  // 速度散度使用第一类边界条件（border_condition = false）
  std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>> 
      corrective_matrices_velocity(fluid_particles.particle_num);
  // 压力拉普拉斯算子使用第二类边界条件（border_condition = true）
  std::vector<Eigen::Matrix<double, CorrectiveMatrix::MATRIX_SIZE, CorrectiveMatrix::MATRIX_SIZE>> 
      corrective_matrices_pressure(fluid_particles.particle_num);
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    // 速度散度：第一类边界条件
    corrective_matrices_velocity[i] = corrective_matrix_calc.ComputeCorrectiveMatrix(
        i, fluid_particles, solid_particles, smoothing_radius, false);
    // 压力拉普拉斯算子：第二类边界条件
    corrective_matrices_pressure[i] = corrective_matrix_calc.ComputeCorrectiveMatrix(
        i, fluid_particles, solid_particles, smoothing_radius, true);
  }
  std::cout << "  修正矩阵计算完成（速度：第一类边界条件，压力：第二类边界条件）" << std::endl;
  
  // 构建PPE矩阵和右端项
  std::cout << "\n构建PPE矩阵..." << std::endl;
  PPEMatrixBuilder matrix_builder;
  Mat A_petsc = NULL;
  Vec b_petsc = NULL;
  
  if (!matrix_builder.BuildPPEMatrixPetsc(
          fluid_particles, solid_particles, corrective_matrices_velocity, corrective_matrices_pressure,
          smoothing_radius, rho0, time_step, particle_spacing, gravity_x, gravity_y, A_petsc, b_petsc)) {
    std::cerr << "错误：构建PPE矩阵失败" << std::endl;
    return 1;
  }
  
  PetscInt A_m, A_n;
  MatGetSize(A_petsc, &A_m, &A_n);
  MatInfo A_info;
  MatGetInfo(A_petsc, MAT_GLOBAL_SUM, &A_info);
  PetscInt A_nnz = static_cast<PetscInt>(A_info.nz_used);
  std::cout << "  A矩阵大小: " << A_m << " x " << A_n << std::endl;
  std::cout << "  A矩阵非零元素数: " << A_nnz << std::endl;
  
  // 构建罚函数系统
  std::cout << "\n构建罚函数系统..." << std::endl;
  Mat K_petsc = NULL;
  Vec f_petsc = NULL;
  double penalty_parameter = 1e3;
  
  if (!matrix_builder.BuildPenaltySystem(
          A_petsc, b_petsc, fluid_particles, penalty_parameter,
          K_petsc, f_petsc)) {
    std::cerr << "错误：构建罚函数系统失败" << std::endl;
    MatDestroy(&A_petsc);
    VecDestroy(&b_petsc);
    return 1;
  }
  
  PetscInt K_m, K_n;
  MatGetSize(K_petsc, &K_m, &K_n);
  MatInfo K_info;
  MatGetInfo(K_petsc, MAT_GLOBAL_SUM, &K_info);
  PetscInt K_nnz = static_cast<PetscInt>(K_info.nz_used);
  std::cout << "  K矩阵大小: " << K_m << " x " << K_n << std::endl;
  std::cout << "  K矩阵非零元素数: " << K_nnz << std::endl;
  
  // 求解PPE方程
  std::cout << "\n求解PPE方程（罚函数方法，CG求解器）..." << std::endl;
  PPESolver::SolverConfig solver_config;
  solver_config.solver_type = PPESolver::SolverType::CG;
  solver_config.max_iterations = 5000;
  solver_config.tolerance = 1e-6;
  solver_config.force_iterative = true;
  solver_config.is_symmetric_positive_definite = true;
  
  PPESolver ppe_solver(solver_config);
  Vec p_petsc = NULL;
  
  auto start_time = std::chrono::high_resolution_clock::now();
  if (!ppe_solver.Solve(K_petsc, f_petsc, p_petsc)) {
    std::cerr << "错误：PPE求解失败" << std::endl;
    MatDestroy(&A_petsc);
    VecDestroy(&b_petsc);
    MatDestroy(&K_petsc);
    VecDestroy(&f_petsc);
    if (p_petsc != NULL) {
      VecDestroy(&p_petsc);
    }
    return 1;
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
  
  std::cout << "  求解完成" << std::endl;
  std::cout << "  迭代次数: " << ppe_solver.GetLastIterations() << std::endl;
  std::cout << "  残差: " << ppe_solver.GetLastResidual() << std::endl;
  std::cout << "  是否收敛: " << (ppe_solver.GetLastConverged() ? "是" : "否") << std::endl;
  std::cout << "  求解时间: " << duration.count() << " 毫秒" << std::endl;
  
  // 提取压力解
  std::cout << "\n提取压力解..." << std::endl;
  std::vector<PetscInt> indices(fluid_particles.particle_num);
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    indices[i] = static_cast<PetscInt>(i);
  }
  std::vector<PetscScalar> pressure_values(fluid_particles.particle_num);
  VecGetValues(p_petsc, fluid_particles.particle_num, indices.data(), pressure_values.data());
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    fluid_particles.pressure[i] = pressure_values[i];
  }
  
  // 验证自由面粒子压力
  std::cout << "\n验证自由面粒子压力..." << std::endl;
  double max_surface_pressure = 0.0;
  double sum_surface_pressure = 0.0;
  int num_surface_particles = 0;
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    if (fluid_particles.surface_type[i] == SurfaceType::SURFACE) {
      double abs_pressure = std::abs(fluid_particles.pressure[i]);
      max_surface_pressure = std::max(max_surface_pressure, abs_pressure);
      sum_surface_pressure += abs_pressure;
      ++num_surface_particles;
    }
  }
  
  double avg_surface_pressure = num_surface_particles > 0 ? 
      sum_surface_pressure / num_surface_particles : 0.0;
  
  std::cout << "  自由面粒子数: " << num_surface_particles << std::endl;
  std::cout << "  最大压力绝对值: " << max_surface_pressure << " Pa" << std::endl;
  std::cout << "  平均压力绝对值: " << avg_surface_pressure << " Pa" << std::endl;
  
  // 计算理论静水压力并比较
  std::cout << "\n验证静水压力分布..." << std::endl;
  double max_error = 0.0;
  double sum_error = 0.0;
  int num_inner = 0;
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    if (fluid_particles.surface_type[i] == SurfaceType::INNER) {
      double y = fluid_particles.position[i].y;
      double theoretical_pressure = rho0 * g * (water_height - y);
      double error = std::abs(fluid_particles.pressure[i] - theoretical_pressure);
      max_error = std::max(max_error, error);
      sum_error += error;
      ++num_inner;
    }
  }
  
  double avg_error = num_inner > 0 ? sum_error / num_inner : 0.0;
  std::cout << "  内部粒子数: " << num_inner << std::endl;
  std::cout << "  最大误差: " << max_error << " Pa" << std::endl;
  std::cout << "  平均误差: " << avg_error << " Pa" << std::endl;
  
  // 计算压力梯度
  std::cout << "\n计算压力梯度..." << std::endl;
  std::vector<mps2D::double2> pressure_gradient(fluid_particles.particle_num);
  
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    const mps2D::double2& pos_i = fluid_particles.position[i];
    double p_i = pressure_values[i];
    
    // 提取压力corrective matrix的前两行（用于梯度计算，第二类边界条件）
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C1 = corrective_matrices_pressure[i].row(0);  // x方向
    Eigen::RowVector<double, CorrectiveMatrix::BASIS_SIZE> C2 = corrective_matrices_pressure[i].row(1);  // y方向
    
    double grad_x = 0.0;
    double grad_y = 0.0;
    
    // 处理流体邻域粒子
    for (int j : fluid_particles.fluid_neighbour_list[i]) {
      const mps2D::double2& pos_j = fluid_particles.position[j];
      double p_j = pressure_values[j];
      
      double dx = pos_j.x - pos_i.x;
      double dy = pos_j.y - pos_i.y;
      double dist = ComputeDistance(pos_i, pos_j);
      
      if (dist < 1e-10 || dist > smoothing_radius) {
        continue;
      }
      
      double d_ij = (p_j - p_i) / dist;
      double weight = WeightFunction(dist, smoothing_radius);
      
      // 计算基函数（标准基函数）
      Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis = 
          corrective_matrix_calc.ComputeBasisFunctions(dx, dy, smoothing_radius);
      
      // 计算梯度贡献
      grad_x += weight * d_ij * (C1 * basis)(0, 0);
      grad_y += weight * d_ij * (C2 * basis)(0, 0);
    }
    
    // 处理壁面邻域粒子（第二类边界条件）
    // 壁面处压力梯度的法向分量：dp/dn = -rho * g * n_y（对于静水压力）
    for (int j : fluid_particles.solid_neighbour_list[i]) {
      const mps2D::double2& pos_j = solid_particles.position[j];
      const mps2D::double2& normal = solid_particles.normal_vector[j];
      
      double dx = pos_j.x - pos_i.x;
      double dy = pos_j.y - pos_i.y;
      double dist = ComputeDistance(pos_i, pos_j);
      
      if (dist < 1e-10 || dist > smoothing_radius) {
        continue;
      }
      
      double weight = WeightFunction(dist, smoothing_radius);
      
      // 使用壁面基函数（第二类边界条件）
      Eigen::Vector<double, CorrectiveMatrix::BASIS_SIZE> basis_wall = 
          corrective_matrix_calc.ComputeBasisFunctionsForWall(
              dx, dy, normal.x, normal.y, smoothing_radius);
      
      // 壁面处压力梯度的法向分量：dp/dn = -rho * g * n_y（对于静水压力）
      // 归一化法向量
      double nn = std::sqrt(normal.x * normal.x + normal.y * normal.y);
      double n_y = (nn > 1e-10) ? normal.y / nn : 0.0;
      double d_ij = -rho0 * g * n_y;
      
      // 计算梯度贡献
      grad_x += weight * d_ij * (C1 * basis_wall)(0, 0);
      grad_y += weight * d_ij * (C2 * basis_wall)(0, 0);
    }
    
    pressure_gradient[i] = {grad_x, grad_y};
  }
  
  // 计算梯度大小
  std::vector<double> grad_magnitude(fluid_particles.particle_num);
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    grad_magnitude[i] = std::sqrt(pressure_gradient[i].x * pressure_gradient[i].x + 
                                   pressure_gradient[i].y * pressure_gradient[i].y);
  }
  std::cout << "  压力梯度计算完成" << std::endl;
  
  // 输出VTK文件
  std::cout << "\n输出VTK文件..." << std::endl;
  std::string output_filename = "output/ppe_result_irregular_" + wall_type + ".vtk";
  FileOperator file_operator;
  
  // 输出流体粒子基础信息
  if (file_operator.writeVTKBase(output_filename, fluid_particles)) {
    // 追加压力标量
    std::vector<double> pressure_vec(fluid_particles.particle_num);
    for (int i = 0; i < fluid_particles.particle_num; ++i) {
      pressure_vec[i] = pressure_values[i];
    }
    if (file_operator.appendVTKScalar(output_filename, "pressure", pressure_vec)) {
      std::cout << "  压力标量已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加压力标量到VTK文件" << std::endl;
    }
    
    // 追加压力梯度向量
    if (file_operator.appendVTKVector(output_filename, "pressure_gradient", pressure_gradient)) {
      std::cout << "  压力梯度向量已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加压力梯度向量到VTK文件" << std::endl;
    }
    
    // 追加梯度大小标量
    if (file_operator.appendVTKScalar(output_filename, "gradient_magnitude", grad_magnitude)) {
      std::cout << "  梯度大小标量已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加梯度大小标量到VTK文件" << std::endl;
    }
    
    // 追加自由面类型标量
    std::vector<int> surface_type_values(fluid_particles.particle_num);
    for (int i = 0; i < fluid_particles.particle_num; ++i) {
      surface_type_values[i] = static_cast<int>(fluid_particles.surface_type[i]);
    }
    if (file_operator.appendVTKScalar(output_filename, "surface_type", surface_type_values)) {
      std::cout << "  自由面类型标量已追加" << std::endl;
    } else {
      std::cerr << "  警告：无法追加自由面类型标量到VTK文件" << std::endl;
    }
    
    std::cout << "  流体粒子结果已保存到: " << output_filename 
              << " (包含 " << fluid_particles.particle_num << " 个流体粒子)" << std::endl;
  } else {
    std::cerr << "  警告：无法写入VTK文件" << std::endl;
  }
  
  // 单独输出壁面粒子到新的VTK文件
  std::string wall_filename = "output/ppe_result_irregular_" + wall_type + "_wall.vtk";
  std::cout << "\n输出壁面粒子到VTK文件..." << std::endl;
  std::ofstream wall_file(wall_filename);
  if (!wall_file.is_open()) {
    std::cerr << "  警告：无法创建壁面粒子VTK文件" << std::endl;
  } else {
    int num_solid = solid_particles.particle_num;
    
    wall_file << std::fixed << std::setprecision(15);
    
    // VTK文件头
    wall_file << "# vtk DataFile Version 3.0\n";
    wall_file << "MPS Particle Data 2D - Wall Particles (" << wall_type << ")\n";
    wall_file << "ASCII\n";
    wall_file << "DATASET POLYDATA\n";
    
    // 写入点坐标
    wall_file << "POINTS " << num_solid << " float\n";
    for (int i = 0; i < num_solid; ++i) {
      const auto& pos = solid_particles.position[i];
      wall_file << pos.x << " " << pos.y << " 0.0\n";
    }
    
    // 写入顶点
    wall_file << "VERTICES " << num_solid << " " << (num_solid * 2) << "\n";
    for (int i = 0; i < num_solid; ++i) {
      wall_file << "1 " << i << "\n";
    }
    
    // 写入点数据
    wall_file << "POINT_DATA " << num_solid << "\n";
    
    // 写入速度向量
    wall_file << "VECTORS velocity float\n";
    for (int i = 0; i < num_solid; ++i) {
      const auto& vel = solid_particles.velocity[i];
      wall_file << vel.x << " " << vel.y << " 0.0\n";
    }
    
    // 写入法向量
    wall_file << "VECTORS normal_vector float\n";
    for (int i = 0; i < num_solid; ++i) {
      const auto& normal = solid_particles.normal_vector[i];
      wall_file << normal.x << " " << normal.y << " 0.0\n";
    }
    
    // 写入法向量X分量
    wall_file << "SCALARS normal_vector_x float\n";
    wall_file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < num_solid; ++i) {
      wall_file << solid_particles.normal_vector[i].x << "\n";
    }
    
    // 写入法向量Y分量
    wall_file << "SCALARS normal_vector_y float\n";
    wall_file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < num_solid; ++i) {
      wall_file << solid_particles.normal_vector[i].y << "\n";
    }
    
    wall_file.close();
    std::cout << "  壁面粒子结果已保存到: " << wall_filename 
              << " (包含 " << num_solid << " 个壁面粒子)" << std::endl;
  }
  
  // 清理
  MatDestroy(&A_petsc);
  VecDestroy(&b_petsc);
  MatDestroy(&K_petsc);
  VecDestroy(&f_petsc);
  if (p_petsc != NULL) {
    VecDestroy(&p_petsc);
  }
  
  PetscFinalize();
  
  std::cout << "\n测试完成！" << std::endl;
  return 0;
}

