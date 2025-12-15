#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
struct double2 {
  double x;
  double y;
};

struct FluidParticle {
  double2 position;
  double pressure;
};

struct SolidParticle {
  double2 position;
  double2 normal_vector;
};

void printMatrix(const Eigen::Matrix<double, 5, 5>& matrix, std::string name) {
  std::cout << name << ":" << std::endl;
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      std::cout << std::setw(10) << matrix(i, j) << " ";
    }
    std::cout << std::endl;
  }
}

Eigen::Vector<double, 5> ComputeBasisFunctions(double dx, double dy, double dist, double smoothing_radius) {
  Eigen::Vector<double, 5> basis;
  basis[0] = dx / dist;
  basis[1] = dy / dist;
  basis[2] = dx * dx / (dist * smoothing_radius);
  basis[3] = dy * dy / (dist * smoothing_radius);
  basis[4] = dx * dy / (dist * smoothing_radius);
  return basis;
}

Eigen::Vector<double, 5> ComputeBasisFunctionsForWall(double dx, double dy, double normal_x, double normal_y, double smoothing_radius) {
  Eigen::Vector<double, 5> basis;
  basis[0] = normal_x;
  basis[1] = normal_y;
  basis[2] = 2.0 * normal_x * dx / smoothing_radius;
  basis[3] = 2.0 * normal_y * dy / smoothing_radius;
  basis[4] = (normal_x * dy + normal_y * dx) / smoothing_radius;
  return basis;
}

// 计算单个粒子的压力梯度
double2 ComputePressureGradient(
    int particle_idx,
    const std::vector<FluidParticle>& fluid_particles,
    const std::vector<SolidParticle>& solid_particles,
    double smoothing_radius,
    double rho,
    double g) {
  
  // 计算系数矩阵
  Eigen::Matrix<double, 5, 5> corrective_matrix = Eigen::Matrix<double, 5, 5>::Zero();
  double2 pos_i = fluid_particles[particle_idx].position;
  
  // 计算流体邻域粒子的系数矩阵
  for (size_t j = 0; j < fluid_particles.size(); j++) {
    double2 pos_j = fluid_particles[j].position;
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = std::sqrt(dx * dx + dy * dy);
    if (dist < 1e-10 || dist > smoothing_radius) continue;
    double weight = pow(1 - dist / smoothing_radius, 2);
    Eigen::Vector<double, 5> basis = ComputeBasisFunctions(dx, dy, dist, smoothing_radius);
    corrective_matrix += weight * basis * basis.transpose();
  }
  
  // 计算固体邻域粒子的系数矩阵
  for (size_t j = 0; j < solid_particles.size(); j++) {
    double2 pos_j = solid_particles[j].position;
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = std::sqrt(dx * dx + dy * dy);
    if (dist < 1e-10 || dist > smoothing_radius) continue;
    double weight = pow(1 - dist / smoothing_radius, 2);
    Eigen::Vector<double, 5> basis = ComputeBasisFunctionsForWall(
        dx, dy, solid_particles[j].normal_vector.x, solid_particles[j].normal_vector.y, smoothing_radius);
    corrective_matrix += weight * basis * basis.transpose();
  }
  
  // 求逆
  Eigen::Matrix<double, 5, 5> inverse_corrective_matrix = corrective_matrix.inverse();
  
  // 计算压力梯度
  double2 pressure_gradient = {0, 0};
  
  // 计算流体邻域粒子的压力梯度贡献
  for (size_t j = 0; j < fluid_particles.size(); j++) {
    double2 pos_j = fluid_particles[j].position;
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = std::sqrt(dx * dx + dy * dy);
    if (dist < 1e-10 || dist > smoothing_radius) continue;
    double weight = pow(1 - dist / smoothing_radius, 2);
    double p_i = fluid_particles[particle_idx].pressure;
    double p_j = fluid_particles[j].pressure;
    double d_ij = (p_j - p_i) / dist;
    Eigen::Vector<double, 5> basis = ComputeBasisFunctions(dx, dy, dist, smoothing_radius);
    pressure_gradient.x += weight * d_ij * (inverse_corrective_matrix.row(0) * basis)(0, 0);
    pressure_gradient.y += weight * d_ij * (inverse_corrective_matrix.row(1) * basis)(0, 0);
  }
  
  // 计算固体邻域粒子的压力梯度贡献
  for (size_t j = 0; j < solid_particles.size(); j++) {
    double2 pos_j = solid_particles[j].position;
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double2 n = solid_particles[j].normal_vector;
    double nn = sqrt(n.x * n.x + n.y * n.y);
    double n_y = n.y / nn;
    double dist = std::sqrt(dx * dx + dy * dy);
    if (dist < 1e-10 || dist > smoothing_radius) continue;
    double d_ij = -rho * g * n_y;
    double weight = pow(1 - dist / smoothing_radius, 2);
    Eigen::Vector<double, 5> basis = ComputeBasisFunctionsForWall(
        dx, dy, solid_particles[j].normal_vector.x, solid_particles[j].normal_vector.y, smoothing_radius);
    pressure_gradient.x += weight * d_ij * (inverse_corrective_matrix.row(0) * basis)(0, 0);
    pressure_gradient.y += weight * d_ij * (inverse_corrective_matrix.row(1) * basis)(0, 0);
  }
  
  return pressure_gradient;
}

// 输出压力梯度到VTK文件
void WritePressureGradientToVTK(
    const std::string& filename,
    const std::vector<FluidParticle>& fluid_particles,
    const std::vector<double2>& pressure_gradients) {
  
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "错误：无法打开文件 " << filename << std::endl;
    return;
  }
  
  file << std::fixed << std::setprecision(15);
  
  // VTK文件头
  file << "# vtk DataFile Version 3.0\n";
  file << "Hydrostatic Pressure Gradient - Manual Test\n";
  file << "ASCII\n";
  file << "DATASET POLYDATA\n";
  
  // 写入点坐标（2D数据添加z=0）
  file << "POINTS " << fluid_particles.size() << " float\n";
  for (size_t i = 0; i < fluid_particles.size(); ++i) {
    const auto& pos = fluid_particles[i].position;
    file << pos.x << " " << pos.y << " 0.0\n";
  }
  
  // 写入顶点（每个点作为一个顶点）
  file << "VERTICES " << fluid_particles.size() << " " << (fluid_particles.size() * 2) << "\n";
  for (size_t i = 0; i < fluid_particles.size(); ++i) {
    file << "1 " << i << "\n";
  }
  
  // 写入点数据
  file << "POINT_DATA " << fluid_particles.size() << "\n";
  
  // 写入压力值
  file << "SCALARS pressure float\n";
  file << "LOOKUP_TABLE default\n";
  for (size_t i = 0; i < fluid_particles.size(); ++i) {
    file << fluid_particles[i].pressure << "\n";
  }
  
  // 写入压力梯度（向量）
  file << "VECTORS pressure_gradient float\n";
  for (size_t i = 0; i < pressure_gradients.size(); ++i) {
    file << pressure_gradients[i].x << " " 
         << pressure_gradients[i].y << " 0.0\n";
  }
  
  // 写入压力梯度X分量
  file << "SCALARS pressure_gradient_x float\n";
  file << "LOOKUP_TABLE default\n";
  for (size_t i = 0; i < pressure_gradients.size(); ++i) {
    file << pressure_gradients[i].x << "\n";
  }
  
  // 写入压力梯度Y分量
  file << "SCALARS pressure_gradient_y float\n";
  file << "LOOKUP_TABLE default\n";
  for (size_t i = 0; i < pressure_gradients.size(); ++i) {
    file << pressure_gradients[i].y << "\n";
  }
  
  // 写入压力梯度大小
  file << "SCALARS pressure_gradient_magnitude float\n";
  file << "LOOKUP_TABLE default\n";
  for (size_t i = 0; i < pressure_gradients.size(); ++i) {
    double mag = std::sqrt(pressure_gradients[i].x * pressure_gradients[i].x + 
                           pressure_gradients[i].y * pressure_gradients[i].y);
    file << mag << "\n";
  }
  
  // 写入Y坐标（用于检查压力分布）
  file << "SCALARS y_coordinate float\n";
  file << "LOOKUP_TABLE default\n";
  for (size_t i = 0; i < fluid_particles.size(); ++i) {
    file << fluid_particles[i].position.y << "\n";
  }
  
  file.close();
  std::cout << "已输出VTK文件: " << filename << std::endl;
}

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <target_fluid_particle_idx>" << std::endl;
    return 1;
  }
  int target_fluid_particle_idx = std::stoi(argv[1]);
  int height = 20;
  int width = 20;
  int depth = 10;
  double particle_spacing = 0.1;
  double smoothing_radius = 2.1 * particle_spacing;
  double rho = 1000;
  double g = 9.8;
  int fluid_particle_num = depth * (width - 2);
  int solid_particle_num = 2 * (height + width - 4);
  std::vector<FluidParticle> fluid_particles;
  std::vector<SolidParticle> solid_particles;
  for (int x = 1; x < width - 1; x++) {
    for (int y = 1; y <= depth; y++) {
        FluidParticle fluid_particle;
        fluid_particle.position = {x * particle_spacing, y * particle_spacing};
        fluid_particle.pressure = rho * g * (depth - y) * particle_spacing;
        fluid_particles.push_back(fluid_particle);
    }
  }
  for (int x = 0; x < width; x++) {
    for (int y = 0; y < height; y++) {
        if (x == 0 || x == width - 1 || y == 0 || y == height - 1) {
            SolidParticle solid_particle;
            solid_particle.position = {x * particle_spacing, y * particle_spacing};
            solid_particle.normal_vector = {0, 0};
            if (x == 0) {
                solid_particle.normal_vector = {1, 0};
            } else if (x == width - 1) {
                solid_particle.normal_vector = {-1, 0};
            } else if (y == 0) {
                solid_particle.normal_vector = {0, 1};
            } else if (y == height - 1) {
                solid_particle.normal_vector = {0, -1};
            }
            solid_particles.push_back(solid_particle);
        }
    }
  }

  std::cout << "target fluid particle idx: " << target_fluid_particle_idx << std::endl;
  std::cout << "target fluid particle position: " << fluid_particles[target_fluid_particle_idx].position.x << " " << fluid_particles[target_fluid_particle_idx].position.y << std::endl;
  std::cout << "target fluid particle pressure: " << fluid_particles[target_fluid_particle_idx].pressure << std::endl;
  
  // 计算目标粒子的压力梯度（用于调试输出）
  double2 target_pressure_gradient = ComputePressureGradient(
      target_fluid_particle_idx, fluid_particles, solid_particles, smoothing_radius, rho, g);
  std::cout << "target particle pressure gradient: " << target_pressure_gradient.x << " " << target_pressure_gradient.y << " Pa/m" << std::endl;
  
  // 计算所有流体粒子的压力梯度
  std::cout << "\n计算所有流体粒子的压力梯度..." << std::endl;
  std::vector<double2> pressure_gradients(fluid_particles.size());
  for (size_t i = 0; i < fluid_particles.size(); ++i) {
    pressure_gradients[i] = ComputePressureGradient(
        i, fluid_particles, solid_particles, smoothing_radius, rho, g);
    if ((i + 1) % 50 == 0 || i == fluid_particles.size() - 1) {
      std::cout << "  已计算 " << (i + 1) << " / " << fluid_particles.size() << " 个粒子" << std::endl;
    }
  }
  
  // 输出到VTK文件
  std::cout << "\n输出压力梯度到VTK文件..." << std::endl;
  WritePressureGradientToVTK("hydrostatic_pressure_gradient_manual.vtk", 
                              fluid_particles, pressure_gradients);
  
  // 统计信息
  std::cout << "\n压力梯度统计信息:" << std::endl;
  double avg_grad_x = 0.0, avg_grad_y = 0.0;
  double max_grad_mag = 0.0, min_grad_mag = 1e10;
  for (size_t i = 0; i < pressure_gradients.size(); ++i) {
    avg_grad_x += pressure_gradients[i].x;
    avg_grad_y += pressure_gradients[i].y;
    double mag = std::sqrt(pressure_gradients[i].x * pressure_gradients[i].x + 
                           pressure_gradients[i].y * pressure_gradients[i].y);
    if (mag > max_grad_mag) max_grad_mag = mag;
    if (mag < min_grad_mag) min_grad_mag = mag;
  }
  avg_grad_x /= fluid_particles.size();
  avg_grad_y /= fluid_particles.size();
  std::cout << "  平均压力梯度: (" << avg_grad_x << ", " << avg_grad_y << ") Pa/m" << std::endl;
  std::cout << "  理论压力梯度: (0.0, -" << (rho * g) << ") Pa/m" << std::endl;
  std::cout << "  压力梯度大小范围: [" << min_grad_mag << ", " << max_grad_mag << "] Pa/m" << std::endl;
  std::cout << "  流体粒子总数: " << fluid_particles.size() << std::endl;
  return 0;
}
