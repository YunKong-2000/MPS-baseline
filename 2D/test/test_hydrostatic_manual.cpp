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
  basis[4] = (normal_x * dx + normal_y * dy) / smoothing_radius;
  return basis;
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
  // 计算系数矩阵
  // 计算流体邻域粒子的系数矩阵
  Eigen::Matrix<double, 5, 5> corrective_matrix = Eigen::Matrix<double, 5, 5>::Zero();
  for (int i = 0; i < fluid_particles.size(); i++) {
    double2 pos_j = fluid_particles[i].position;
    double2 pos_i = fluid_particles[target_fluid_particle_idx].position;
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = std::sqrt(dx * dx + dy * dy);
    if (dist < 1e-10 || dist > smoothing_radius) continue;
    double weight = sqrt(1 - dist / smoothing_radius);
    Eigen::Vector<double, 5> basis = ComputeBasisFunctions(dx, dy, dist, smoothing_radius);
    corrective_matrix += weight * basis * basis.transpose();
  }
  // 计算固体邻域粒子的系数矩阵
  for (int j = 0; j < solid_particles.size(); j++) {
    double2 pos_j = solid_particles[j].position;
    double2 pos_i = fluid_particles[target_fluid_particle_idx].position;
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = std::sqrt(dx * dx + dy * dy);
    if (dist < 1e-10 || dist > smoothing_radius) continue;
    double weight = sqrt(1 - dist / smoothing_radius);
    Eigen::Vector<double, 5> basis = ComputeBasisFunctionsForWall(dx, dy, solid_particles[j].normal_vector.x, solid_particles[j].normal_vector.y, smoothing_radius);
    corrective_matrix += weight * basis * basis.transpose();
  }
  // 求逆
  printMatrix(corrective_matrix, "corrective_matrix");
  Eigen::Matrix<double, 5, 5> inverse_corrective_matrix = corrective_matrix.inverse();
  printMatrix(inverse_corrective_matrix, "inverse_corrective_matrix");
  // 计算流体邻域粒子的压力梯度
  double2 pressure_gradient = {0, 0};
  for (int j = 0; j < fluid_particles.size(); j++) {
    double2 pos_j = fluid_particles[j].position;
    double2 pos_i = fluid_particles[target_fluid_particle_idx].position;
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double dist = std::sqrt(dx * dx + dy * dy);
    if (dist < 1e-10 || dist > smoothing_radius) continue;
    double weight = sqrt(1 - dist / smoothing_radius);
    double p_i = fluid_particles[target_fluid_particle_idx].pressure;
    double p_j = fluid_particles[j].pressure;
    double d_ij = (p_j - p_i) / dist;
    Eigen::Vector<double, 5> basis = ComputeBasisFunctions(dx, dy, dist, smoothing_radius);
    pressure_gradient.x += weight * d_ij * (inverse_corrective_matrix.row(0) * basis)(0, 0);
    pressure_gradient.y += weight * d_ij * (inverse_corrective_matrix.row(1) * basis)(0, 0);
  }
  // 计算固体邻域粒子的压力梯度
  for (int j = 0; j < solid_particles.size(); j++) {
    double2 pos_j = solid_particles[j].position;
    double2 pos_i = fluid_particles[target_fluid_particle_idx].position;
    double dx = pos_j.x - pos_i.x;
    double dy = pos_j.y - pos_i.y;
    double2 n = solid_particles[j].normal_vector;
    double nn = sqrt(n.x * n.x + n.y * n.y);
    double n_y = n.y / nn;
    double dist = std::sqrt(dx * dx + dy * dy);
    if (dist < 1e-10 || dist > smoothing_radius) continue;
    double d_ij = -rho * g * n_y;
    double weight = sqrt(1 - dist / smoothing_radius);
    Eigen::Vector<double, 5> basis = ComputeBasisFunctionsForWall(dx, dy, solid_particles[j].normal_vector.x, solid_particles[j].normal_vector.y, smoothing_radius);
    pressure_gradient.x += weight * d_ij * (inverse_corrective_matrix.row(0) * basis)(0, 0);
    pressure_gradient.y += weight * d_ij * (inverse_corrective_matrix.row(1) * basis)(0, 0);
  }
  std::cout << "pressure gradient: " << pressure_gradient.x << " " << pressure_gradient.y << " Pa/m" << std::endl;
  return 0;
}
