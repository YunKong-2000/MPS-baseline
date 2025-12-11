#include "NeighborListSearcher.hpp"
#include <limits>
#include <algorithm>

namespace mps {

void NeighborListSearcher::BuildNeighborList(FluidParticle& fluid_particles,
                                             SolidParticle& solid_particles,
                                             const double particle_radius,
                                             const double r_e,
                                             const double r_cell) {
  // 步骤1：计算域获取
  std::array<double, 6> domain = ComputeDomain(fluid_particles, solid_particles);
  double3 domain_min = {domain[0], domain[1], domain[2]};
  double3 domain_max = {domain[3], domain[4], domain[5]};

  // 步骤2：背景网格划分
  // 计算各方向的网格数量（添加一个边界网格用于安全）
  int num_cells_x = static_cast<int>(std::ceil((domain_max[0] - domain_min[0]) / r_cell)) + 2;
  int num_cells_y = static_cast<int>(std::ceil((domain_max[1] - domain_min[1]) / r_cell)) + 2;
  int num_cells_z = static_cast<int>(std::ceil((domain_max[2] - domain_min[2]) / r_cell)) + 2;

  // 调整domain_min以包含边界网格
  domain_min[0] -= r_cell;
  domain_min[1] -= r_cell;
  domain_min[2] -= r_cell;

  // 步骤3：背景网格索引获取
  std::vector<int> fluid_cell_indices(fluid_particles.particle_num);
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    int3 grid_coords = GetGridCoordinates(fluid_particles.position[i], domain_min, r_cell);
    fluid_cell_indices[i] = GetCellIndex(grid_coords[0], grid_coords[1], grid_coords[2],
                                         num_cells_x, num_cells_y, num_cells_z);
  }

  std::vector<int> solid_cell_indices(solid_particles.particle_num);
  for (int i = 0; i < solid_particles.particle_num; ++i) {
    int3 grid_coords = GetGridCoordinates(solid_particles.position[i], domain_min, r_cell);
    solid_cell_indices[i] = GetCellIndex(grid_coords[0], grid_coords[1], grid_coords[2],
                                         num_cells_x, num_cells_y, num_cells_z);
  }

  // 步骤4：排序 - 根据网格索引对流体粒子和固体粒子进行排序
  SortFluidParticlesByCellIndex(fluid_particles, fluid_cell_indices);
  SortSolidParticlesByCellIndex(solid_particles, solid_cell_indices);
  
  // 初始化邻居列表
  fluid_particles.fluid_neighbour_list.resize(fluid_particles.particle_num);
  fluid_particles.solid_neighbour_list.resize(fluid_particles.particle_num);
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    fluid_particles.fluid_neighbour_list[i].clear();
    fluid_particles.solid_neighbour_list[i].clear();
  }

  // 步骤5：邻域搜索
  // 为每个网格建立粒子索引映射，提高搜索效率
  int total_cells = num_cells_x * num_cells_y * num_cells_z;
  std::vector<std::vector<int>> fluid_cell_particles(total_cells);
  std::vector<std::vector<int>> solid_cell_particles(total_cells);

  // 建立流体粒子的网格映射（排序后需要重新计算网格索引）
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    int3 grid_coords = GetGridCoordinates(fluid_particles.position[i], domain_min, r_cell);
    int cell_idx = GetCellIndex(grid_coords[0], grid_coords[1], grid_coords[2],
                                num_cells_x, num_cells_y, num_cells_z);
    if (cell_idx >= 0 && cell_idx < total_cells) {
      fluid_cell_particles[cell_idx].push_back(i);
    }
  }

  // 建立固体粒子的网格映射（排序后需要重新计算网格索引）
  for (int i = 0; i < solid_particles.particle_num; ++i) {
    int3 grid_coords = GetGridCoordinates(solid_particles.position[i], domain_min, r_cell);
    int cell_idx = GetCellIndex(grid_coords[0], grid_coords[1], grid_coords[2],
                                num_cells_x, num_cells_y, num_cells_z);
    if (cell_idx >= 0 && cell_idx < total_cells) {
      solid_cell_particles[cell_idx].push_back(i);
    }
  }

  // 对于每个流体粒子，在其所在网格及相邻网格中搜索邻居
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    int3 grid_coords = GetGridCoordinates(fluid_particles.position[i], domain_min, r_cell);
    
    // 在当前网格及其26个相邻网格中搜索
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          int search_grid_x = grid_coords[0] + dx;
          int search_grid_y = grid_coords[1] + dy;
          int search_grid_z = grid_coords[2] + dz;

          // 检查网格坐标是否有效
          if (search_grid_x < 0 || search_grid_x >= num_cells_x ||
              search_grid_y < 0 || search_grid_y >= num_cells_y ||
              search_grid_z < 0 || search_grid_z >= num_cells_z) {
            continue;
          }

          int search_cell_index = GetCellIndex(search_grid_x, search_grid_y, search_grid_z,
                                               num_cells_x, num_cells_y, num_cells_z);

          // 搜索流体粒子邻居（只搜索该网格中的粒子）
          const auto& fluid_in_cell = fluid_cell_particles[search_cell_index];
          for (int j : fluid_in_cell) {
            if (i == j) {
              continue;  // 跳过自己
            }
            double dist = mps::ComputeDistance(fluid_particles.position[i],
                                         fluid_particles.position[j]);
            if (dist < r_e) {
              fluid_particles.fluid_neighbour_list[i].push_back(j);
            }
          }

          // 搜索固体粒子邻居（只搜索该网格中的粒子）
          const auto& solid_in_cell = solid_cell_particles[search_cell_index];
          for (int j : solid_in_cell) {
            double dist = mps::ComputeDistance(fluid_particles.position[i],
                                         solid_particles.position[j]);
            if (dist < r_e) {
              fluid_particles.solid_neighbour_list[i].push_back(j);
            }
          }
        }
      }
    }
  }
}

std::array<double, 6> NeighborListSearcher::ComputeDomain(
    const FluidParticle& fluid_particles,
    const SolidParticle& solid_particles) const {
  double min_x = std::numeric_limits<double>::max();
  double min_y = std::numeric_limits<double>::max();
  double min_z = std::numeric_limits<double>::max();
  double max_x = std::numeric_limits<double>::lowest();
  double max_y = std::numeric_limits<double>::lowest();
  double max_z = std::numeric_limits<double>::lowest();

  // 遍历流体粒子
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    const auto& pos = fluid_particles.position[i];
    min_x = std::min(min_x, pos[0]);
    min_y = std::min(min_y, pos[1]);
    min_z = std::min(min_z, pos[2]);
    max_x = std::max(max_x, pos[0]);
    max_y = std::max(max_y, pos[1]);
    max_z = std::max(max_z, pos[2]);
  }

  // 遍历固体粒子
  for (int i = 0; i < solid_particles.particle_num; ++i) {
    const auto& pos = solid_particles.position[i];
    min_x = std::min(min_x, pos[0]);
    min_y = std::min(min_y, pos[1]);
    min_z = std::min(min_z, pos[2]);
    max_x = std::max(max_x, pos[0]);
    max_y = std::max(max_y, pos[1]);
    max_z = std::max(max_z, pos[2]);
  }

  // 如果没有任何粒子，返回默认域
  if (fluid_particles.particle_num == 0 && solid_particles.particle_num == 0) {
    return {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
  }

  return {min_x, min_y, min_z, max_x, max_y, max_z};
}

int NeighborListSearcher::GetCellIndex(int grid_x, int grid_y, int grid_z,
                                       int num_cells_x, int num_cells_y, int num_cells_z) const {
  return grid_z * num_cells_x * num_cells_y + grid_y * num_cells_x + grid_x;
}

int3 NeighborListSearcher::GetGridCoordinates(const double3& pos,
                                              const double3& domain_min,
                                              double r_cell) const {
  int grid_x = static_cast<int>(std::floor((pos[0] - domain_min[0]) / r_cell));
  int grid_y = static_cast<int>(std::floor((pos[1] - domain_min[1]) / r_cell));
  int grid_z = static_cast<int>(std::floor((pos[2] - domain_min[2]) / r_cell));
  
  // 确保索引非负
  grid_x = std::max(0, grid_x);
  grid_y = std::max(0, grid_y);
  grid_z = std::max(0, grid_z);
  
  return {grid_x, grid_y, grid_z};
}

void NeighborListSearcher::SortFluidParticlesByCellIndex(
    FluidParticle& fluid_particles,
    const std::vector<int>& cell_indices) {
  // 创建索引数组用于排序
  std::vector<int> indices(fluid_particles.particle_num);
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    indices[i] = i;
  }

  // 根据网格索引排序
  std::sort(indices.begin(), indices.end(),
            [&cell_indices](int a, int b) {
              return cell_indices[a] < cell_indices[b];
            });

  // 创建临时存储
  std::vector<double3> temp_position(fluid_particles.particle_num);
  std::vector<double3> temp_velocity(fluid_particles.particle_num);
  std::vector<double> temp_density(fluid_particles.particle_num);
  std::vector<double> temp_pressure(fluid_particles.particle_num);
  std::vector<SurfaceType> temp_surface_type(fluid_particles.particle_num);

  // 按照排序后的顺序重新排列所有属性
  for (int i = 0; i < fluid_particles.particle_num; ++i) {
    int old_idx = indices[i];
    temp_position[i] = fluid_particles.position[old_idx];
    temp_velocity[i] = fluid_particles.velocity[old_idx];
    
    // 安全地拷贝其他属性（如果存在）
    if (old_idx < static_cast<int>(fluid_particles.density.size()) &&
        i < static_cast<int>(temp_density.size())) {
      temp_density[i] = fluid_particles.density[old_idx];
    }
    if (old_idx < static_cast<int>(fluid_particles.pressure.size()) &&
        i < static_cast<int>(temp_pressure.size())) {
      temp_pressure[i] = fluid_particles.pressure[old_idx];
    }
    if (old_idx < static_cast<int>(fluid_particles.surface_type.size()) &&
        i < static_cast<int>(temp_surface_type.size())) {
      temp_surface_type[i] = fluid_particles.surface_type[old_idx];
    }
  }

  // 将排序后的数据复制回原向量
  fluid_particles.position = std::move(temp_position);
  fluid_particles.velocity = std::move(temp_velocity);
  fluid_particles.density = std::move(temp_density);
  fluid_particles.pressure = std::move(temp_pressure);
  fluid_particles.surface_type = std::move(temp_surface_type);
}

void NeighborListSearcher::SortSolidParticlesByCellIndex(
    SolidParticle& solid_particles,
    const std::vector<int>& cell_indices) {
  // 创建索引数组用于排序
  std::vector<int> indices(solid_particles.particle_num);
  for (int i = 0; i < solid_particles.particle_num; ++i) {
    indices[i] = i;
  }

  // 根据网格索引排序
  std::sort(indices.begin(), indices.end(),
            [&cell_indices](int a, int b) {
              return cell_indices[a] < cell_indices[b];
            });

  // 创建临时存储
  std::vector<double3> temp_position(solid_particles.particle_num);
  std::vector<double3> temp_velocity(solid_particles.particle_num);
  std::vector<double3> temp_normal_vector(solid_particles.particle_num);

  // 按照排序后的顺序重新排列所有属性
  for (int i = 0; i < solid_particles.particle_num; ++i) {
    int old_idx = indices[i];
    temp_position[i] = solid_particles.position[old_idx];
    temp_velocity[i] = solid_particles.velocity[old_idx];
    
    // 安全地拷贝法向向量（如果存在）
    if (old_idx < static_cast<int>(solid_particles.normal_vector.size()) &&
        i < static_cast<int>(temp_normal_vector.size())) {
      temp_normal_vector[i] = solid_particles.normal_vector[old_idx];
    }
  }

  // 将排序后的数据复制回原向量
  solid_particles.position = std::move(temp_position);
  solid_particles.velocity = std::move(temp_velocity);
  solid_particles.normal_vector = std::move(temp_normal_vector);
}

} // namespace mps

