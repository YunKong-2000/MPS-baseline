#pragma once
#include "../core/Particle.hpp"
#include "core/Types.h"
#include "core/MPSUtils.h"
#include <vector>
#include <algorithm>
#include <cmath>

namespace mps2D {

// 邻居列表搜索器类（2D版本）
// 功能：为每个流体粒子寻找周围一定范围内的邻居粒子
class NeighborListSearcher {
public:
  NeighborListSearcher() = default;
  ~NeighborListSearcher() = default;

  // 构建邻居列表
  // 参数：
  //   fluid_particles: 流体粒子对象（会被排序）
  //   solid_particles: 固体粒子对象（会被排序）
  //   particle_radius: 粒子初始间距
  //   r_e: 流体粒子相邻域半径（smoothing_radius）
  //   r_cell: 背景网格尺寸（cell_size）
  void BuildNeighborList(FluidParticle& fluid_particles,
                         SolidParticle& solid_particles,
                         const double particle_radius,
                         const double r_e,
                         const double r_cell);

private:
  // 计算包含所有粒子的二维空间域
  // 返回：{min_x, min_y, max_x, max_y}
  std::array<double, 4> ComputeDomain(const FluidParticle& fluid_particles,
                                      const SolidParticle& solid_particles) const;

  // 将二维网格坐标转换为一维索引
  // 参数：
  //   grid_x, grid_y: 网格坐标
  //   num_cells_x, num_cells_y: 各方向的网格数量
  int GetCellIndex(int grid_x, int grid_y,
                   int num_cells_x, int num_cells_y) const;

  // 获取粒子所在网格的二维坐标
  // 参数：
  //   pos: 粒子位置
  //   domain_min: 域的最小坐标 {min_x, min_y}
  //   r_cell: 网格尺寸
  int2 GetGridCoordinates(const double2& pos,
                         const double2& domain_min,
                         double r_cell) const;

  // 根据网格索引对流体粒子进行排序（包括所有属性）
  void SortFluidParticlesByCellIndex(FluidParticle& fluid_particles,
                                     const std::vector<int>& cell_indices);

  // 根据网格索引对固体粒子进行排序（包括所有属性）
  void SortSolidParticlesByCellIndex(SolidParticle& solid_particles,
                                     const std::vector<int>& cell_indices);
};

} // namespace mps2D

