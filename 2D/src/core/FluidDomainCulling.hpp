#pragma once

#include "Particle.hpp"

namespace mps2D {

// 删除位置落在闭域 [xmin,xmax]×[ymin,ymax] 之外的流体粒子。
// 邻域列表会在下一时间步开始重建；本函数与并行兼容（纯数据压缩，后续可改为并行筛选）。
// 返回值：被删除的粒子数量。
int RemoveFluidParticlesOutsideDomain(FluidParticle& fluid,
                                      double xmin, double xmax,
                                      double ymin, double ymax);

}  // namespace mps2D
