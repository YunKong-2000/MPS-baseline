#include "FluidDomainCulling.hpp"

#include <vector>

namespace mps2D {

namespace {

template <typename T>
void CompactByKeepIndices(std::vector<T>* vec, const std::vector<int>& keep) {
  std::vector<T> out;
  out.reserve(keep.size());
  for (int idx : keep) {
    out.push_back((*vec)[static_cast<size_t>(idx)]);
  }
  vec->swap(out);
}

}  // namespace

int RemoveFluidParticlesOutsideDomain(FluidParticle& fluid,
                                      double xmin, double xmax,
                                      double ymin, double ymax) {
  const int n = fluid.particle_num;
  if (n <= 0) {
    return 0;
  }

  std::vector<int> keep;
  keep.reserve(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) {
    const double2& p = fluid.position[static_cast<size_t>(i)];
    if (p.x >= xmin && p.x <= xmax && p.y >= ymin && p.y <= ymax) {
      keep.push_back(i);
    }
  }

  const int removed = n - static_cast<int>(keep.size());
  if (removed == 0) {
    return 0;
  }

  CompactByKeepIndices(&fluid.position, keep);
  CompactByKeepIndices(&fluid.velocity, keep);
  CompactByKeepIndices(&fluid.density, keep);
  CompactByKeepIndices(&fluid.pressure, keep);
  CompactByKeepIndices(&fluid.surface_type, keep);
  CompactByKeepIndices(&fluid.fluid_neighbour_list, keep);
  CompactByKeepIndices(&fluid.solid_neighbour_list, keep);

  fluid.particle_num = static_cast<int>(keep.size());
  return removed;
}

}  // namespace mps2D
