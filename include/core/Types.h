#pragma once
#include<array>
#include<vector>
namespace mps 
{
  using double3 = std::array<double, 3>;
  using int3 = std::array<int, 3>;
  enum class SurfaceType {
    INNER,
    NEAR_SURFACE,
    SURFACE,
    SPLASH
  };
} // namespace mps