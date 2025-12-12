#pragma once

namespace mps2D {
  struct double2 {
    double x;
    double y;
  };
  struct int2 {
    int x;
    int y;
  };
  enum class SurfaceType {
    INNER,
    NEAR_SURFACE,
    SURFACE,
    SPLASH
  };
} // namespace mps2D
