#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Particle2D {
  double x;
  double y;
};

struct SolidParticle2D {
  double x;
  double y;
  double nx;
  double ny;
};

int GridCount(double length, double spacing) {
  return static_cast<int>(std::round(length / spacing)) + 1;
}

void WriteFluidFile(const std::string& file_path,
                    const std::vector<Particle2D>& fluid_particles) {
  std::ofstream out(file_path);
  if (!out.is_open()) {
    throw std::runtime_error("无法打开流体粒子文件: " + file_path);
  }

  out << std::fixed << std::setprecision(15);
  out << "# Fluid Particle Data File\n";
  out << "# Format: position_x position_y velocity_x velocity_y\n";
  for (const auto& p : fluid_particles) {
    out << p.x << " " << p.y << " 0.0 0.0\n";
  }
}

void WriteSolidFile(const std::string& file_path,
                    const std::vector<SolidParticle2D>& solid_particles) {
  std::ofstream out(file_path);
  if (!out.is_open()) {
    throw std::runtime_error("无法打开固体粒子文件: " + file_path);
  }

  out << std::fixed << std::setprecision(15);
  out << "# Solid Particle Data File\n";
  out << "# Format: position_x position_y velocity_x velocity_y normal_x normal_y\n";
  for (const auto& p : solid_particles) {
    out << p.x << " " << p.y << " 0.0 0.0 " << p.nx << " " << p.ny << "\n";
  }
}

}  // namespace

int main() {
  constexpr double kParticleSpacing = 0.0025;    // m
  constexpr double kContainerWidth = 1.2;     // m
  constexpr double kContainerHeight = 1.2;    // m
  constexpr double kFluidLength = 0.68;         // m
  constexpr double kFluidHeight = 0.12;         // m
  constexpr double kInvSqrt2 = 0.7071067811865475244;

  const int nx_fluid = GridCount(kFluidLength, kParticleSpacing);
  const int ny_fluid = GridCount(kFluidHeight, kParticleSpacing);
  const int nx_wall = GridCount(kContainerWidth, kParticleSpacing);
  const int ny_wall = GridCount(kContainerHeight, kParticleSpacing);

  std::vector<Particle2D> fluid_particles;
  fluid_particles.reserve(static_cast<size_t>(nx_fluid) * static_cast<size_t>(ny_fluid));
  for (int j = 0; j < ny_fluid; ++j) {
    for (int i = 0; i < nx_fluid; ++i) {
      fluid_particles.push_back({i * kParticleSpacing, j * kParticleSpacing});
    }
  }

  std::vector<SolidParticle2D> solid_particles;
  solid_particles.reserve(static_cast<size_t>(2 * nx_wall + 2 * ny_wall + 4));

  const double left_x = -kParticleSpacing;
  const double right_x = kContainerWidth + kParticleSpacing;
  const double bottom_y = -kParticleSpacing;
  const double top_y = kContainerHeight + kParticleSpacing;

  // 底壁（法向量向上）
  for (int i = 0; i < nx_wall; ++i) {
    solid_particles.push_back({i * kParticleSpacing, bottom_y, 0.0, 1.0});
  }

  // 顶壁（法向量向下）
  for (int i = 0; i < nx_wall; ++i) {
    solid_particles.push_back({i * kParticleSpacing, top_y, 0.0, -1.0});
  }

  // 左壁（法向量向右）
  for (int j = 0; j < ny_wall; ++j) {
    solid_particles.push_back({left_x, j * kParticleSpacing, 1.0, 0.0});
  }

  // 右壁（法向量向左）
  for (int j = 0; j < ny_wall; ++j) {
    solid_particles.push_back({right_x, j * kParticleSpacing, -1.0, 0.0});
  }

  // 四角补点，防止壁面角部不连续
  solid_particles.push_back({left_x, bottom_y, kInvSqrt2, kInvSqrt2});    // 左下角
  solid_particles.push_back({left_x, top_y, kInvSqrt2, -kInvSqrt2});      // 左上角
  solid_particles.push_back({right_x, bottom_y, -kInvSqrt2, kInvSqrt2});  // 右下角
  solid_particles.push_back({right_x, top_y, -kInvSqrt2, -kInvSqrt2});    // 右上角

  std::filesystem::path output_dir;
  const std::filesystem::path cwd = std::filesystem::current_path();
  if (cwd.filename() == "tools") {
    output_dir = cwd / "data";
  } else if (std::filesystem::exists(cwd / "tools")) {
    output_dir = cwd / "tools" / "data";
  } else {
    output_dir = cwd / "data";
  }

  std::filesystem::create_directories(output_dir);
  const std::string fluid_file =
      (output_dir / "fluid_particles_dam_break_Hu.txt").string();
  const std::string solid_file =
      (output_dir / "solid_particles_dam_break_Hu.txt").string();

  try {
    WriteFluidFile(fluid_file, fluid_particles);
    WriteSolidFile(solid_file, solid_particles);
  } catch (const std::exception& e) {
    std::cerr << "错误: " << e.what() << std::endl;
    return 1;
  }

  std::cout << "=== 溃坝算例前处理完成 ===\n";
  std::cout << "粒子间距: " << kParticleSpacing << " m\n";
  std::cout << "容器尺寸: " << kContainerWidth << " m x " << kContainerHeight << " m\n";
  std::cout << "流体初始尺寸: " << kFluidLength << " m x " << kFluidHeight << " m\n";
  std::cout << "流体粒子数: " << fluid_particles.size() << "\n";
  std::cout << "固体粒子数: " << solid_particles.size() << "\n";
  std::cout << "输出文件:\n";
  std::cout << "  " << fluid_file << "\n";
  std::cout << "  " << solid_file << "\n";

  return 0;
}
