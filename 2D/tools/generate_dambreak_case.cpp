// 生成二维溃坝算例的工具程序
// 用途：根据固定几何参数输出流体粒子与壁面粒子文件

#include <cmath>
#include <cerrno>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

namespace {

struct DamBreakConfig {
  double particle_spacing = 0.02;  // 粒子间距 (m)
  double container_width = 3.22;    // 容器长度 (m)
  double container_height = 3.22;   // 容器高度 (m)
  double fluid_width = 1.2;        // 初始液柱宽度 (m)
  double fluid_height = 0.6;       // 初始液柱高度 (m)
  std::string output_dir = "data";
  std::string fluid_file_name = "fluid_particles_dambreak_2d.txt";
  std::string solid_file_name = "solid_particles_dambreak_2d.txt";
};

int CountNodes(double length, double spacing) {
  return static_cast<int>(std::round(length / spacing)) + 1;
}

}  // namespace

int main() {
  const DamBreakConfig config;

  const int nx_fluid = CountNodes(config.fluid_width, config.particle_spacing);
  const int ny_fluid = CountNodes(config.fluid_height, config.particle_spacing);
  const int num_fluid = nx_fluid * ny_fluid;

  const int nx_bottom = CountNodes(config.container_width, config.particle_spacing);
  const int ny_wall = CountNodes(config.container_height, config.particle_spacing);

  // 固体粒子由：底壁 + 左壁 + 右壁 + 四个角点构成
  const int num_solid = nx_bottom + 2 * ny_wall + 4;
  const int total_num = num_fluid + num_solid;

  std::cout << "=== 生成二维溃坝算例 ===" << std::endl;
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "容器尺寸: " << config.container_width << " m x "
            << config.container_height << " m" << std::endl;
  std::cout << "液柱尺寸: " << config.fluid_width << " m x " << config.fluid_height
            << " m (靠左壁)" << std::endl;
  std::cout << "粒子间距: " << config.particle_spacing << " m" << std::endl;
  std::cout << "流体粒子数: " << num_fluid << std::endl;
  std::cout << "壁面粒子数: " << num_solid << std::endl;
  std::cout << "总粒子数: " << total_num << std::endl;

  if (mkdir(config.output_dir.c_str(), 0755) != 0 && errno != EEXIST) {
    std::cerr << "错误：创建目录失败 " << config.output_dir << std::endl;
    return 1;
  }

  const std::string fluid_path = config.output_dir + "/" + config.fluid_file_name;
  std::ofstream fluid_out(fluid_path);
  if (!fluid_out.is_open()) {
    std::cerr << "错误：无法打开文件 " << fluid_path << std::endl;
    return 1;
  }
  fluid_out << std::fixed << std::setprecision(15);
  fluid_out << "# Fluid Particle Data File\n";
  fluid_out << "# Format: position_x position_y velocity_x velocity_y\n";

  // 流体区域：左下角靠近容器左壁，宽 1m、高 1m
  for (int j = 0; j < ny_fluid; ++j) {
    for (int i = 0; i < nx_fluid; ++i) {
      const double x = i * config.particle_spacing;
      const double y = j * config.particle_spacing;
      fluid_out << x << " " << y << " 0.0 0.0\n";
    }
  }
  fluid_out.close();

  const std::string solid_path = config.output_dir + "/" + config.solid_file_name;
  std::ofstream solid_out(solid_path);
  if (!solid_out.is_open()) {
    std::cerr << "错误：无法打开文件 " << solid_path << std::endl;
    return 1;
  }
  solid_out << std::fixed << std::setprecision(15);
  solid_out << "# Solid Particle Data File\n";
  solid_out << "# Format: position_x position_y velocity_x velocity_y normal_x normal_y\n";

  // 底部壁面（法向量向上）
  for (int i = 0; i < nx_bottom; ++i) {
    const double x = i * config.particle_spacing;
    const double y = -config.particle_spacing;
    solid_out << x << " " << y << " 0.0 0.0 0.0 1.0\n";
  }

  // 左侧壁面（法向量向右）
  const double left_x = -config.particle_spacing;
  for (int j = 0; j < ny_wall; ++j) {
    const double y = j * config.particle_spacing;
    solid_out << left_x << " " << y << " 0.0 0.0 1.0 0.0\n";
  }

  // 右侧壁面（法向量向左）
  const double right_x = config.container_width + config.particle_spacing;
  for (int j = 0; j < ny_wall; ++j) {
    const double y = j * config.particle_spacing;
    solid_out << right_x << " " << y << " 0.0 0.0 -1.0 0.0\n";
  }

  // 显式补齐四个角点，避免角部粒子漏洞
  const double inv_sqrt2 = 1.0 / std::sqrt(2.0);

  // 左下角：底壁(0,1) + 左壁(1,0)
  solid_out << -config.particle_spacing << " " << -config.particle_spacing
            << " 0.0 0.0 " << inv_sqrt2 << " " << inv_sqrt2 << "\n";

  // 右下角：底壁(0,1) + 右壁(-1,0)
  solid_out << config.container_width + config.particle_spacing << " "
            << -config.particle_spacing << " 0.0 0.0 "
            << -inv_sqrt2 << " " << inv_sqrt2 << "\n";

  // 左上角：顶部向下(0,-1) + 左壁(1,0)
  solid_out << -config.particle_spacing << " "
            << config.container_height + config.particle_spacing
            << " 0.0 0.0 " << inv_sqrt2 << " " << -inv_sqrt2 << "\n";

  // 右上角：顶部向下(0,-1) + 右壁(-1,0)
  solid_out << config.container_width + config.particle_spacing << " "
            << config.container_height + config.particle_spacing
            << " 0.0 0.0 " << -inv_sqrt2 << " " << -inv_sqrt2 << "\n";

  solid_out.close();

  std::cout << "流体粒子文件: " << fluid_path << std::endl;
  std::cout << "壁面粒子文件: " << solid_path << std::endl;
  std::cout << "算例生成完成。" << std::endl;
  return 0;
}

