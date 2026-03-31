// 简单的静水问题粒子生成工具
// 独立编译，不依赖项目其他模块

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

int main(int argc, char* argv[]) {
  std::cout << "=== 生成静水问题算例（约2万粒子）===" << std::endl;
  
  // 几何参数
  double container_width = 1.0;   // 容器宽度 (m)
  double container_height = 2.0;  // 容器高度 (m)
  double water_height = 1.0;      // 水位高度 (m)
  double particle_spacing = 0.02; // 粒子间距 (m)
  
  // 如果通过命令行参数指定粒子间距
  if (argc > 1) {
    particle_spacing = std::stod(argv[1]);
  }
  
  // 计算粒子数
  int nx_fluid = static_cast<int>(std::round(container_width / particle_spacing)) + 1;
  int ny_fluid = static_cast<int>(std::round(water_height / particle_spacing)) + 1;
  int num_fluid = nx_fluid * ny_fluid;
  
  int nx_bottom = static_cast<int>(std::round(container_width / particle_spacing)) + 1;
  int ny_wall = static_cast<int>(std::round(container_height / particle_spacing)) + 1;
  // 壁面当前由：底部一条直线 + 左右两条竖直边构成，为了保证矩形容器四角都有粒子，
  // 下面在统计时会额外加入 4 个角点粒子
  int num_solid = nx_bottom + 2 * ny_wall + 4;
  int total_particles = num_fluid + num_solid;
  
  std::cout << "\n几何参数:" << std::endl;
  std::cout << "  容器宽度: " << container_width << " m" << std::endl;
  std::cout << "  容器高度: " << container_height << " m" << std::endl;
  std::cout << "  水位高度: " << water_height << " m" << std::endl;
  std::cout << "  粒子间距: " << particle_spacing << " m" << std::endl;
  
  std::cout << "\n粒子数:" << std::endl;
  std::cout << "  流体粒子数: " << num_fluid << std::endl;
  std::cout << "  固体粒子数: " << num_solid << std::endl;
  std::cout << "  总粒子数: " << total_particles << std::endl;
  
  // 确保data目录存在
  system("mkdir -p data");
  
  // 输出流体粒子文件
  std::string fluid_file = "data/fluid_particles_2d.txt";
  std::ofstream fluid_out(fluid_file);
  if (!fluid_out.is_open()) {
    std::cerr << "错误：无法打开文件 " << fluid_file << std::endl;
    return 1;
  }
  fluid_out << std::fixed << std::setprecision(15);
  fluid_out << "# Fluid Particle Data File\n";
  fluid_out << "# Format: position_x position_y velocity_x velocity_y\n";
  
  for (int j = 0; j < ny_fluid; ++j) {
    for (int i = 0; i < nx_fluid; ++i) {
      double x = i * particle_spacing;
      double y = j * particle_spacing;
      fluid_out << x << " " << y << " 0.0 0.0\n";
    }
  }
  fluid_out.close();
  std::cout << "\n已输出流体粒子文件: " << fluid_file << std::endl;
  
  // 输出固体粒子文件
  std::string solid_file = "data/solid_particles_2d.txt";
  std::ofstream solid_out(solid_file);
  if (!solid_out.is_open()) {
    std::cerr << "错误：无法打开文件 " << solid_file << std::endl;
    return 1;
  }
  solid_out << std::fixed << std::setprecision(15);
  solid_out << "# Solid Particle Data File\n";
  solid_out << "# Format: position_x position_y velocity_x velocity_y normal_x normal_y\n";
  
  // 底部壁面（法向量向上）
  for (int i = 0; i < nx_bottom; ++i) {
    double x = i * particle_spacing;
    double y = -particle_spacing;
    solid_out << x << " " << y << " 0.0 0.0 0.0 1.0\n";
  }
  
  // 左侧壁面（法向量向右）
  double left_x = -particle_spacing;
  for (int j = 0; j < ny_wall; ++j) {
    double x = left_x;
    double y = j * particle_spacing;
    solid_out << x << " " << y << " 0.0 0.0 1.0 0.0\n";
  }
  
  // 右侧壁面（法向量向左）
  double right_x = container_width + particle_spacing;
  for (int j = 0; j < ny_wall; ++j) {
    double x = right_x;
    double y = j * particle_spacing;
    solid_out << x << " " << y << " 0.0 0.0 -1.0 0.0\n";
  }

  // 容器四角补点：
  // 采用与相邻两条壁面法向量的平均作为角点法向量（再单位化），
  // 使得角点在力学上更接近真实几何角。
  const double inv_sqrt2 = 1.0 / std::sqrt(2.0);

  // 左下角：由底部向上 (0,1) 与左壁向右 (1,0) 组成 -> (1,1)/sqrt(2)
  {
    double x = -particle_spacing;
    double y = -particle_spacing;
    double nx = 1.0 * inv_sqrt2;
    double ny = 1.0 * inv_sqrt2;
    solid_out << x << " " << y << " 0.0 0.0 " << nx << " " << ny << "\n";
  }

  // 右下角：底部向上 (0,1) 与右壁向左 (-1,0) 组成 -> (-1,1)/sqrt(2)
  {
    double x = container_width + particle_spacing;
    double y = -particle_spacing;
    double nx = -1.0 * inv_sqrt2;
    double ny = 1.0 * inv_sqrt2;
    solid_out << x << " " << y << " 0.0 0.0 " << nx << " " << ny << "\n";
  }

  // 左上角：假想顶部壁面向下 (0,-1) 与左壁向右 (1,0) -> (1,-1)/sqrt(2)
  {
    double x = -particle_spacing;
    double y = container_height + particle_spacing;
    double nx = 1.0 * inv_sqrt2;
    double ny = -1.0 * inv_sqrt2;
    solid_out << x << " " << y << " 0.0 0.0 " << nx << " " << ny << "\n";
  }

  // 右上角：假想顶部壁面向下 (0,-1) 与右壁向左 (-1,0) -> (-1,-1)/sqrt(2)
  {
    double x = container_width + particle_spacing;
    double y = container_height + particle_spacing;
    double nx = -1.0 * inv_sqrt2;
    double ny = -1.0 * inv_sqrt2;
    solid_out << x << " " << y << " 0.0 0.0 " << nx << " " << ny << "\n";
  }
  solid_out.close();
  std::cout << "已输出固体粒子文件: " << solid_file << std::endl;
  
  std::cout << "\n算例生成完成！" << std::endl;
  std::cout << "总粒子数: " << total_particles << std::endl;
  
  return 0;
}

