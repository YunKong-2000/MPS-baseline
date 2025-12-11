#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <limits>
#include "../src/surface_detection/SurfaceDetector.hpp"
#include "../src/neighbour_list/NeighborListSearcher.hpp"
#include "../src/core/Particle.hpp"
#include "../src/core/FileOperator.hpp"
#include "core/MPSUtils.h"

using namespace mps;

// 辅助函数：打印分隔线
void printSeparator(const std::string& title) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "  " << title << "\n";
    std::cout << std::string(60, '=') << "\n";
}

// 创建均匀分布的正方体粒子（可以有自由面）
// 参数：
//   cube_size: 正方体边长
//   particle_spacing: 粒子间距
//   origin: 正方体原点位置
//   open_top: 如果为true，顶部开放（有自由面）
FluidParticle createUniformCubeParticles(double cube_size, 
                                        double particle_spacing,
                                        const double3& origin = {0.0, 0.0, 0.0},
                                        bool open_top = true) {
    FluidParticle fluid("cube_fluid");
    
    // 计算每个方向的粒子数量
    int num_x = static_cast<int>(std::ceil(cube_size / particle_spacing)) + 1;
    int num_y = num_x;
    int num_z = open_top ? num_x - 1 : num_x;  // 如果顶部开放，减少一层
    
    // 先计算总粒子数
    int total_particles = num_x * num_y * num_z;
    
    fluid.particle_num = total_particles;
    fluid.position.resize(total_particles);
    fluid.velocity.resize(total_particles);
    fluid.density.resize(total_particles);
    fluid.pressure.resize(total_particles);
    fluid.surface_type.resize(total_particles);
    fluid.fluid_neighbour_list.resize(total_particles);
    fluid.solid_neighbour_list.resize(total_particles);
    
    // 初始化所有粒子
    int idx = 0;
    for (int k = 0; k < num_z; ++k) {
        for (int j = 0; j < num_y; ++j) {
            for (int i = 0; i < num_x; ++i) {
                fluid.position[idx] = {
                    origin[0] + i * particle_spacing,
                    origin[1] + j * particle_spacing,
                    origin[2] + k * particle_spacing
                };
                fluid.velocity[idx] = {0.0, 0.0, 0.0};
                fluid.density[idx] = 1000.0;
                fluid.pressure[idx] = 0.0;
                fluid.surface_type[idx] = SurfaceType::INNER;
                ++idx;
            }
        }
    }
    
    std::cout << "创建了 " << total_particles << " 个均匀分布的粒子\n";
    std::cout << "  正方体尺寸: " << cube_size << " x " << cube_size << " x " 
              << (open_top ? "部分高度" : std::to_string(cube_size)) << "\n";
    std::cout << "  粒子间距: " << particle_spacing << "\n";
    std::cout << "  每个方向粒子数: " << num_x << " x " << num_y << " x " << num_z << "\n";
    if (open_top) {
        std::cout << "  注意: 顶部开放，应该有自由面粒子\n";
    }
    
    return fluid;
}

// 统计自由面粒子类型
void printSurfaceTypeStatistics(const FluidParticle& fluid) {
    int inner_count = 0;
    int near_surface_count = 0;
    int surface_count = 0;
    int splash_count = 0;
    
    for (int i = 0; i < fluid.particle_num; ++i) {
        switch (fluid.surface_type[i]) {
            case SurfaceType::INNER:
                ++inner_count;
                break;
            case SurfaceType::NEAR_SURFACE:
                ++near_surface_count;
                break;
            case SurfaceType::SURFACE:
                ++surface_count;
                break;
            case SurfaceType::SPLASH:
                ++splash_count;
                break;
        }
    }
    
    std::cout << "\n自由面粒子类型统计:\n";
    std::cout << "  内部粒子 (INNER): " << inner_count << "\n";
    std::cout << "  近自由面粒子 (NEAR_SURFACE): " << near_surface_count << "\n";
    std::cout << "  自由面粒子 (SURFACE): " << surface_count << "\n";
    std::cout << "  飞溅粒子 (SPLASH): " << splash_count << "\n";
    std::cout << "  总计: " << fluid.particle_num << "\n";
}

// 将SurfaceType转换为int用于VTK输出
std::vector<int> convertSurfaceTypeToInt(const FluidParticle& fluid) {
    std::vector<int> result(fluid.particle_num);
    for (int i = 0; i < fluid.particle_num; ++i) {
        result[i] = static_cast<int>(fluid.surface_type[i]);
    }
    return result;
}

// 创建矩形容器的壁面（固体粒子）
SolidParticle createRectangularContainer(double length, double width, double height,
                                        double wall_thickness, double particle_spacing,
                                        const double3& origin = {0.0, 0.0, 0.0}) {
    SolidParticle solid("container");
    std::vector<double3> positions;
    
    // 底部
    int num_x_bottom = static_cast<int>(std::ceil(length / particle_spacing)) + 1;
    int num_y_bottom = static_cast<int>(std::ceil(width / particle_spacing)) + 1;
    for (int j = 0; j < num_y_bottom; ++j) {
        for (int i = 0; i < num_x_bottom; ++i) {
            double x = origin[0] + i * particle_spacing;
            double y = origin[1] + j * particle_spacing;
            if (x <= origin[0] + length && y <= origin[1] + width) {
                positions.push_back({x, y, origin[2]});
            }
        }
    }
    
    // 四个侧面
    int num_layers = static_cast<int>(std::ceil(height / particle_spacing)) + 1;
    
    // 前后面 (y = 0 和 y = width)
    for (int k = 0; k < num_layers; ++k) {
        double z = origin[2] + k * particle_spacing;
        if (z > origin[2] + height) break;
        for (int i = 0; i < num_x_bottom; ++i) {
            double x = origin[0] + i * particle_spacing;
            if (x <= origin[0] + length) {
                positions.push_back({x, origin[1], z});  // 前面
                positions.push_back({x, origin[1] + width, z});  // 后面
            }
        }
    }
    
    // 左右面 (x = 0 和 x = length)
    for (int k = 0; k < num_layers; ++k) {
        double z = origin[2] + k * particle_spacing;
        if (z > origin[2] + height) break;
        for (int j = 0; j < num_y_bottom; ++j) {
            double y = origin[1] + j * particle_spacing;
            if (y <= origin[1] + width) {
                positions.push_back({origin[0], y, z});  // 左面
                positions.push_back({origin[0] + length, y, z});  // 右面
            }
        }
    }
    
    solid.particle_num = positions.size();
    solid.position = positions;
    solid.velocity.resize(solid.particle_num);
    solid.normal_vector.resize(solid.particle_num);
    
    for (int i = 0; i < solid.particle_num; ++i) {
        solid.velocity[i] = {0.0, 0.0, 0.0};
        solid.normal_vector[i] = {0.0, 0.0, 1.0};  // 默认法向量向上
    }
    
    std::cout << "创建了 " << solid.particle_num << " 个矩形容器壁面粒子\n";
    std::cout << "  容器尺寸: " << length << " x " << width << " x " << height << "\n";
    
    return solid;
}

// 创建不规则形状容器（L形或U形）
SolidParticle createIrregularContainer(double base_length, double base_width, double height,
                                      double particle_spacing,
                                      const double3& origin = {0.0, 0.0, 0.0}) {
    SolidParticle solid("irregular_container");
    std::vector<double3> positions;
    
    // 创建L形容器
    // 主区域
    double main_length = base_length;
    double main_width = base_width * 0.6;
    // 侧翼区域
    double wing_length = base_length * 0.4;
    double wing_width = base_width * 0.4;
    double wing_offset = base_width * 0.6;
    
    int num_layers = static_cast<int>(std::ceil(height / particle_spacing)) + 1;
    
    // 底部
    for (double y = origin[1]; y <= origin[1] + base_width; y += particle_spacing) {
        for (double x = origin[0]; x <= origin[0] + base_length; x += particle_spacing) {
            // 主区域或侧翼区域
            bool in_main = (x <= origin[0] + main_length && y <= origin[1] + main_width);
            bool in_wing = (x <= origin[0] + wing_length && y >= origin[1] + wing_offset && y <= origin[1] + base_width);
            if (in_main || in_wing) {
                positions.push_back({x, y, origin[2]});
            }
        }
    }
    
    // 侧面壁面
    for (int k = 0; k < num_layers; ++k) {
        double z = origin[2] + k * particle_spacing;
        if (z > origin[2] + height) break;
        
        // 主区域边界
        for (double x = origin[0]; x <= origin[0] + main_length; x += particle_spacing) {
            positions.push_back({x, origin[1], z});  // 前面
            positions.push_back({x, origin[1] + main_width, z});  // 主区域后面
        }
        for (double y = origin[1]; y <= origin[1] + main_width; y += particle_spacing) {
            positions.push_back({origin[0], y, z});  // 左面
            positions.push_back({origin[0] + main_length, y, z});  // 主区域右面
        }
        
        // 侧翼区域边界
        for (double x = origin[0]; x <= origin[0] + wing_length; x += particle_spacing) {
            positions.push_back({x, origin[1] + wing_offset, z});  // 侧翼前面
            positions.push_back({x, origin[1] + base_width, z});  // 侧翼后面
        }
        for (double y = origin[1] + wing_offset; y <= origin[1] + base_width; y += particle_spacing) {
            positions.push_back({origin[0], y, z});  // 侧翼左面
            positions.push_back({origin[0] + wing_length, y, z});  // 侧翼右面
        }
    }
    
    solid.particle_num = positions.size();
    solid.position = positions;
    solid.velocity.resize(solid.particle_num);
    solid.normal_vector.resize(solid.particle_num);
    
    for (int i = 0; i < solid.particle_num; ++i) {
        solid.velocity[i] = {0.0, 0.0, 0.0};
        solid.normal_vector[i] = {0.0, 0.0, 1.0};
    }
    
    std::cout << "创建了 " << solid.particle_num << " 个不规则容器壁面粒子\n";
    std::cout << "  容器类型: L形\n";
    
    return solid;
}

// 创建容器内的液体（流体粒子，确保在容器内部）
FluidParticle createFluidInContainer(const SolidParticle& container,
                                    double fluid_height, double particle_spacing,
                                    double margin = 0.02) {
    FluidParticle fluid("fluid_in_container");
    std::vector<double3> positions;
    
    // 找到容器的边界
    double min_x = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double min_y = std::numeric_limits<double>::max();
    double max_y = std::numeric_limits<double>::lowest();
    double min_z = std::numeric_limits<double>::max();
    double max_z = std::numeric_limits<double>::lowest();
    
    for (int i = 0; i < container.particle_num; ++i) {
        const auto& pos = container.position[i];
        min_x = std::min(min_x, pos[0]);
        max_x = std::max(max_x, pos[0]);
        min_y = std::min(min_y, pos[1]);
        max_y = std::max(max_y, pos[1]);
        min_z = std::min(min_z, pos[2]);
        max_z = std::max(max_z, pos[2]);
    }
    
    double container_bottom = min_z;
    double fluid_top = container_bottom + fluid_height;
    
    // 在容器内生成流体粒子
    int num_x = static_cast<int>(std::ceil((max_x - min_x - 2*margin) / particle_spacing)) + 1;
    int num_y = static_cast<int>(std::ceil((max_y - min_y - 2*margin) / particle_spacing)) + 1;
    int num_z = static_cast<int>(std::ceil(fluid_height / particle_spacing)) + 1;
    
    for (int k = 0; k < num_z; ++k) {
        double z = container_bottom + margin + k * particle_spacing;
        if (z > fluid_top - margin) break;
        
        for (int j = 0; j < num_y; ++j) {
            double y = min_y + margin + j * particle_spacing;
            if (y > max_y - margin) break;
            
            for (int i = 0; i < num_x; ++i) {
                double x = min_x + margin + i * particle_spacing;
                if (x > max_x - margin) break;
                
                // 检查是否在容器内（简单检查，实际应该用更复杂的方法）
                positions.push_back({x, y, z});
            }
        }
    }
    
    fluid.particle_num = positions.size();
    fluid.position = positions;
    fluid.velocity.resize(fluid.particle_num);
    fluid.density.resize(fluid.particle_num);
    fluid.pressure.resize(fluid.particle_num);
    fluid.surface_type.resize(fluid.particle_num);
    fluid.fluid_neighbour_list.resize(fluid.particle_num);
    fluid.solid_neighbour_list.resize(fluid.particle_num);
    
    for (int i = 0; i < fluid.particle_num; ++i) {
        fluid.velocity[i] = {0.0, 0.0, 0.0};
        fluid.density[i] = 1000.0;
        fluid.pressure[i] = 0.0;
        fluid.surface_type[i] = SurfaceType::INNER;
    }
    
    std::cout << "创建了 " << fluid.particle_num << " 个容器内的流体粒子\n";
    std::cout << "  液体高度: " << fluid_height << "\n";
    
    return fluid;
}

// 创建具有大变形自由面的液体（波浪、飞溅等）
FluidParticle createDeformedFluidSurface(const SolidParticle& container,
                                        double base_fluid_height, double particle_spacing,
                                        double wave_amplitude, double wave_frequency,
                                        double splash_height = 0.0,
                                        double margin = 0.02) {
    FluidParticle fluid("deformed_fluid");
    std::vector<double3> positions;
    
    // 找到容器边界
    double min_x = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double min_y = std::numeric_limits<double>::max();
    double max_y = std::numeric_limits<double>::lowest();
    double min_z = std::numeric_limits<double>::max();
    
    for (int i = 0; i < container.particle_num; ++i) {
        const auto& pos = container.position[i];
        min_x = std::min(min_x, pos[0]);
        max_x = std::max(max_x, pos[0]);
        min_y = std::min(min_y, pos[1]);
        max_y = std::max(max_y, pos[1]);
        min_z = std::min(min_z, pos[2]);
    }
    
    double container_bottom = min_z;
    double max_fluid_height = base_fluid_height + wave_amplitude + splash_height;
    
    // 生成流体粒子
    int num_x = static_cast<int>(std::ceil((max_x - min_x - 2*margin) / particle_spacing)) + 1;
    int num_y = static_cast<int>(std::ceil((max_y - min_y - 2*margin) / particle_spacing)) + 1;
    int num_z = static_cast<int>(std::ceil(max_fluid_height / particle_spacing)) + 1;
    
    for (int j = 0; j < num_y; ++j) {
        double y = min_y + margin + j * particle_spacing;
        if (y > max_y - margin) break;
        
        for (int i = 0; i < num_x; ++i) {
            double x = min_x + margin + i * particle_spacing;
            if (x > max_x - margin) break;
            
            // 计算该位置的自由面高度（波浪）
            double wave_z = wave_amplitude * std::sin(wave_frequency * x) * 
                           std::cos(wave_frequency * y);
            double local_surface_height = container_bottom + base_fluid_height + wave_z;
            
            // 生成从底部到自由面的粒子
            int num_z_local = static_cast<int>(std::ceil((local_surface_height - container_bottom) / particle_spacing)) + 1;
            
            for (int k = 0; k < num_z_local; ++k) {
                double z = container_bottom + margin + k * particle_spacing;
                if (z > local_surface_height + margin) break;
                
                // 添加飞溅粒子（在表面附近随机分布）
                if (splash_height > 0 && z > local_surface_height - particle_spacing) {
                    double splash_prob = (z - local_surface_height + particle_spacing) / (particle_spacing * 2);
                    if (splash_prob > 0.3 && splash_prob < 0.7) {
                        // 随机添加一些飞溅粒子
                        double splash_x = x + (std::rand() % 100 - 50) * 0.001;
                        double splash_y = y + (std::rand() % 100 - 50) * 0.001;
                        double splash_z = z + splash_height * (std::rand() % 100) * 0.01;
                        if (splash_z < local_surface_height + splash_height) {
                            positions.push_back({splash_x, splash_y, splash_z});
                        }
                    }
                }
                
                positions.push_back({x, y, z});
            }
        }
    }
    
    fluid.particle_num = positions.size();
    fluid.position = positions;
    fluid.velocity.resize(fluid.particle_num);
    fluid.density.resize(fluid.particle_num);
    fluid.pressure.resize(fluid.particle_num);
    fluid.surface_type.resize(fluid.particle_num);
    fluid.fluid_neighbour_list.resize(fluid.particle_num);
    fluid.solid_neighbour_list.resize(fluid.particle_num);
    
    for (int i = 0; i < fluid.particle_num; ++i) {
        fluid.velocity[i] = {0.0, 0.0, 0.0};
        fluid.density[i] = 1000.0;
        fluid.pressure[i] = 0.0;
        fluid.surface_type[i] = SurfaceType::INNER;
    }
    
    std::cout << "创建了 " << fluid.particle_num << " 个大变形自由面的流体粒子\n";
    std::cout << "  基础高度: " << base_fluid_height << ", 波浪幅度: " << wave_amplitude << "\n";
    
    return fluid;
}

// 创建圆柱形容器中的流体粒子（复杂壁面）
FluidParticle createCylinderFluidParticles(double radius, double height,
                                          double particle_spacing,
                                          const double3& center = {0.0, 0.0, 0.0}) {
    FluidParticle fluid("cylinder_fluid");
    std::vector<double3> positions;
    
    // 在圆柱体内生成粒子
    int num_layers = static_cast<int>(std::ceil(height / particle_spacing)) + 1;
    double z_start = center[2];
    
    for (int layer = 0; layer < num_layers; ++layer) {
        double z = z_start + layer * particle_spacing;
        if (z > center[2] + height) break;
        
        // 在当前层生成圆形分布的粒子
        double r_max = radius - particle_spacing * 0.5;  // 留出边界
        int num_rings = static_cast<int>(std::ceil(r_max / particle_spacing)) + 1;
        
        for (int ring = 0; ring < num_rings; ++ring) {
            double r = ring * particle_spacing;
            if (r > r_max) break;
            
            if (ring == 0) {
                // 中心点
                positions.push_back({center[0], center[1], z});
            } else {
                // 圆周上的点
                int num_points_on_ring = static_cast<int>(std::ceil(2.0 * M_PI * r / particle_spacing));
                if (num_points_on_ring < 1) num_points_on_ring = 1;
                
                for (int i = 0; i < num_points_on_ring; ++i) {
                    double angle = 2.0 * M_PI * i / num_points_on_ring;
                    double x = center[0] + r * std::cos(angle);
                    double y = center[1] + r * std::sin(angle);
                    double dist_from_center = std::sqrt((x - center[0]) * (x - center[0]) + 
                                                       (y - center[1]) * (y - center[1]));
                    if (dist_from_center <= r_max) {
                        positions.push_back({x, y, z});
                    }
                }
            }
        }
    }
    
    fluid.particle_num = positions.size();
    fluid.position = positions;
    fluid.velocity.resize(fluid.particle_num);
    fluid.density.resize(fluid.particle_num);
    fluid.pressure.resize(fluid.particle_num);
    fluid.surface_type.resize(fluid.particle_num);
    fluid.fluid_neighbour_list.resize(fluid.particle_num);
    fluid.solid_neighbour_list.resize(fluid.particle_num);
    
    for (int i = 0; i < fluid.particle_num; ++i) {
        fluid.velocity[i] = {0.0, 0.0, 0.0};
        fluid.density[i] = 1000.0;
        fluid.pressure[i] = 0.0;
        fluid.surface_type[i] = SurfaceType::INNER;
    }
    
    std::cout << "创建了 " << fluid.particle_num << " 个圆柱形容器中的流体粒子\n";
    std::cout << "  圆柱半径: " << radius << ", 高度: " << height << "\n";
    std::cout << "  粒子间距: " << particle_spacing << "\n";
    
    return fluid;
}

// 创建波浪形状的自由面流体（复杂自由面）
FluidParticle createWaveSurfaceParticles(double domain_x, double domain_y, double domain_z,
                                        double particle_spacing,
                                        double wave_amplitude, double wave_frequency,
                                        const double3& origin = {0.0, 0.0, 0.0}) {
    FluidParticle fluid("wave_fluid");
    std::vector<double3> positions;
    
    int num_x = static_cast<int>(std::ceil(domain_x / particle_spacing)) + 1;
    int num_y = static_cast<int>(std::ceil(domain_y / particle_spacing)) + 1;
    int num_z = static_cast<int>(std::ceil(domain_z / particle_spacing)) + 1;
    
    double base_z = origin[2] + domain_z * 0.3;  // 基础高度
    
    for (int k = 0; k < num_z; ++k) {
        double z_base = origin[2] + k * particle_spacing;
        if (z_base > origin[2] + domain_z) break;
        
        for (int j = 0; j < num_y; ++j) {
            double y = origin[1] + j * particle_spacing;
            if (y > origin[1] + domain_y) break;
            
            for (int i = 0; i < num_x; ++i) {
                double x = origin[0] + i * particle_spacing;
                if (x > origin[0] + domain_x) break;
                
                // 计算波浪高度
                double wave_z = wave_amplitude * std::sin(wave_frequency * x) * 
                               std::cos(wave_frequency * y);
                double z = base_z + wave_z + k * particle_spacing;
                
                // 只保留在域内的粒子
                if (z >= origin[2] && z <= origin[2] + domain_z) {
                    positions.push_back({x, y, z});
                }
            }
        }
    }
    
    fluid.particle_num = positions.size();
    fluid.position = positions;
    fluid.velocity.resize(fluid.particle_num);
    fluid.density.resize(fluid.particle_num);
    fluid.pressure.resize(fluid.particle_num);
    fluid.surface_type.resize(fluid.particle_num);
    fluid.fluid_neighbour_list.resize(fluid.particle_num);
    fluid.solid_neighbour_list.resize(fluid.particle_num);
    
    for (int i = 0; i < fluid.particle_num; ++i) {
        fluid.velocity[i] = {0.0, 0.0, 0.0};
        fluid.density[i] = 1000.0;
        fluid.pressure[i] = 0.0;
        fluid.surface_type[i] = SurfaceType::INNER;
    }
    
    std::cout << "创建了 " << fluid.particle_num << " 个波浪形状的流体粒子\n";
    std::cout << "  域尺寸: " << domain_x << " x " << domain_y << " x " << domain_z << "\n";
    std::cout << "  波浪幅度: " << wave_amplitude << ", 频率: " << wave_frequency << "\n";
    
    return fluid;
}

// 创建大规模粒子测试（简单立方体，但粒子数多）
FluidParticle createLargeScaleParticles(double cube_size, double particle_spacing,
                                       const double3& origin = {0.0, 0.0, 0.0}) {
    FluidParticle fluid("large_scale_fluid");
    
    int num_x = static_cast<int>(std::ceil(cube_size / particle_spacing)) + 1;
    int num_y = num_x;
    int num_z = num_x;
    
    int total_particles = num_x * num_y * num_z;
    
    fluid.particle_num = total_particles;
    fluid.position.resize(total_particles);
    fluid.velocity.resize(total_particles);
    fluid.density.resize(total_particles);
    fluid.pressure.resize(total_particles);
    fluid.surface_type.resize(total_particles);
    fluid.fluid_neighbour_list.resize(total_particles);
    fluid.solid_neighbour_list.resize(total_particles);
    
    int idx = 0;
    for (int k = 0; k < num_z; ++k) {
        for (int j = 0; j < num_y; ++j) {
            for (int i = 0; i < num_x; ++i) {
                fluid.position[idx] = {
                    origin[0] + i * particle_spacing,
                    origin[1] + j * particle_spacing,
                    origin[2] + k * particle_spacing
                };
                fluid.velocity[idx] = {0.0, 0.0, 0.0};
                fluid.density[idx] = 1000.0;
                fluid.pressure[idx] = 0.0;
                fluid.surface_type[idx] = SurfaceType::INNER;
                ++idx;
            }
        }
    }
    
    std::cout << "创建了 " << total_particles << " 个大尺度粒子\n";
    std::cout << "  立方体尺寸: " << cube_size << " x " << cube_size << " x " << cube_size << "\n";
    std::cout << "  粒子间距: " << particle_spacing << "\n";
    std::cout << "  每个方向粒子数: " << num_x << " x " << num_y << " x " << num_z << "\n";
    
    return fluid;
}

// 输出流体粒子到VTK文件
bool writeFluidVTK(const std::string& filename,
                   const FluidParticle& fluid,
                   const std::vector<double>* shadow_ratios = nullptr) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    if (fluid.particle_num == 0) {
        file.close();
        return false;
    }
    
    file << std::fixed << std::setprecision(15);
    
    // VTK文件头
    file << "# vtk DataFile Version 3.0\n";
    file << "MPS Particle Data - Fluid Particles\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";
    
    // 写入点坐标
    file << "POINTS " << fluid.particle_num << " float\n";
    for (int i = 0; i < fluid.particle_num; ++i) {
        const auto& pos = fluid.position[i];
        file << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
    }
    
    // 写入顶点（每个点作为一个顶点）
    file << "VERTICES " << fluid.particle_num << " " << (fluid.particle_num * 2) << "\n";
    for (int i = 0; i < fluid.particle_num; ++i) {
        file << "1 " << i << "\n";
    }
    
    // 写入点数据
    file << "POINT_DATA " << fluid.particle_num << "\n";
    
    // 写入速度向量
    file << "VECTORS velocity float\n";
    for (int i = 0; i < fluid.particle_num; ++i) {
        const auto& vel = fluid.velocity[i];
        file << vel[0] << " " << vel[1] << " " << vel[2] << "\n";
    }
    
    // 写入surface_type
    file << "SCALARS surface_type int\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < fluid.particle_num; ++i) {
        file << static_cast<int>(fluid.surface_type[i]) << "\n";
    }
    
    // 写入密度
    if (fluid.density.size() >= static_cast<size_t>(fluid.particle_num)) {
        file << "SCALARS density float\n";
        file << "LOOKUP_TABLE default\n";
        for (int i = 0; i < fluid.particle_num; ++i) {
            file << fluid.density[i] << "\n";
        }
    }
    
    // 写入阴影面积比例（仅在使用虚拟光源法时有效）
    if (shadow_ratios != nullptr && shadow_ratios->size() >= static_cast<size_t>(fluid.particle_num)) {
        file << "SCALARS shadow_area_ratio float\n";
        file << "LOOKUP_TABLE default\n";
        for (int i = 0; i < fluid.particle_num; ++i) {
            file << (*shadow_ratios)[i] << "\n";
        }
    }
    
    file.close();
    return true;
}

// 输出固体粒子到VTK文件
bool writeSolidVTK(const std::string& filename,
                   const SolidParticle& solid) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    if (solid.particle_num == 0) {
        file.close();
        return false;
    }
    
    file << std::fixed << std::setprecision(15);
    
    // VTK文件头
    file << "# vtk DataFile Version 3.0\n";
    file << "MPS Particle Data - Solid Particles\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";
    
    // 写入点坐标
    file << "POINTS " << solid.particle_num << " float\n";
    for (int i = 0; i < solid.particle_num; ++i) {
        const auto& pos = solid.position[i];
        file << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
    }
    
    // 写入顶点（每个点作为一个顶点）
    file << "VERTICES " << solid.particle_num << " " << (solid.particle_num * 2) << "\n";
    for (int i = 0; i < solid.particle_num; ++i) {
        file << "1 " << i << "\n";
    }
    
    // 写入点数据
    file << "POINT_DATA " << solid.particle_num << "\n";
    
    // 写入速度向量
    file << "VECTORS velocity float\n";
    for (int i = 0; i < solid.particle_num; ++i) {
        const auto& vel = solid.velocity[i];
        file << vel[0] << " " << vel[1] << " " << vel[2] << "\n";
    }
    
    file.close();
    return true;
}

// 将流体和固体粒子合并输出到VTK文件（保留用于兼容性）
bool writeCombinedVTK(const std::string& filename,
                      const FluidParticle& fluid,
                      const SolidParticle& solid,
                      const std::vector<double>* shadow_ratios = nullptr) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    int total_particles = fluid.particle_num + solid.particle_num;
    if (total_particles == 0) {
        file.close();
        return false;
    }
    
    file << std::fixed << std::setprecision(15);
    
    // VTK文件头
    file << "# vtk DataFile Version 3.0\n";
    file << "MPS Particle Data - Combined Fluid and Solid\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";
    
    // 写入所有点坐标（先流体，后固体）
    file << "POINTS " << total_particles << " float\n";
    for (int i = 0; i < fluid.particle_num; ++i) {
        const auto& pos = fluid.position[i];
        file << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
    }
    for (int i = 0; i < solid.particle_num; ++i) {
        const auto& pos = solid.position[i];
        file << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
    }
    
    // 写入顶点（每个点作为一个顶点）
    file << "VERTICES " << total_particles << " " << (total_particles * 2) << "\n";
    for (int i = 0; i < total_particles; ++i) {
        file << "1 " << i << "\n";
    }
    
    // 写入点数据
    file << "POINT_DATA " << total_particles << "\n";
    
    // 写入速度向量（先流体，后固体）
    file << "VECTORS velocity float\n";
    for (int i = 0; i < fluid.particle_num; ++i) {
        const auto& vel = fluid.velocity[i];
        file << vel[0] << " " << vel[1] << " " << vel[2] << "\n";
    }
    for (int i = 0; i < solid.particle_num; ++i) {
        const auto& vel = solid.velocity[i];
        file << vel[0] << " " << vel[1] << " " << vel[2] << "\n";
    }
    
    // 写入粒子类型标量（0=流体，1=固体）
    file << "SCALARS particle_type int\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < fluid.particle_num; ++i) {
        file << "0\n";  // 流体粒子
    }
    for (int i = 0; i < solid.particle_num; ++i) {
        file << "1\n";  // 固体粒子
    }
    
    // 写入surface_type（流体粒子有值，固体粒子设为-1表示无效）
    file << "SCALARS surface_type int\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < fluid.particle_num; ++i) {
        file << static_cast<int>(fluid.surface_type[i]) << "\n";
    }
    for (int i = 0; i < solid.particle_num; ++i) {
        file << "-1\n";  // 固体粒子没有surface_type
    }
    
    // 写入密度（流体粒子有值，固体粒子设为0）
    if (fluid.density.size() >= static_cast<size_t>(fluid.particle_num)) {
        file << "SCALARS density float\n";
        file << "LOOKUP_TABLE default\n";
        for (int i = 0; i < fluid.particle_num; ++i) {
            file << fluid.density[i] << "\n";
        }
        for (int i = 0; i < solid.particle_num; ++i) {
            file << "0.0\n";  // 固体粒子密度为0
        }
    }
    
    // 写入阴影面积比例（仅在使用虚拟光源法时有效）
    if (shadow_ratios != nullptr && shadow_ratios->size() >= static_cast<size_t>(fluid.particle_num)) {
        file << "SCALARS shadow_area_ratio float\n";
        file << "LOOKUP_TABLE default\n";
        for (int i = 0; i < fluid.particle_num; ++i) {
            file << (*shadow_ratios)[i] << "\n";
        }
        for (int i = 0; i < solid.particle_num; ++i) {
            file << "-1.0\n";  // 固体粒子没有阴影面积比例
        }
    }
    
    file.close();
    return true;
}

// 执行单个测试案例
void runTestCase(const std::string& test_name,
                FluidParticle& fluid, SolidParticle& solid,
                double particle_radius, double smoothing_radius, double cell_size,
                const std::string& output_filename,
                bool use_virtual_light = false) {
    printSeparator(test_name);
    
    std::cout << "\n测试参数:\n";
    std::cout << "  粒子半径: " << particle_radius << "\n";
    std::cout << "  平滑半径 (r_e): " << smoothing_radius << "\n";
    std::cout << "  背景网格尺寸: " << cell_size << "\n";
    std::cout << "  流体粒子数: " << fluid.particle_num << "\n";
    std::cout << "  固体粒子数: " << solid.particle_num << "\n";
    
    // 构建邻域列表
    printSeparator("构建邻域列表");
    NeighborListSearcher neighbor_searcher;
    neighbor_searcher.BuildNeighborList(fluid, solid, particle_radius, 
                                      smoothing_radius, cell_size);
    
    // 统计邻域粒子数量
    int total_neighbors = 0;
    int max_neighbors = 0;
    int min_neighbors = std::numeric_limits<int>::max();
    for (int i = 0; i < fluid.particle_num; ++i) {
        int neighbor_count = fluid.fluid_neighbour_list[i].size();
        total_neighbors += neighbor_count;
        max_neighbors = std::max(max_neighbors, neighbor_count);
        min_neighbors = std::min(min_neighbors, neighbor_count);
    }
    double avg_neighbors = static_cast<double>(total_neighbors) / fluid.particle_num;
    
    std::cout << "\n邻域列表统计:\n";
    std::cout << "  平均邻域粒子数: " << std::fixed << std::setprecision(2) 
              << avg_neighbors << "\n";
    std::cout << "  最大邻域粒子数: " << max_neighbors << "\n";
    std::cout << "  最小邻域粒子数: " << min_neighbors << "\n";
    
    // 判定自由面粒子
    printSeparator("判定自由面粒子");
    SurfaceDetector surface_detector;
    std::cout << "\n使用的方法: " << (use_virtual_light ? "虚拟光源法" : "锥形区域法") << "\n";
    
    if (use_virtual_light) {
        std::cout << "  虚拟光源法参数:\n";
        std::cout << "    球面网格分辨率: " << 36 << " x " << 72 << "\n";
        std::cout << "    阴影面积阈值（内部粒子）: " << 0.75 << "\n";
        std::cout << "    角落/棱边阴影面积上限: " << 0.60 << "\n";
        std::cout << "    自由面阴影面积上限: " << 0.65 << "\n";
        std::cout << "    最小邻域粒子数（角落/棱边）: " << 8 << "\n";
    }
    
    surface_detector.DetectSurfaceParticles(fluid, solid, smoothing_radius, particle_radius, use_virtual_light);
    
    // 打印统计信息
    printSurfaceTypeStatistics(fluid);
    
    // 如果使用虚拟光源法，计算并存储阴影面积比例
    std::vector<double> shadow_ratios;
    if (use_virtual_light) {
        std::cout << "\n计算阴影面积比例用于VTK输出...\n";
        shadow_ratios.resize(fluid.particle_num);
        for (int i = 0; i < fluid.particle_num; ++i) {
            // 跳过飞溅粒子（它们已经在DetectSurfaceParticles中被处理）
            if (fluid.surface_type[i] == SurfaceType::SPLASH) {
                shadow_ratios[i] = 0.0;  // 飞溅粒子阴影面积为0
                continue;
            }
            
            // 计算阴影面积比例
            int total_neighbors = fluid.fluid_neighbour_list[i].size() + 
                                 fluid.solid_neighbour_list[i].size();
            if (total_neighbors == 0) {
                shadow_ratios[i] = 0.0;  // 没有邻域粒子，阴影面积为0
            } else {
                // 使用SurfaceDetector的公开方法计算阴影面积比例
                shadow_ratios[i] = surface_detector.ComputeShadowAreaRatio(
                    i, fluid, solid, smoothing_radius, particle_radius);
            }
        }
        std::cout << "  已计算 " << fluid.particle_num << " 个粒子的阴影面积比例\n";
    }
    
    // 输出结果到VTK文件（分别输出流体和固体粒子）
    printSeparator("输出结果到VTK文件");
    
    // 生成输出文件名（流体和固体分别）
    std::string fluid_filename = output_filename;
    std::string solid_filename = output_filename;
    
    // 在文件名中插入"_fluid"和"_solid"标识
    size_t dot_pos = fluid_filename.find_last_of(".");
    if (dot_pos != std::string::npos) {
        fluid_filename.insert(dot_pos, "_fluid");
        solid_filename.insert(dot_pos, "_solid");
    } else {
        fluid_filename += "_fluid";
        solid_filename += "_solid";
    }
    
    // 输出流体粒子
    bool success_fluid = writeFluidVTK(fluid_filename, fluid, 
                                      use_virtual_light ? &shadow_ratios : nullptr);
    if (!success_fluid) {
        std::cerr << "错误: 无法写入流体粒子VTK文件: " << fluid_filename << "\n";
        return;
    }
    std::cout << "已写入流体粒子VTK文件: " << fluid_filename << "\n";
    std::cout << "  流体粒子数: " << fluid.particle_num << "\n";
    
    // 输出固体粒子
    bool success_solid = writeSolidVTK(solid_filename, solid);
    if (!success_solid) {
        std::cerr << "错误: 无法写入固体粒子VTK文件: " << solid_filename << "\n";
        return;
    }
    std::cout << "已写入固体粒子VTK文件: " << solid_filename << "\n";
    std::cout << "  固体粒子数: " << solid.particle_num << "\n";
    
    std::cout << "\n测试案例完成！结果已保存到:\n";
    std::cout << "  流体粒子: " << fluid_filename << "\n";
    std::cout << "  固体粒子: " << solid_filename << "\n";
}

int main(int argc, char* argv[]) {
    printSeparator("自由面粒子判定测试套件");
    
    // 解析命令行参数，选择测试案例
    int test_case = 0;
    if (argc > 1) {
        test_case = std::atoi(argv[1]);
    }
    
    // 解析是否使用虚拟光源法（默认使用虚拟光源法）
    bool use_virtual_light = true;
    if (argc > 2) {
        std::string method = argv[2];
        if (method == "cone" || method == "c") {
            use_virtual_light = false;  // 显式指定使用锥形区域法
        }
    }
    
    std::cout << "\n可用测试案例:\n";
    std::cout << "  0 - 基础测试：矩形容器中的液体\n";
    std::cout << "  1 - 复杂壁面：L形不规则容器中的液体\n";
    std::cout << "  2 - 大变形自由面：矩形容器中具有波浪和飞溅的液体\n";
    std::cout << "  3 - 大规模测试：大尺寸容器中的液体\n";
    std::cout << "  4 - 运行所有测试案例\n";
    std::cout << "  5 - 调试测试：小规模正方体（每条边10个粒子），打印虚拟光源法网格标记过程\n";
    std::cout << "\n使用方法: ./test_surface_detection <test_case> [cone|virtual_light]\n";
    std::cout << "  默认使用虚拟光源法，添加 'cone' 或 'c' 使用锥形区域法\n";
    std::cout << "\n当前选择: 测试案例 " << test_case << ", 方法: " 
              << (use_virtual_light ? "虚拟光源法" : "锥形区域法") << "\n";
    
    if (test_case == 0 || test_case == 4) {
        // 测试案例0：基础测试 - 矩形容器中的液体
        double length = 1.0;
        double width = 1.0;
        double height = 1.2;
        double wall_thickness = 0.05;
        double particle_spacing = 0.05;
        double fluid_height = 0.7;
        double particle_radius = particle_spacing * 0.5;
        double smoothing_radius = particle_spacing * 2.5;
        double cell_size = smoothing_radius * 1.5;
        
        SolidParticle container = createRectangularContainer(length, width, height, 
                                                           wall_thickness, particle_spacing);
        FluidParticle fluid = createFluidInContainer(container, fluid_height, particle_spacing);
        
        runTestCase("测试案例0: 基础测试 - 矩形容器中的液体",
                   fluid, container, particle_radius, smoothing_radius, cell_size,
                   "data/test_case0_container.vtk", use_virtual_light);
    }
    
    if (test_case == 1 || test_case == 4) {
        // 测试案例1：复杂壁面 - L形不规则容器
        double base_length = 1.2;
        double base_width = 1.0;
        double height = 1.0;
        double particle_spacing = 0.045;
        double fluid_height = 0.6;
        double particle_radius = particle_spacing * 0.5;
        double smoothing_radius = particle_spacing * 2.5;
        double cell_size = smoothing_radius * 1.5;
        
        SolidParticle container = createIrregularContainer(base_length, base_width, height, 
                                                          particle_spacing);
        FluidParticle fluid = createFluidInContainer(container, fluid_height, particle_spacing);
        
        runTestCase("测试案例1: 复杂壁面 - L形不规则容器中的液体",
                   fluid, container, particle_radius, smoothing_radius, cell_size,
                   "data/test_case1_irregular.vtk", use_virtual_light);
    }
    
    if (test_case == 2 || test_case == 4) {
        // 测试案例2：大变形自由面 - 波浪和飞溅
        double length = 1.5;
        double width = 1.0;
        double height = 1.2;
        double wall_thickness = 0.05;
        double particle_spacing = 0.05;
        double base_fluid_height = 0.6;
        double wave_amplitude = 0.12;
        double wave_frequency = 3.0;
        double splash_height = 0.08;
        double particle_radius = particle_spacing * 0.5;
        double smoothing_radius = particle_spacing * 2.5;
        double cell_size = smoothing_radius * 1.5;
        
        SolidParticle container = createRectangularContainer(length, width, height,
                                                             wall_thickness, particle_spacing);
        FluidParticle fluid = createDeformedFluidSurface(container, base_fluid_height,
                                                         particle_spacing,
                                                         wave_amplitude, wave_frequency,
                                                         splash_height);
        
        runTestCase("测试案例2: 大变形自由面 - 波浪和飞溅液体",
                   fluid, container, particle_radius, smoothing_radius, cell_size,
                   "data/test_case2_deformed.vtk", use_virtual_light);
    }
    
    if (test_case == 3 || test_case == 4) {
        // 测试案例3：大规模测试 - 大尺寸容器
        double length = 2.0;
        double width = 1.5;
        double height = 1.5;
        double wall_thickness = 0.05;
        double particle_spacing = 0.035;  // 更小的间距以获得更多粒子
        double fluid_height = 1.0;
        double particle_radius = particle_spacing * 0.5;
        double smoothing_radius = particle_spacing * 2.5;
        double cell_size = smoothing_radius * 1.5;
        
        SolidParticle container = createRectangularContainer(length, width, height,
                                                             wall_thickness, particle_spacing);
        FluidParticle fluid = createFluidInContainer(container, fluid_height, particle_spacing);
        
        runTestCase("测试案例3: 大规模测试 - 大尺寸容器中的液体",
                   fluid, container, particle_radius, smoothing_radius, cell_size,
                   "data/test_case3_large_scale.vtk", use_virtual_light);
    }
    
    // 调试测试：小规模正方体测例，打印虚拟光源法标记网格过程
    if (test_case == 5 || test_case == 0) {
        printSeparator("调试测试: 小规模正方体 - 虚拟光源法网格标记过程");
        
        // 创建小规模正方体：每条边10个粒子
        double cube_size = 0.9;
        double particle_spacing = cube_size / 9.0;  // 10个粒子，9个间隔
        double particle_radius = particle_spacing * 0.5;
        double smoothing_radius = particle_spacing * 2.5;
        double cell_size = smoothing_radius * 1.5;
        
        // 创建流体粒子（完全填充的正方体）
        FluidParticle fluid = createUniformCubeParticles(cube_size, particle_spacing, 
                                                        {0.0, 0.0, 0.0}, false);
        
        // 创建空的固体粒子（用于测试）
        SolidParticle solid("debug_solid");
        solid.particle_num = 0;
        
        std::cout << "\n测试参数:\n";
        std::cout << "  正方体尺寸: " << cube_size << " x " << cube_size << " x " << cube_size << "\n";
        std::cout << "  粒子间距: " << particle_spacing << "\n";
        std::cout << "  粒子半径: " << particle_radius << "\n";
        std::cout << "  平滑半径: " << smoothing_radius << "\n";
        std::cout << "  流体粒子数: " << fluid.particle_num << "\n";
        
        // 构建邻域列表
        NeighborListSearcher neighbor_searcher;
        neighbor_searcher.BuildNeighborList(fluid, solid, particle_radius, 
                                          smoothing_radius, cell_size);
        
        // 选取位于面上的一个粒子（例如：x=0面上的粒子）
        // 找到x坐标最小（接近0）的粒子
        int surface_particle_idx = -1;
        double min_x = std::numeric_limits<double>::max();
        for (int i = 0; i < fluid.particle_num; ++i) {
            if (fluid.position[i][0] < min_x) {
                min_x = fluid.position[i][0];
                surface_particle_idx = i;
            }
        }
        
        if (surface_particle_idx >= 0) {
            std::cout << "\n选取的面粒子索引: " << surface_particle_idx << "\n";
            std::cout << "粒子位置: (" << std::fixed << std::setprecision(4)
                      << fluid.position[surface_particle_idx][0] << ", "
                      << fluid.position[surface_particle_idx][1] << ", "
                      << fluid.position[surface_particle_idx][2] << ")\n";
            
            // 计算阴影面积比例
            SurfaceDetector surface_detector;
            double shadow_ratio = surface_detector.ComputeShadowAreaRatio(
                surface_particle_idx, fluid, solid, smoothing_radius, particle_radius);
            
            std::cout << "\n最终阴影面积比例: " << std::fixed << std::setprecision(6) 
                      << shadow_ratio << " (" << shadow_ratio * 100.0 << "%)\n";
        } else {
            std::cout << "错误: 未找到合适的面粒子\n";
        }
    }
    
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "所有测试完成！\n";
    std::cout << "\n提示: surface_type值含义:\n";
    std::cout << "  0 = INNER (内部粒子)\n";
    std::cout << "  1 = NEAR_SURFACE (近自由面粒子)\n";
    std::cout << "  2 = SURFACE (自由面粒子)\n";
    std::cout << "  3 = SPLASH (飞溅粒子)\n";
    std::cout << "\n可以使用ParaView查看VTK文件结果\n";
    
    return 0;
}
