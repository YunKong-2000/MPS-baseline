#pragma once
#include<vector>
#include<string>
#include "core/Types.h"
namespace mps 
{
// 前向声明
class FileOperator;

class Particle {
  public:
  std::vector<double3> position;
  std::vector<double3> velocity;
  int particle_num;
  std::string name;

  Particle(std::string name) : name(name) {
    particle_num = 0;
    position.resize(particle_num);
    velocity.resize(particle_num);
  }
  ~Particle() {}
  
  // 从文件读取粒子数据（在基类中实现，派生类共享）
  // 返回读取的粒子数量
  int getParticleFromFile(std::string filename);
  
  // 从同类型粒子中拷贝所有数据（在基类中实现，派生类共享）
  // 当粒子数不一致时，以源粒子数为准
  void copyParticle(const Particle& other);
  
  // 获取粒子数量
  int getParticleNum() const { return particle_num; }
};

class FluidParticle : public Particle {
  public:
  std::vector<double> density;
  std::vector<double> pressure;
  std::vector<SurfaceType> surface_type;
  std::vector<std::vector<int>> fluid_neighbour_list;   // 每个粒子的相邻流体粒子索引列表
  std::vector<std::vector<int>> solid_neighbour_list;   // 每个粒子的相邻固体粒子索引列表
  
  FluidParticle(std::string name) : Particle(name) {
    particle_num = 0;
    density.resize(particle_num);
    pressure.resize(particle_num);
    surface_type.resize(particle_num);
    fluid_neighbour_list.resize(particle_num);
    solid_neighbour_list.resize(particle_num);
  }
  ~FluidParticle() {}
};

class SolidParticle : public Particle {
  public:
  std::vector<double3> normal_vector;
  // 声明FileOperator为友元类，以便访问private成员
  SolidParticle(std::string name) : Particle(name) {
    particle_num = 0;
    normal_vector.resize(particle_num);
  }
  ~SolidParticle() {}
};
} // namespace mps