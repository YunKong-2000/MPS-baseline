#include "Particle.hpp"
#include "FileOperator.hpp"
#include <typeinfo>

namespace mps {

int Particle::getParticleFromFile(std::string filename) {
  FileOperator file_op;
  
  // 使用运行时类型信息（RTTI）来确定具体类型
  // 然后调用相应的模板实例化
  const std::type_info& type = typeid(*this);
  
  if (type == typeid(SolidParticle)) {
    // 转换为 SolidParticle 并调用模板函数
    SolidParticle& solid = static_cast<SolidParticle&>(*this);
    return file_op.getParticleFromFile(filename, solid);
  } else if (type == typeid(FluidParticle)) {
    // 转换为 FluidParticle 并调用模板函数
    FluidParticle& fluid = static_cast<FluidParticle&>(*this);
    return file_op.getParticleFromFile(filename, fluid);
  } else {
    // 对于基类 Particle，直接调用（不处理派生类特有数据）
    return file_op.getParticleFromFile(filename, *this);
  }
}

void Particle::copyParticle(const Particle& other) {
  // 检查类型是否匹配
  const std::type_info& this_type = typeid(*this);
  const std::type_info& other_type = typeid(other);
  
  if (this_type != other_type) {
    // 类型不匹配，无法拷贝
    return;
  }
  
  // 以源粒子数为准，调整目标粒子大小
  particle_num = other.particle_num;
  
  // 拷贝基类数据
  position = other.position;
  velocity = other.velocity;
  name = other.name;
  
  // 根据具体类型拷贝派生类特有数据
  if (this_type == typeid(SolidParticle)) {
    const SolidParticle& other_solid = static_cast<const SolidParticle&>(other);
    SolidParticle& this_solid = static_cast<SolidParticle&>(*this);
    
    // 拷贝法向向量（赋值操作符会自动调整大小以匹配源向量）
    this_solid.normal_vector = other_solid.normal_vector;
  } else if (this_type == typeid(FluidParticle)) {
    const FluidParticle& other_fluid = static_cast<const FluidParticle&>(other);
    FluidParticle& this_fluid = static_cast<FluidParticle&>(*this);
    
    // 拷贝密度、压力和表面类型（赋值操作符会自动调整大小以匹配源向量）
    this_fluid.density = other_fluid.density;
    this_fluid.pressure = other_fluid.pressure;
    this_fluid.surface_type = other_fluid.surface_type;
  }
  // 如果是基类 Particle，只拷贝基类数据即可
}

} // namespace mps