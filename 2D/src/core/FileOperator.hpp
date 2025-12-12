#pragma once
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <set>
#include <iomanip>
#include "Particle.hpp"
#include "core/Types.h"

namespace mps2D {

enum class FileType {
  TXT,
  CSV,
  HDF5,
  VTK,
};

class FileOperator {
public:
  FileOperator() = default;
  ~FileOperator() = default;

  // 粒子数据读取操作
  // 参数：文件名、粒子引用
  // 返回：读取的粒子数量（文件类型从文件名自动推断）
  template<typename ParticleType>
  int getParticleFromFile(const std::string& filename, 
                          ParticleType& particle);

  // 粒子数据写入操作
  // 参数：文件名、粒子引用、自定义字段集合（可选）
  // 固定写入位置和速度，custom_fields指定其他要写入的字段（文件类型从文件名自动推断）
  template<typename ParticleType>
  bool writeParticleToFile(const std::string& filename,
                           const ParticleType& particle,
                           const std::set<std::string>& custom_fields = {});

  // 基础向量写入操作（用于debug）
  // 文件类型从文件名自动推断
  template<typename T>
  bool writeVectorToFile(const std::string& filename,
                        const std::vector<T>& data,
                        const std::string& delimiter = " ");

  // VTK基础输出函数：写入位置和速度信息（模板函数）
  // 支持类型：FluidParticle, SolidParticle
  // 注意：2D数据在VTK中会添加z=0作为第三维
  template<typename ParticleType>
  bool writeVTKBase(const std::string& filename, const ParticleType& particles);

  // VTK追加函数：向已存在的VTK文件追加标量数据（模板函数）
  // 支持类型：int, double
  template<typename ScalarType>
  bool appendVTKScalar(const std::string& filename,
                       const std::string& scalar_name,
                       const std::vector<ScalarType>& scalar_data);

  // VTK追加函数：向已存在的VTK文件追加向量数据（模板函数）
  // 支持类型：int2, double2（在VTK中会添加z=0作为第三维）
  template<typename VectorType>
  bool appendVTKVector(const std::string& filename,
                       const std::string& vector_name,
                       const std::vector<VectorType>& vector_data);

private:
  // TXT文件读取实现
  template<typename ParticleType>
  int readFromTXT(const std::string& filename, ParticleType& particle);

  // CSV文件读取实现
  template<typename ParticleType>
  int readFromCSV(const std::string& filename, ParticleType& particle);

  // HDF5文件读取实现（暂未实现，需要HDF5库支持）
  template<typename ParticleType>
  int readFromHDF5(const std::string& filename, ParticleType& particle);

  // TXT文件写入实现
  template<typename ParticleType>
  bool writeToTXT(const std::string& filename,
                  const ParticleType& particle,
                  const std::set<std::string>& custom_fields);

  // CSV文件写入实现
  template<typename ParticleType>
  bool writeToCSV(const std::string& filename,
                  const ParticleType& particle,
                  const std::set<std::string>& custom_fields);

  // HDF5文件写入实现（暂未实现，需要HDF5库支持）
  template<typename ParticleType>
  bool writeToHDF5(const std::string& filename,
                   const ParticleType& particle,
                   const std::set<std::string>& custom_fields);

  // 从文件名推断文件类型
  FileType getFileTypeFromFilename(const std::string& filename) const;

  // 解析一行数据（位置、速度，可能还有法向向量）
  // 2D版本：基本粒子需要4个值（位置x,y + 速度x,y）
  // SolidParticle需要6个值（位置x,y + 速度x,y + 法向向量x,y）
  bool parseParticleLine(const std::string& line,
                        double2& position,
                        double2& velocity,
                        double2& normal_vector,
                        bool has_normal);

  // 写入粒子数据到流（辅助函数）
  template<typename ParticleType>
  void writeParticleData(std::ostream& os,
                         const ParticleType& particle,
                         int index,
                         const std::set<std::string>& custom_fields,
                         const std::string& delimiter);

  // VTK文件解析辅助函数：读取VTK文件内容
  bool readVTKFile(const std::string& filename,
                   std::vector<std::string>& header_lines,
                   std::vector<std::string>& point_data_lines,
                   int& num_points);

  // VTK文件快速读取：只读取点数量（用于验证）
  bool readVTKPointCount(const std::string& filename, int& num_points);
};

// double2向量写入特化声明
template<>
bool FileOperator::writeVectorToFile<double2>(const std::string& filename,
                                              const std::vector<mps2D::double2>& data,
                                              const std::string& delimiter);

} // namespace mps2D
