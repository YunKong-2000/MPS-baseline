#include "FileOperator.hpp"
#include <algorithm>
#include <cctype>
#include <cstring>
#include <type_traits>

namespace mps2D {

// 从文件名推断文件类型
FileType FileOperator::getFileTypeFromFilename(const std::string& filename) const {
  // 查找最后一个点号的位置
  size_t dot_pos = filename.find_last_of('.');
  
  if (dot_pos == std::string::npos) {
    // 没有扩展名，默认返回TXT
    return FileType::TXT;
  }
  
  // 获取扩展名（转换为小写进行比较）
  std::string extension = filename.substr(dot_pos + 1);
  std::transform(extension.begin(), extension.end(), extension.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  
  if (extension == "csv") {
    return FileType::CSV;
  } else if (extension == "h5" || extension == "hdf5") {
    return FileType::HDF5;
  } else if (extension == "vtk") {
    return FileType::VTK;
  } else {
    // 默认返回TXT（包括.txt和其他未知扩展名）
    return FileType::TXT;
  }
}

// 解析一行数据（2D版本）
bool FileOperator::parseParticleLine(const std::string& line,
                                     double2& position,
                                     double2& velocity,
                                     double2& normal_vector,
                                     bool has_normal) {
  std::istringstream iss(line);
  std::vector<double> values;
  double value;
  
  while (iss >> value) {
    values.push_back(value);
  }
  
  // 2D基本粒子：位置(2) + 速度(2) = 4个值
  // 2D Solid粒子：位置(2) + 速度(2) + 法向向量(2) = 6个值
  size_t required_values = has_normal ? 6 : 4;
  if (values.size() < required_values) {
    return false;
  }
  
  // 读取位置（x, y）
  position.x = values[0];
  position.y = values[1];
  
  // 读取速度（x, y）
  velocity.x = values[2];
  velocity.y = values[3];
  
  // 读取法向向量（如果需要）
  if (has_normal) {
    normal_vector.x = values[4];
    normal_vector.y = values[5];
  }
  
  return true;
}

// TXT文件读取实现
template<typename ParticleType>
int FileOperator::readFromTXT(const std::string& filename, ParticleType& particle) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    return -1;  // 文件打开失败
  }
  
  std::string line;
  int count = 0;
  bool is_solid = std::is_same_v<ParticleType, SolidParticle>;
  
  while (std::getline(file, line)) {
    // 跳过空行和注释行
    if (line.empty() || line[0] == '#') {
      continue;
    }
    
    double2 position = {0.0, 0.0};
    double2 velocity = {0.0, 0.0};
    double2 normal_vector = {0.0, 0.0};
    
    if (parseParticleLine(line, position, velocity, normal_vector, is_solid)) {
      particle.position.push_back(position);
      particle.velocity.push_back(velocity);
      
      if constexpr (std::is_same_v<ParticleType, SolidParticle>) {
        particle.normal_vector.push_back(normal_vector);
      }
      
      count++;
    }
  }
  
  particle.particle_num = count;
  
  // 为FluidParticle调整neighbour_list大小
  if constexpr (std::is_same_v<ParticleType, FluidParticle>) {
    particle.fluid_neighbour_list.resize(count);
    particle.solid_neighbour_list.resize(count);
  }
  
  file.close();
  return count;
}

// CSV文件读取实现
template<typename ParticleType>
int FileOperator::readFromCSV(const std::string& filename, ParticleType& particle) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    return -1;  // 文件打开失败
  }
  
  std::string line;
  int count = 0;
  bool is_solid = std::is_same_v<ParticleType, SolidParticle>;
  
  while (std::getline(file, line)) {
    // 跳过空行和注释行
    if (line.empty() || line[0] == '#') {
      continue;
    }
    
    // CSV格式：用逗号分隔
    std::replace(line.begin(), line.end(), ',', ' ');
    
    double2 position = {0.0, 0.0};
    double2 velocity = {0.0, 0.0};
    double2 normal_vector = {0.0, 0.0};
    
    if (parseParticleLine(line, position, velocity, normal_vector, is_solid)) {
      particle.position.push_back(position);
      particle.velocity.push_back(velocity);
      
      if constexpr (std::is_same_v<ParticleType, SolidParticle>) {
        particle.normal_vector.push_back(normal_vector);
      }
      
      count++;
    }
  }
  
  particle.particle_num = count;
  
  // 为FluidParticle调整neighbour_list大小
  if constexpr (std::is_same_v<ParticleType, FluidParticle>) {
    particle.fluid_neighbour_list.resize(count);
    particle.solid_neighbour_list.resize(count);
  }
  
  file.close();
  return count;
}

// HDF5文件读取实现（暂未实现）
template<typename ParticleType>
int FileOperator::readFromHDF5(const std::string& filename, ParticleType& particle) {
  // TODO: 需要HDF5库支持
  // 当前返回错误
  (void)filename;
  (void)particle;
  return -1;
}

// 通用读取函数
template<typename ParticleType>
int FileOperator::getParticleFromFile(const std::string& filename,
                                      ParticleType& particle) {
  FileType file_type = getFileTypeFromFilename(filename);
  switch (file_type) {
    case FileType::TXT:
      return readFromTXT(filename, particle);
    case FileType::CSV:
      return readFromCSV(filename, particle);
    case FileType::HDF5:
      return readFromHDF5(filename, particle);
    default:
      return -1;
  }
}

// 写入粒子数据到流（辅助函数）
template<typename ParticleType>
void FileOperator::writeParticleData(std::ostream& os,
                                     const ParticleType& particle,
                                     int index,
                                     const std::set<std::string>& custom_fields,
                                     const std::string& delimiter) {
  // 边界检查
  if (index < 0 || index >= static_cast<int>(particle.position.size()) ||
      index >= static_cast<int>(particle.velocity.size())) {
    return;
  }
  
  // 固定写入位置和速度
  const auto& pos = particle.position[index];
  const auto& vel = particle.velocity[index];
  
  os << std::fixed << std::setprecision(15);
  os << pos.x << delimiter << pos.y << delimiter;
  os << vel.x << delimiter << vel.y;
  
  // 写入自定义字段
  if constexpr (std::is_same_v<ParticleType, SolidParticle>) {
    if (custom_fields.find("normal_vector") != custom_fields.end() &&
        index < static_cast<int>(particle.normal_vector.size())) {
      const auto& normal = particle.normal_vector[index];
      os << delimiter << normal.x << delimiter << normal.y;
    }
  } else if constexpr (std::is_same_v<ParticleType, FluidParticle>) {
    if (custom_fields.find("density") != custom_fields.end() && 
        index < static_cast<int>(particle.density.size())) {
      os << delimiter << particle.density[index];
    }
    if (custom_fields.find("pressure") != custom_fields.end() && 
        index < static_cast<int>(particle.pressure.size())) {
      os << delimiter << particle.pressure[index];
    }
    if (custom_fields.find("surface_type") != custom_fields.end() && 
        index < static_cast<int>(particle.surface_type.size())) {
      os << delimiter << static_cast<int>(particle.surface_type[index]);
    }
  }
}

// TXT文件写入实现
template<typename ParticleType>
bool FileOperator::writeToTXT(const std::string& filename,
                              const ParticleType& particle,
                              const std::set<std::string>& custom_fields) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    return false;
  }
  
  file << std::fixed << std::setprecision(15);
  
  // 使用实际的向量大小而不是 particle_num，确保安全
  int actual_count = std::min(particle.particle_num, 
                               static_cast<int>(particle.position.size()));
  for (int i = 0; i < actual_count; ++i) {
    writeParticleData(file, particle, i, custom_fields, " ");
    file << "\n";
  }
  
  file.close();
  return true;
}

// CSV文件写入实现
template<typename ParticleType>
bool FileOperator::writeToCSV(const std::string& filename,
                              const ParticleType& particle,
                              const std::set<std::string>& custom_fields) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    return false;
  }
  
  file << std::fixed << std::setprecision(15);
  
  // 使用实际的向量大小而不是 particle_num，确保安全
  int actual_count = std::min(particle.particle_num, 
                               static_cast<int>(particle.position.size()));
  for (int i = 0; i < actual_count; ++i) {
    writeParticleData(file, particle, i, custom_fields, ",");
    file << "\n";
  }
  
  file.close();
  return true;
}

// HDF5文件写入实现（暂未实现）
template<typename ParticleType>
bool FileOperator::writeToHDF5(const std::string& filename,
                               const ParticleType& particle,
                               const std::set<std::string>& custom_fields) {
  // TODO: 需要HDF5库支持
  (void)filename;
  (void)particle;
  (void)custom_fields;
  return false;
}

// 通用写入函数
template<typename ParticleType>
bool FileOperator::writeParticleToFile(const std::string& filename,
                                      const ParticleType& particle,
                                      const std::set<std::string>& custom_fields) {
  FileType file_type = getFileTypeFromFilename(filename);
  switch (file_type) {
    case FileType::TXT:
      return writeToTXT(filename, particle, custom_fields);
    case FileType::CSV:
      return writeToCSV(filename, particle, custom_fields);
    case FileType::HDF5:
      return writeToHDF5(filename, particle, custom_fields);
    case FileType::VTK:
      // VTK格式不支持通过writeParticleToFile写入，请使用writeVTKBase
      return false;
    default:
      return false;
  }
}

// 基础向量写入操作（用于debug）
template<typename T>
bool FileOperator::writeVectorToFile(const std::string& filename,
                                    const std::vector<T>& data,
                                    const std::string& delimiter) {
  FileType file_type = getFileTypeFromFilename(filename);
  std::ofstream file(filename);
  if (!file.is_open()) {
    return false;
  }
  
  file << std::fixed << std::setprecision(15);
  
  std::string sep = (file_type == FileType::CSV) ? "," : delimiter;
  
  for (size_t i = 0; i < data.size(); ++i) {
    file << data[i];
    if (i < data.size() - 1) {
      file << sep;
    }
  }
  file << "\n";
  
  file.close();
  return true;
}

// 特化：double2向量写入
template<>
bool FileOperator::writeVectorToFile<double2>(const std::string& filename,
                                              const std::vector<double2>& data,
                                              const std::string& delimiter) {
  FileType file_type = getFileTypeFromFilename(filename);
  std::ofstream file(filename);
  if (!file.is_open()) {
    return false;
  }
  
  file << std::fixed << std::setprecision(15);
  
  std::string sep = (file_type == FileType::CSV) ? "," : delimiter;
  
  // 每行写入一个double2向量
  for (size_t i = 0; i < data.size(); ++i) {
    file << data[i].x << sep << data[i].y;
    file << "\n";
  }
  
  file.close();
  return true;
}

// VTK基础输出函数：写入粒子的位置和速度（模板函数）
// 支持类型：FluidParticle, SolidParticle
// 注意：2D数据在VTK中会添加z=0作为第三维
template<typename ParticleType>
bool FileOperator::writeVTKBase(const std::string& filename, const ParticleType& particles) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    return false;
  }
  
  int actual_count = std::min(particles.particle_num, 
                               static_cast<int>(particles.position.size()));
  
  if (actual_count == 0) {
    file.close();
    return false;
  }
  
  file << std::fixed << std::setprecision(15);
  
  // VTK文件头
  file << "# vtk DataFile Version 3.0\n";
  file << "MPS Particle Data 2D - " << particles.name << "\n";
  file << "ASCII\n";
  file << "DATASET POLYDATA\n";
  
  // 写入点坐标（2D数据添加z=0）
  file << "POINTS " << actual_count << " float\n";
  for (int i = 0; i < actual_count; ++i) {
    if (i < static_cast<int>(particles.position.size())) {
      const auto& pos = particles.position[i];
      file << pos.x << " " << pos.y << " 0.0\n";
    }
  }
  
  // 写入顶点（每个点作为一个顶点）
  file << "VERTICES " << actual_count << " " << (actual_count * 2) << "\n";
  for (int i = 0; i < actual_count; ++i) {
    file << "1 " << i << "\n";
  }
  
  // 写入点数据
  file << "POINT_DATA " << actual_count << "\n";
  
  // 写入速度向量（2D数据添加z=0）
  if (actual_count <= static_cast<int>(particles.velocity.size())) {
    file << "VECTORS velocity float\n";
    for (int i = 0; i < actual_count; ++i) {
      const auto& vel = particles.velocity[i];
      file << vel.x << " " << vel.y << " 0.0\n";
    }
  }
  
  file.close();
  return true;
}

// VTK文件解析辅助函数：读取VTK文件内容
bool FileOperator::readVTKFile(const std::string& filename,
                               std::vector<std::string>& header_lines,
                               std::vector<std::string>& point_data_lines,
                               int& num_points) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    return false;
  }
  
  std::string line;
  bool in_point_data = false;
  bool reading_points = false;
  bool reading_vertices = false;
  bool reading_vector_data = false;
  bool reading_scalar_data = false;
  int points_read = 0;
  int vertices_read = 0;
  int vertices_total = 0;
  int vector_data_read = 0;
  int scalar_data_read = 0;
  int expected_scalar_count = 0;
  
  header_lines.clear();
  point_data_lines.clear();
  num_points = 0;
  
  while (std::getline(file, line)) {
    // 解析POINTS行
    if (line.find("POINTS") == 0) {
      std::istringstream iss(line);
      std::string token;
      iss >> token;  // "POINTS"
      iss >> num_points;
      iss >> token;  // "float"
      reading_points = true;
      points_read = 0;
      header_lines.push_back(line);
      continue;
    }
    
    // 解析VERTICES行
    if (line.find("VERTICES") == 0) {
      std::istringstream iss(line);
      std::string token;
      iss >> token;  // "VERTICES"
      int num_vertices;
      iss >> num_vertices;
      iss >> vertices_total;
      reading_vertices = true;
      vertices_read = 0;
      header_lines.push_back(line);
      continue;
    }
    
    // 解析POINT_DATA行
    if (line.find("POINT_DATA") == 0) {
      in_point_data = true;
      header_lines.push_back(line);
      continue;
    }
    
    // 解析VECTORS声明行
    if (in_point_data && line.find("VECTORS") == 0) {
      reading_vector_data = true;
      vector_data_read = 0;
      point_data_lines.push_back(line);
      continue;
    }
    
    // 解析SCALARS声明行
    if (in_point_data && line.find("SCALARS") == 0) {
      reading_scalar_data = true;
      scalar_data_read = 0;
      expected_scalar_count = num_points;
      point_data_lines.push_back(line);
      continue;
    }
    
    // 解析LOOKUP_TABLE行
    if (in_point_data && line.find("LOOKUP_TABLE") == 0) {
      point_data_lines.push_back(line);
      continue;
    }
    
    // 读取点坐标
    if (reading_points && points_read < num_points) {
      header_lines.push_back(line);
      ++points_read;
      if (points_read >= num_points) {
        reading_points = false;
      }
      continue;
    }
    
    // 读取顶点数据
    if (reading_vertices && vertices_read < vertices_total) {
      header_lines.push_back(line);
      ++vertices_read;
      if (vertices_read >= vertices_total) {
        reading_vertices = false;
      }
      continue;
    }
    
    // 读取向量数据
    if (reading_vector_data && vector_data_read < num_points) {
      point_data_lines.push_back(line);
      ++vector_data_read;
      if (vector_data_read >= num_points) {
        reading_vector_data = false;
      }
      continue;
    }
    
    // 读取标量数据
    if (reading_scalar_data && scalar_data_read < expected_scalar_count) {
      point_data_lines.push_back(line);
      ++scalar_data_read;
      if (scalar_data_read >= expected_scalar_count) {
        reading_scalar_data = false;
      }
      continue;
    }
    
    // 其他行
    if (in_point_data) {
      point_data_lines.push_back(line);
    } else {
      header_lines.push_back(line);
    }
  }
  
  file.close();
  return (num_points > 0);
}

// VTK文件快速读取：只读取点数量（用于验证）
bool FileOperator::readVTKPointCount(const std::string& filename, int& num_points) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    return false;
  }
  
  std::string line;
  num_points = 0;
  
  // 只读取POINTS行来获取点数量
  while (std::getline(file, line)) {
    if (line.find("POINTS") == 0) {
      std::istringstream iss(line);
      std::string token;
      iss >> token;  // "POINTS"
      iss >> num_points;
      file.close();
      return (num_points > 0);
    }
  }
  
  file.close();
  return false;
}

// VTK追加函数：向已存在的VTK文件追加标量数据（模板函数）
// 支持类型：int, double
template<typename ScalarType>
bool FileOperator::appendVTKScalar(const std::string& filename,
                                   const std::string& scalar_name,
                                   const std::vector<ScalarType>& scalar_data) {
  // 快速读取点数量进行验证
  int num_points = 0;
  if (!readVTKPointCount(filename, num_points)) {
    return false;
  }
  
  // 检查数据大小是否匹配
  if (static_cast<int>(scalar_data.size()) != num_points) {
    return false;
  }
  
  // 直接以追加模式打开文件
  std::ofstream file(filename, std::ios::app);
  if (!file.is_open()) {
    return false;
  }
  
  file << std::fixed << std::setprecision(15);
  
  // 根据类型确定VTK数据类型
  std::string vtk_type;
  if constexpr (std::is_same_v<ScalarType, int>) {
    vtk_type = "int";
  } else if constexpr (std::is_same_v<ScalarType, double>) {
    vtk_type = "float";
  } else {
    // 默认使用float
    vtk_type = "float";
  }
  
  // 追加新的标量数据
  file << "SCALARS " << scalar_name << " " << vtk_type << "\n";
  file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_points; ++i) {
    file << scalar_data[i] << "\n";
  }
  
  file.close();
  return true;
}

// VTK追加函数：向已存在的VTK文件追加向量数据（模板函数）
// 支持类型：int2, double2（在VTK中会添加z=0作为第三维）
template<typename VectorType>
bool FileOperator::appendVTKVector(const std::string& filename,
                                  const std::string& vector_name,
                                  const std::vector<VectorType>& vector_data) {
  // 快速读取点数量进行验证
  int num_points = 0;
  if (!readVTKPointCount(filename, num_points)) {
    return false;
  }
  
  // 检查数据大小是否匹配
  if (static_cast<int>(vector_data.size()) != num_points) {
    return false;
  }
  
  // 直接以追加模式打开文件
  std::ofstream file(filename, std::ios::app);
  if (!file.is_open()) {
    return false;
  }
  
  file << std::fixed << std::setprecision(15);
  
  // 根据类型确定VTK数据类型
  std::string vtk_type;
  if constexpr (std::is_same_v<VectorType, int2>) {
    vtk_type = "int";
  } else if constexpr (std::is_same_v<VectorType, double2>) {
    vtk_type = "float";
  } else {
    // 默认使用float
    vtk_type = "float";
  }
  
  // 追加新的向量数据（2D向量添加z=0）
  file << "VECTORS " << vector_name << " " << vtk_type << "\n";
  for (int i = 0; i < num_points; ++i) {
    const auto& vec = vector_data[i];
    file << vec.x << " " << vec.y << " 0.0\n";
  }
  
  file.close();
  return true;
}

// 显式模板实例化
template int FileOperator::getParticleFromFile<FluidParticle>(
    const std::string& filename, FluidParticle& particle);
template int FileOperator::getParticleFromFile<SolidParticle>(
    const std::string& filename, SolidParticle& particle);
template int FileOperator::getParticleFromFile<Particle>(
    const std::string& filename, Particle& particle);

template bool FileOperator::writeParticleToFile<FluidParticle>(
    const std::string& filename,
    const FluidParticle& particle,
    const std::set<std::string>& custom_fields);
template bool FileOperator::writeParticleToFile<SolidParticle>(
    const std::string& filename,
    const SolidParticle& particle,
    const std::set<std::string>& custom_fields);
template bool FileOperator::writeParticleToFile<Particle>(
    const std::string& filename,
    const Particle& particle,
    const std::set<std::string>& custom_fields);

// 显式实例化 writeVTKBase
template bool FileOperator::writeVTKBase<FluidParticle>(
    const std::string& filename,
    const FluidParticle& particles);
template bool FileOperator::writeVTKBase<SolidParticle>(
    const std::string& filename,
    const SolidParticle& particles);

// 显式实例化 appendVTKScalar
template bool FileOperator::appendVTKScalar<int>(
    const std::string& filename,
    const std::string& scalar_name,
    const std::vector<int>& scalar_data);
template bool FileOperator::appendVTKScalar<double>(
    const std::string& filename,
    const std::string& scalar_name,
    const std::vector<double>& scalar_data);

// 显式实例化 appendVTKVector
template bool FileOperator::appendVTKVector<int2>(
    const std::string& filename,
    const std::string& vector_name,
    const std::vector<int2>& vector_data);
template bool FileOperator::appendVTKVector<double2>(
    const std::string& filename,
    const std::string& vector_name,
    const std::vector<double2>& vector_data);

// 显式实例化 writeVectorToFile
template bool FileOperator::writeVectorToFile<double>(
    const std::string& filename,
    const std::vector<double>& data,
    const std::string& delimiter);

} // namespace mps2D

