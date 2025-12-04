#include "FileOperator.hpp"
#include <algorithm>
#include <cctype>
#include <cstring>

namespace mps {

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
  } else {
    // 默认返回TXT（包括.txt和其他未知扩展名）
    return FileType::TXT;
  }
}

// 解析一行数据
bool FileOperator::parseParticleLine(const std::string& line,
                                     double3& position,
                                     double3& velocity,
                                     double3& normal_vector,
                                     bool has_normal) {
  std::istringstream iss(line);
  std::vector<double> values;
  double value;
  
  while (iss >> value) {
    values.push_back(value);
  }
  
  // 基本粒子：位置(3) + 速度(3) = 6个值
  // Solid粒子：位置(3) + 速度(3) + 法向向量(3) = 9个值
  int required_values = has_normal ? 9 : 6;
  if (values.size() < required_values) {
    return false;
  }
  
  // 读取位置
  for (int i = 0; i < 3; ++i) {
    position[i] = values[i];
  }
  
  // 读取速度
  for (int i = 0; i < 3; ++i) {
    velocity[i] = values[i + 3];
  }
  
  // 读取法向向量（如果需要）
  if (has_normal) {
    for (int i = 0; i < 3; ++i) {
      normal_vector[i] = values[i + 6];
    }
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
    
    double3 position = {0.0, 0.0, 0.0};
    double3 velocity = {0.0, 0.0, 0.0};
    double3 normal_vector = {0.0, 0.0, 0.0};
    
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
    
    double3 position = {0.0, 0.0, 0.0};
    double3 velocity = {0.0, 0.0, 0.0};
    double3 normal_vector = {0.0, 0.0, 0.0};
    
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
  file.close();
  return count;
}

// HDF5文件读取实现（暂未实现）
template<typename ParticleType>
int FileOperator::readFromHDF5(const std::string& filename, ParticleType& particle) {
  // TODO: 需要HDF5库支持
  // 当前返回错误
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
  os << pos[0] << delimiter << pos[1] << delimiter << pos[2] << delimiter;
  os << vel[0] << delimiter << vel[1] << delimiter << vel[2];
  
  // 写入自定义字段
  if constexpr (std::is_same_v<ParticleType, SolidParticle>) {
    if (custom_fields.find("normal_vector") != custom_fields.end() &&
        index < static_cast<int>(particle.normal_vector.size())) {
      const auto& normal = particle.normal_vector[index];
      os << delimiter << normal[0] << delimiter << normal[1] << delimiter << normal[2];
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

// 特化：double3向量写入
template<>
bool FileOperator::writeVectorToFile<double3>(const std::string& filename,
                                              const std::vector<double3>& data,
                                              const std::string& delimiter) {
  FileType file_type = getFileTypeFromFilename(filename);
  std::ofstream file(filename);
  if (!file.is_open()) {
    return false;
  }
  
  file << std::fixed << std::setprecision(15);
  
  std::string sep = (file_type == FileType::CSV) ? "," : delimiter;
  
  // 每行写入一个double3向量
  for (size_t i = 0; i < data.size(); ++i) {
    file << data[i][0] << sep << data[i][1] << sep << data[i][2];
    file << "\n";
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

// 显式实例化 writeVectorToFile
template bool FileOperator::writeVectorToFile<double>(
    const std::string& filename,
    const std::vector<double>& data,
    const std::string& delimiter);

} // namespace mps
