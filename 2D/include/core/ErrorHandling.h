#pragma once

#include <stdexcept>
#include <string>
#include <fmt/format.h>

namespace mps2D {

/**
 * @brief 自定义运行时异常类，用于错误处理
 */
class MPSException : public std::runtime_error {
 public:
  explicit MPSException(const std::string& message)
      : std::runtime_error(message) {}
  
  explicit MPSException(const char* message)
      : std::runtime_error(message) {}
};

/**
 * @brief 检查条件，如果失败则抛出异常
 * @param condition 要检查的条件
 * @param message 错误消息（支持fmt格式化）
 * @param args fmt格式化参数
 */
template<typename... Args>
inline void CHECK(bool condition, const std::string& message, Args&&... args) {
  if (!condition) {
    std::string formatted_msg = fmt::format(message, std::forward<Args>(args)...);
    throw MPSException(formatted_msg);
  }
}

/**
 * @brief 检查条件，如果失败则抛出异常（无格式化版本）
 * @param condition 要检查的条件
 * @param message 错误消息
 */
inline void CHECK(bool condition, const std::string& message) {
  if (!condition) {
    throw MPSException(message);
  }
}

/**
 * @brief 检查指针非空
 * @param ptr 要检查的指针
 * @param message 错误消息（支持fmt格式化）
 * @param args fmt格式化参数
 */
template<typename T, typename... Args>
inline void CHECK_NOT_NULL(T* ptr, const std::string& message, Args&&... args) {
  if (ptr == nullptr) {
    std::string formatted_msg = fmt::format(message, std::forward<Args>(args)...);
    throw MPSException(formatted_msg);
  }
}

/**
 * @brief 检查指针非空（无格式化版本）
 * @param ptr 要检查的指针
 * @param message 错误消息
 */
template<typename T>
inline void CHECK_NOT_NULL(T* ptr, const std::string& message) {
  if (ptr == nullptr) {
    throw MPSException(message);
  }
}

/**
 * @brief 检查数组边界
 * @param index 索引值
 * @param size 数组大小
 * @param array_name 数组名称（用于错误消息）
 */
template<typename... Args>
inline void CHECK_BOUNDS(size_t index, size_t size, const std::string& array_name, Args&&... args) {
  if (index >= size) {
    std::string formatted_msg = fmt::format(
        "数组边界检查失败: {} 索引 {} 超出范围 [0, {})", 
        array_name, index, size, std::forward<Args>(args)...);
    throw MPSException(formatted_msg);
  }
}

/**
 * @brief 检查数组边界（简化版本）
 * @param index 索引值
 * @param size 数组大小
 */
inline void CHECK_BOUNDS(size_t index, size_t size) {
  if (index >= size) {
    std::string formatted_msg = fmt::format(
        "数组边界检查失败: 索引 {} 超出范围 [0, {})", index, size);
    throw MPSException(formatted_msg);
  }
}

/**
 * @brief 检查文件操作是否成功
 * @param file_stream 文件流对象
 * @param filename 文件名（用于错误消息）
 * @param operation 操作类型（如"打开"、"读取"、"写入"）
 */
template<typename FileStream>
inline void CHECK_FILE(FileStream& file_stream, const std::string& filename, 
                       const std::string& operation) {
  if (!file_stream) {
    std::string formatted_msg = fmt::format(
        "文件操作失败: 无法{}文件 '{}'", operation, filename);
    throw MPSException(formatted_msg);
  }
}

/**
 * @brief 检查文件是否成功打开
 * @param file_stream 文件流对象
 * @param filename 文件名
 */
template<typename FileStream>
inline void CHECK_FILE_OPEN(FileStream& file_stream, const std::string& filename) {
  CHECK_FILE(file_stream, filename, "打开");
}

/**
 * @brief 检查数值范围
 * @param value 要检查的数值
 * @param min_val 最小值（包含）
 * @param max_val 最大值（不包含）
 * @param value_name 数值名称（用于错误消息）
 */
template<typename T>
inline void CHECK_RANGE(T value, T min_val, T max_val, const std::string& value_name) {
  if (value < min_val || value >= max_val) {
    std::string formatted_msg = fmt::format(
        "数值范围检查失败: {} = {} 超出范围 [{}, {})", 
        value_name, value, min_val, max_val);
    throw MPSException(formatted_msg);
  }
}

/**
 * @brief 检查数值是否为正数
 * @param value 要检查的数值
 * @param value_name 数值名称（用于错误消息）
 */
template<typename T>
inline void CHECK_POSITIVE(T value, const std::string& value_name) {
  if (value <= 0) {
    std::string formatted_msg = fmt::format(
        "数值检查失败: {} = {} 必须为正数", value_name, value);
    throw MPSException(formatted_msg);
  }
}

/**
 * @brief 检查数值是否为非负数
 * @param value 要检查的数值
 * @param value_name 数值名称（用于错误消息）
 */
template<typename T>
inline void CHECK_NON_NEGATIVE(T value, const std::string& value_name) {
  if (value < 0) {
    std::string formatted_msg = fmt::format(
        "数值检查失败: {} = {} 必须为非负数", value_name, value);
    throw MPSException(formatted_msg);
  }
}

} // namespace mps2D
