#pragma once
#include "core/Types.h"
#include <cmath>
#include <iostream>
#include <iomanip>

namespace mps {

// MPS工具函数命名空间
// 提供粒子间距离计算和权重函数等通用工具函数

// ============================================================================
// double3 基本运算
// ============================================================================

// 向量加法：对位相加
// 参数：
//   vec1: 第一个向量
//   vec2: 第二个向量
// 返回：两个向量对位相加的结果
inline double3 operator+(const double3& vec1, const double3& vec2) {
  return {vec1[0] + vec2[0], vec1[1] + vec2[1], vec1[2] + vec2[2]};
}

// 向量减法：对位相减
// 参数：
//   vec1: 第一个向量（被减数）
//   vec2: 第二个向量（减数）
// 返回：两个向量对位相减的结果
inline double3 operator-(const double3& vec1, const double3& vec2) {
  return {vec1[0] - vec2[0], vec1[1] - vec2[1], vec1[2] - vec2[2]};
}

// 向量除以标量
// 参数：
//   vec: 待除的向量
//   scalar: 标量值
// 返回：向量每个分量除以标量的结果
// 注意：如果 scalar 为 0，行为未定义
inline double3 operator/(const double3& vec, double scalar) {
  return {vec[0] / scalar, vec[1] / scalar, vec[2] / scalar};
}

// 向量内积（点积）
// 参数：
//   vec1: 第一个向量
//   vec2: 第二个向量
// 返回：两个向量的内积（标量）
inline double operator*(const double3& vec1, const double3& vec2) {
  return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

// 向量赋值：从三个 double 值构造 double3
// 参数：
//   x: x 分量
//   y: y 分量
//   z: z 分量
// 返回：构造的 double3 向量
inline double3 MakeDouble3(double x, double y, double z) {
  return {x, y, z};
}

// 向量赋值：从另一个 double3 复制
// 参数：
//   vec: 源向量
// 返回：复制后的向量
inline double3 AssignDouble3(const double3& vec) {
  return {vec[0], vec[1], vec[2]};
}

// 打印 double3 向量到输出流
// 参数：
//   os: 输出流（如 std::cout）
//   vec: 待打印的向量
// 返回：输出流引用，支持链式调用
// 格式：{x, y, z}
inline std::ostream& operator<<(std::ostream& os, const double3& vec) {
  os << "{" << vec[0] << ", " << vec[1] << ", " << vec[2] << "}";
  return os;
}

// ============================================================================
// 其他工具函数
// ============================================================================

// 计算两点之间的欧氏距离
// 参数：
//   pos1: 第一个粒子的位置
//   pos2: 第二个粒子的位置
// 返回：两点之间的距离
inline double ComputeDistance(const double3& pos1, const double3& pos2) {
  double dx = pos1[0] - pos2[0];
  double dy = pos1[1] - pos2[1];
  double dz = pos1[2] - pos2[2];
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// MPS标准权重函数
// 参数：
//   distance: 两点之间的距离
//   smoothing_radius: 平滑半径（r_e）
// 返回：权重值，当distance >= smoothing_radius时返回0.0
inline double WeightFunction(double distance, double smoothing_radius) {
  // MPS标准权重函数: w(r) = (1 - r/r_e)^2, r < r_e
  if (distance >= smoothing_radius) {
    return 0.0;
  }
  
  double r = distance / smoothing_radius;
  return (1.0 - r) * (1.0 - r);
}

// 计算向量的模长
// 参数：
//   vec: 三维向量
// 返回：向量的模长
inline double ComputeVectorMagnitude(const double3& vec) {
  return std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

// 归一化向量
// 参数：
//   vec: 待归一化的三维向量
// 返回：归一化后的向量，如果输入向量为零向量则返回零向量
inline double3 NormalizeVector(const double3& vec) {
  double magnitude = ComputeVectorMagnitude(vec);
  if (magnitude > 1e-10) {
    return {vec[0] / magnitude, vec[1] / magnitude, vec[2] / magnitude};
  } else {
    return {0.0, 0.0, 0.0};
  }
}

// 计算两个向量的点积
// 参数：
//   vec1: 第一个向量
//   vec2: 第二个向量
// 返回：两个向量的点积
inline double DotProduct(const double3& vec1, const double3& vec2) {
  return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

} // namespace mps
