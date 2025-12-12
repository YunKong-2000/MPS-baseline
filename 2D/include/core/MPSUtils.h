#pragma once
#include "core/Types.h"
#include <cmath>
#include <iostream>
#include <iomanip>

namespace mps2D {

// MPS工具函数命名空间（2D版本）
// 提供粒子间距离计算和权重函数等通用工具函数

// ============================================================================
// double2 基本运算
// ============================================================================

// 向量加法：对位相加
// 参数：
//   vec1: 第一个向量
//   vec2: 第二个向量
// 返回：两个向量对位相加的结果
inline double2 operator+(const double2& vec1, const double2& vec2) {
  return {vec1.x + vec2.x, vec1.y + vec2.y};
}

// 向量减法：对位相减
// 参数：
//   vec1: 第一个向量（被减数）
//   vec2: 第二个向量（减数）
// 返回：两个向量对位相减的结果
inline double2 operator-(const double2& vec1, const double2& vec2) {
  return {vec1.x - vec2.x, vec1.y - vec2.y};
}

// 向量除以标量
// 参数：
//   vec: 待除的向量
//   scalar: 标量值
// 返回：向量每个分量除以标量的结果
// 注意：如果 scalar 为 0，行为未定义
inline double2 operator/(const double2& vec, double scalar) {
  return {vec.x / scalar, vec.y / scalar};
}

// 向量内积（点积）
// 参数：
//   vec1: 第一个向量
//   vec2: 第二个向量
// 返回：两个向量的内积（标量）
inline double operator*(const double2& vec1, const double2& vec2) {
  return vec1.x * vec2.x + vec1.y * vec2.y;
}

// 向量赋值：从两个 double 值构造 double2
// 参数：
//   x: x 分量
//   y: y 分量
// 返回：构造的 double2 向量
inline double2 MakeDouble2(double x, double y) {
  return {x, y};
}

// 向量赋值：从另一个 double2 复制
// 参数：
//   vec: 源向量
// 返回：复制后的向量
inline double2 AssignDouble2(const double2& vec) {
  return {vec.x, vec.y};
}

// 打印 double2 向量到输出流
// 参数：
//   os: 输出流（如 std::cout）
//   vec: 待打印的向量
// 返回：输出流引用，支持链式调用
// 格式：{x, y}
inline std::ostream& operator<<(std::ostream& os, const double2& vec) {
  os << "{" << vec.x << ", " << vec.y << "}";
  return os;
}

// ============================================================================
// int2 基本运算
// ============================================================================

// 向量加法：对位相加
// 参数：
//   vec1: 第一个向量
//   vec2: 第二个向量
// 返回：两个向量对位相加的结果
inline int2 operator+(const int2& vec1, const int2& vec2) {
  return {vec1.x + vec2.x, vec1.y + vec2.y};
}

// 向量减法：对位相减
// 参数：
//   vec1: 第一个向量（被减数）
//   vec2: 第二个向量（减数）
// 返回：两个向量对位相减的结果
inline int2 operator-(const int2& vec1, const int2& vec2) {
  return {vec1.x - vec2.x, vec1.y - vec2.y};
}

// 向量内积（点积）
// 参数：
//   vec1: 第一个向量
//   vec2: 第二个向量
// 返回：两个向量的内积（标量）
inline int operator*(const int2& vec1, const int2& vec2) {
  return vec1.x * vec2.x + vec1.y * vec2.y;
}

// 向量赋值：从两个 int 值构造 int2
// 参数：
//   x: x 分量
//   y: y 分量
// 返回：构造的 int2 向量
inline int2 MakeInt2(int x, int y) {
  return {x, y};
}

// 向量赋值：从另一个 int2 复制
// 参数：
//   vec: 源向量
// 返回：复制后的向量
inline int2 AssignInt2(const int2& vec) {
  return {vec.x, vec.y};
}

// 打印 int2 向量到输出流
// 参数：
//   os: 输出流（如 std::cout）
//   vec: 待打印的向量
// 返回：输出流引用，支持链式调用
// 格式：{x, y}
inline std::ostream& operator<<(std::ostream& os, const int2& vec) {
  os << "{" << vec.x << ", " << vec.y << "}";
  return os;
}

// ============================================================================
// 其他工具函数
// ============================================================================

// 计算两点之间的欧氏距离（二维）
// 参数：
//   pos1: 第一个粒子的位置
//   pos2: 第二个粒子的位置
// 返回：两点之间的距离
inline double ComputeDistance(const double2& pos1, const double2& pos2) {
  double dx = pos1.x - pos2.x;
  double dy = pos1.y - pos2.y;
  return std::sqrt(dx * dx + dy * dy);
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

// 计算向量的模长（二维）
// 参数：
//   vec: 二维向量
// 返回：向量的模长
inline double ComputeVectorMagnitude(const double2& vec) {
  return std::sqrt(vec.x * vec.x + vec.y * vec.y);
}

// 归一化向量（二维）
// 参数：
//   vec: 待归一化的二维向量
// 返回：归一化后的向量，如果输入向量为零向量则返回零向量
inline double2 NormalizeVector(const double2& vec) {
  double magnitude = ComputeVectorMagnitude(vec);
  if (magnitude > 1e-10) {
    return {vec.x / magnitude, vec.y / magnitude};
  } else {
    return {0.0, 0.0};
  }
}

// 计算两个向量的点积（二维）
// 参数：
//   vec1: 第一个向量
//   vec2: 第二个向量
// 返回：两个向量的点积
inline double DotProduct(const double2& vec1, const double2& vec2) {
  return vec1.x * vec2.x + vec1.y * vec2.y;
}

} // namespace mps2D

