#pragma once
#include "core/Types.h"
#include <cmath>

namespace mps {

// MPS工具函数命名空间
// 提供粒子间距离计算和权重函数等通用工具函数

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
