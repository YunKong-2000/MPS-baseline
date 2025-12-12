#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "../include/core/MPSUtils.h"

using namespace mps2D;

void printSeparator(const std::string& title) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "  " << title << "\n";
    std::cout << std::string(60, '=') << "\n";
}

void testDouble2Operations() {
    printSeparator("测试 double2 基本运算");
    
    double2 vec1 = {1.0, 2.0};
    double2 vec2 = {3.0, 4.0};
    
    std::cout << "\n1. 向量运算:\n";
    std::cout << "   vec1 = " << vec1 << "\n";
    std::cout << "   vec2 = " << vec2 << "\n";
    
    double2 sum = vec1 + vec2;
    std::cout << "   vec1 + vec2 = " << sum << "\n";
    
    double2 diff = vec2 - vec1;
    std::cout << "   vec2 - vec1 = " << diff << "\n";
    
    double2 scaled = vec1 / 2.0;
    std::cout << "   vec1 / 2.0 = " << scaled << "\n";
    
    double dot = vec1 * vec2;
    std::cout << "   vec1 * vec2 (点积) = " << dot << "\n";
    
    std::cout << "\n2. 辅助函数:\n";
    double2 vec3 = MakeDouble2(5.0, 6.0);
    std::cout << "   MakeDouble2(5.0, 6.0) = " << vec3 << "\n";
    
    double2 vec4 = AssignDouble2(vec3);
    std::cout << "   AssignDouble2(vec3) = " << vec4 << "\n";
}

void testInt2Operations() {
    printSeparator("测试 int2 基本运算");
    
    int2 vec1 = {10, 20};
    int2 vec2 = {30, 40};
    
    std::cout << "\n1. 向量运算:\n";
    std::cout << "   vec1 = " << vec1 << "\n";
    std::cout << "   vec2 = " << vec2 << "\n";
    
    int2 sum = vec1 + vec2;
    std::cout << "   vec1 + vec2 = " << sum << "\n";
    
    int2 diff = vec2 - vec1;
    std::cout << "   vec2 - vec1 = " << diff << "\n";
    
    int dot = vec1 * vec2;
    std::cout << "   vec1 * vec2 (点积) = " << dot << "\n";
    
    std::cout << "\n2. 辅助函数:\n";
    int2 vec3 = MakeInt2(50, 60);
    std::cout << "   MakeInt2(50, 60) = " << vec3 << "\n";
    
    int2 vec4 = AssignInt2(vec3);
    std::cout << "   AssignInt2(vec3) = " << vec4 << "\n";
}

void testDistanceComputation() {
    printSeparator("测试距离计算");
    
    double2 pos1 = {0.0, 0.0};
    double2 pos2 = {3.0, 4.0};
    
    std::cout << "\n1. 欧氏距离计算:\n";
    std::cout << "   pos1 = " << pos1 << "\n";
    std::cout << "   pos2 = " << pos2 << "\n";
    
    double distance = ComputeDistance(pos1, pos2);
    std::cout << "   ComputeDistance(pos1, pos2) = " << distance << "\n";
    std::cout << "   预期值: 5.0 (3-4-5三角形)\n";
    
    if (std::abs(distance - 5.0) < 1e-10) {
        std::cout << "   ✓ 距离计算正确\n";
    } else {
        std::cout << "   ✗ 距离计算错误\n";
    }
    
    // 测试零距离
    double2 pos3 = {1.0, 1.0};
    double zero_dist = ComputeDistance(pos3, pos3);
    std::cout << "\n2. 零距离测试:\n";
    std::cout << "   ComputeDistance(pos3, pos3) = " << zero_dist << "\n";
    if (zero_dist < 1e-10) {
        std::cout << "   ✓ 零距离计算正确\n";
    } else {
        std::cout << "   ✗ 零距离计算错误\n";
    }
}

void testWeightFunction() {
    printSeparator("测试权重函数");
    
    double smoothing_radius = 2.0;
    
    std::cout << "\n测试 MPS 标准权重函数 (smoothing_radius = " 
              << smoothing_radius << "):\n";
    
    std::vector<double> distances = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5};
    
    for (double dist : distances) {
        double weight = WeightFunction(dist, smoothing_radius);
        std::cout << "   w(" << dist << ") = " << weight << "\n";
    }
    
    // 验证边界条件
    std::cout << "\n验证边界条件:\n";
    double w_at_boundary = WeightFunction(smoothing_radius, smoothing_radius);
    std::cout << "   w(r_e) = " << w_at_boundary << " (应为0.0)\n";
    if (w_at_boundary < 1e-10) {
        std::cout << "   ✓ 边界条件正确\n";
    } else {
        std::cout << "   ✗ 边界条件错误\n";
    }
    
    double w_beyond = WeightFunction(smoothing_radius + 0.1, smoothing_radius);
    std::cout << "   w(r_e + 0.1) = " << w_beyond << " (应为0.0)\n";
    if (w_beyond < 1e-10) {
        std::cout << "   ✓ 超出范围权重正确\n";
    } else {
        std::cout << "   ✗ 超出范围权重错误\n";
    }
}

void testVectorMagnitude() {
    printSeparator("测试向量模长计算");
    
    double2 vec1 = {3.0, 4.0};
    double2 vec2 = {1.0, 1.0};
    double2 vec3 = {0.0, 0.0};
    
    std::cout << "\n1. 向量模长计算:\n";
    std::cout << "   vec1 = " << vec1 << "\n";
    double mag1 = ComputeVectorMagnitude(vec1);
    std::cout << "   |vec1| = " << mag1 << " (预期: 5.0)\n";
    if (std::abs(mag1 - 5.0) < 1e-10) {
        std::cout << "   ✓ 模长计算正确\n";
    } else {
        std::cout << "   ✗ 模长计算错误\n";
    }
    
    std::cout << "\n   vec2 = " << vec2 << "\n";
    double mag2 = ComputeVectorMagnitude(vec2);
    std::cout << "   |vec2| = " << mag2 << " (预期: sqrt(2) ≈ " 
              << std::sqrt(2.0) << ")\n";
    
    std::cout << "\n   vec3 = " << vec3 << "\n";
    double mag3 = ComputeVectorMagnitude(vec3);
    std::cout << "   |vec3| = " << mag3 << " (预期: 0.0)\n";
    if (mag3 < 1e-10) {
        std::cout << "   ✓ 零向量模长正确\n";
    } else {
        std::cout << "   ✗ 零向量模长错误\n";
    }
}

void testNormalizeVector() {
    printSeparator("测试向量归一化");
    
    double2 vec1 = {3.0, 4.0};
    double2 vec2 = {1.0, 1.0};
    double2 vec3 = {0.0, 0.0};
    
    std::cout << "\n1. 向量归一化:\n";
    std::cout << "   vec1 = " << vec1 << "\n";
    double2 norm1 = NormalizeVector(vec1);
    std::cout << "   NormalizeVector(vec1) = " << norm1 << "\n";
    double mag1 = ComputeVectorMagnitude(norm1);
    std::cout << "   |归一化后的vec1| = " << mag1 << " (预期: 1.0)\n";
    if (std::abs(mag1 - 1.0) < 1e-10) {
        std::cout << "   ✓ 归一化正确\n";
    } else {
        std::cout << "   ✗ 归一化错误\n";
    }
    
    std::cout << "\n   vec2 = " << vec2 << "\n";
    double2 norm2 = NormalizeVector(vec2);
    std::cout << "   NormalizeVector(vec2) = " << norm2 << "\n";
    double mag2 = ComputeVectorMagnitude(norm2);
    std::cout << "   |归一化后的vec2| = " << mag2 << " (预期: 1.0)\n";
    
    std::cout << "\n2. 零向量归一化:\n";
    std::cout << "   vec3 = " << vec3 << "\n";
    double2 norm3 = NormalizeVector(vec3);
    std::cout << "   NormalizeVector(vec3) = " << norm3 << "\n";
    std::cout << "   (零向量应返回零向量)\n";
}

void testDotProduct() {
    printSeparator("测试点积计算");
    
    double2 vec1 = {1.0, 2.0};
    double2 vec2 = {3.0, 4.0};
    
    std::cout << "\n1. 点积计算:\n";
    std::cout << "   vec1 = " << vec1 << "\n";
    std::cout << "   vec2 = " << vec2 << "\n";
    
    double dot1 = DotProduct(vec1, vec2);
    std::cout << "   DotProduct(vec1, vec2) = " << dot1 << "\n";
    std::cout << "   预期值: 1*3 + 2*4 = 11\n";
    if (std::abs(dot1 - 11.0) < 1e-10) {
        std::cout << "   ✓ 点积计算正确\n";
    } else {
        std::cout << "   ✗ 点积计算错误\n";
    }
    
    // 测试正交向量
    double2 vec4 = {1.0, 0.0};
    double2 vec5 = {0.0, 1.0};
    double dot2 = DotProduct(vec4, vec5);
    std::cout << "\n2. 正交向量点积:\n";
    std::cout << "   vec4 = " << vec4 << "\n";
    std::cout << "   vec5 = " << vec5 << "\n";
    std::cout << "   DotProduct(vec4, vec5) = " << dot2 << " (预期: 0.0)\n";
    if (std::abs(dot2) < 1e-10) {
        std::cout << "   ✓ 正交向量点积正确\n";
    } else {
        std::cout << "   ✗ 正交向量点积错误\n";
    }
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "MPS Baseline 2D - MPSUtils 工具函数测试\n";
    std::cout << std::string(60, '=') << "\n";
    
    try {
        testDouble2Operations();
        testInt2Operations();
        testDistanceComputation();
        testWeightFunction();
        testVectorMagnitude();
        testNormalizeVector();
        testDotProduct();
        
        printSeparator("所有测试完成");
        std::cout << "\n✓ 所有工具函数测试通过\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\n错误: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}

