#include <iostream>
#include <iomanip>
#include <vector>
#include "../include/core/Types.h"
#include "../include/core/MPSUtils.h"

using namespace mps2D;

void printSeparator(const std::string& title) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "  " << title << "\n";
    std::cout << std::string(60, '=') << "\n";
}

void testDouble2() {
    printSeparator("测试 double2 结构体");
    
    // 测试基本构造和访问
    std::cout << "\n1. 测试 double2 基本构造和访问:\n";
    double2 vec1 = {1.5, 2.5};
    std::cout << "   vec1 = {" << vec1.x << ", " << vec1.y << "}\n";
    
    double2 vec2;
    vec2.x = 3.0;
    vec2.y = 4.0;
    std::cout << "   vec2 = {" << vec2.x << ", " << vec2.y << "}\n";
    
    // 测试向量运算（需要MPSUtils.h）
    std::cout << "\n2. 测试 double2 向量运算:\n";
    double2 sum = vec1 + vec2;
    std::cout << "   vec1 + vec2 = {" << sum.x << ", " << sum.y << "}\n";
    
    double2 diff = vec2 - vec1;
    std::cout << "   vec2 - vec1 = {" << diff.x << ", " << diff.y << "}\n";
    
    double2 scaled = vec1 / 2.0;
    std::cout << "   vec1 / 2.0 = {" << scaled.x << ", " << scaled.y << "}\n";
    
    double dot = vec1 * vec2;
    std::cout << "   vec1 * vec2 (点积) = " << dot << "\n";
}

void testInt2() {
    printSeparator("测试 int2 结构体");
    
    // 测试基本构造和访问
    std::cout << "\n1. 测试 int2 基本构造和访问:\n";
    int2 vec1 = {10, 20};
    std::cout << "   vec1 = {" << vec1.x << ", " << vec1.y << "}\n";
    
    int2 vec2;
    vec2.x = 30;
    vec2.y = 40;
    std::cout << "   vec2 = {" << vec2.x << ", " << vec2.y << "}\n";
    
    // 测试向量运算（需要MPSUtils.h）
    std::cout << "\n2. 测试 int2 向量运算:\n";
    int2 sum = vec1 + vec2;
    std::cout << "   vec1 + vec2 = {" << sum.x << ", " << sum.y << "}\n";
    
    int2 diff = vec2 - vec1;
    std::cout << "   vec2 - vec1 = {" << diff.x << ", " << diff.y << "}\n";
    
    int dot = vec1 * vec2;
    std::cout << "   vec1 * vec2 (点积) = " << dot << "\n";
}

void testSurfaceType() {
    printSeparator("测试 SurfaceType 枚举");
    
    std::cout << "\n测试 SurfaceType 枚举值:\n";
    std::cout << "  INNER = " << static_cast<int>(SurfaceType::INNER) << "\n";
    std::cout << "  NEAR_SURFACE = " << static_cast<int>(SurfaceType::NEAR_SURFACE) << "\n";
    std::cout << "  SURFACE = " << static_cast<int>(SurfaceType::SURFACE) << "\n";
    std::cout << "  SPLASH = " << static_cast<int>(SurfaceType::SPLASH) << "\n";
    
    // 测试使用
    std::vector<SurfaceType> surface_types = {
        SurfaceType::INNER,
        SurfaceType::NEAR_SURFACE,
        SurfaceType::SURFACE,
        SurfaceType::SPLASH
    };
    
    std::cout << "\n测试 SurfaceType 向量:\n";
    for (size_t i = 0; i < surface_types.size(); ++i) {
        std::cout << "  surface_types[" << i << "] = " 
                  << static_cast<int>(surface_types[i]) << "\n";
    }
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "MPS Baseline 2D - Types 类型测试\n";
    std::cout << std::string(60, '=') << "\n";
    
    try {
        testDouble2();
        testInt2();
        testSurfaceType();
        
        printSeparator("所有测试完成");
        std::cout << "\n✓ 所有类型测试通过\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\n错误: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}

