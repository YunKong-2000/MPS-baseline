# 2D子项目构建说明

本文档说明如何独立构建2D子项目。

## 独立构建2D子项目

2D子项目可以独立构建，拥有自己的构建目录，但会共享主项目的第三方库。

### 构建步骤

```bash
# 在项目根目录下
cd 2D

# 创建构建目录
mkdir build
cd build

# 配置CMake
cmake ..

# 编译
make

# 运行测试程序
cd bin
./test_types
./test_mps_utils
./test_particle
```

### 构建目录结构

```
2D/
├── CMakeLists.txt          # 主CMake配置文件
├── cmake/                  # CMake模块文件
│   ├── Common.cmake        # 通用配置
│   └── TestProgram.cmake   # 测试程序配置
├── build/                  # 构建目录（独立）
│   ├── bin/                # 可执行文件输出目录
│   │   ├── test_types
│   │   ├── test_mps_utils
│   │   └── test_particle
│   └── lib/                # 库文件输出目录
├── include/                # 头文件目录
│   └── core/
│       ├── Types.h
│       └── MPSUtils.h
├── src/                    # 源代码目录
│   └── core/
│       ├── Particle.hpp
│       └── Particle.cpp
└── test/                   # 测试代码目录
    ├── test_types.cpp
    ├── test_mps_utils.cpp
    └── test_particle.cpp
```

## 第三方库共享

2D子项目会自动使用主项目的第三方库：

- **third_party目录**: 通过相对路径 `../third_party` 访问主项目的第三方库
- **fmt库**: 如果系统未安装，会自动下载（与主项目使用相同的版本）

## CMake配置说明

### CMakeLists.txt
- 定义项目名称为 `MPSBaseline2D`
- 设置独立的输出目录
- 包含通用配置和测试程序配置

### cmake/Common.cmake
- 设置C++17标准
- 配置编译选项
- 设置包含目录（2D项目的include和主项目的third_party）
- 配置fmt库（与主项目相同）

### cmake/TestProgram.cmake
- 定义三个测试程序：test_types、test_mps_utils、test_particle
- 设置每个测试程序的包含目录和源文件
- 配置调试选项

## 与主项目的区别

| 特性 | 主项目 | 2D子项目 |
|------|--------|----------|
| 项目名称 | MPSBaseline | MPSBaseline2D |
| 构建目录 | build/ | 2D/build/ |
| 命名空间 | mps | mps2D |
| 向量类型 | double3, int3 | double2, int2 |
| 源文件 | src/core/ | 2D/src/core/ |
| 头文件 | include/core/ | 2D/include/core/ |
| 第三方库 | third_party/ | 共享主项目的third_party/ |

## 清理构建

```bash
cd 2D/build
rm -rf *
```

## 注意事项

1. **独立构建**: 2D子项目可以完全独立构建，不依赖主项目的构建
2. **共享库**: 第三方库（如fmt）会与主项目共享，避免重复下载
3. **路径**: 确保在2D目录下运行cmake，这样相对路径才能正确解析
4. **C++标准**: 使用C++17标准，与主项目一致
5. **FileOperator**: 当前 `FileOperator` 仅提供占位符实现，文件I/O功能尚未实现。`getParticleFromFile()` 方法会返回 -1 表示功能未实现

## 已知限制

- **FileOperator**: 2D版本的 `FileOperator` 类尚未完全实现，仅提供占位符接口用于编译通过
- **文件I/O**: 粒子数据文件读取/写入功能暂不可用
- **测试数据**: 当前测试程序使用程序内部生成的测试数据，不依赖外部文件

