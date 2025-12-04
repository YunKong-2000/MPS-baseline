# MPS Baseline

MPS (Moving Particle Semi-implicit) 方法的基线实现。

## 项目结构

```
mps-baseline/
├── CMakeLists.txt          # 主 CMake 配置文件
├── cmake/                  # CMake 模块文件
│   ├── Common.cmake        # 通用配置（C++标准、编译选项、核心源文件）
│   ├── MainProgram.cmake   # 主程序配置
│   └── TestProgram.cmake   # 测试程序配置
├── src/                    # 源代码目录
│   ├── main.cpp           # 主程序入口
│   ├── core/              # 核心模块
│   │   ├── Particle.hpp
│   │   ├── Particle.cpp
│   │   ├── FileOperator.hpp
│   │   └── FileOperator.cpp
│   └── config/            # 配置模块
│       ├── ConfigManager.hpp
│       ├── ConfigManager.cpp
│       ├── MPSConfig.hpp
│       └── MPSConfig.cpp
├── include/               # 头文件目录
│   └── core/
│       ├── Types.h
│       └── Marcos.h
├── third_party/           # 第三方库目录
│   └── SimpleIni.h        # INI 解析库
├── test/                  # 测试代码目录
│   ├── test_particle.cpp # 粒子类测试
│   ├── test_config.cpp   # 配置读取测试
│   └── README.md         # 测试说明
├── data/                  # 测试数据目录
│   ├── fluid_particles.txt
│   ├── fluid_particles.csv
│   ├── solid_particles.txt
│   └── solid_particles.csv
└── config.ini             # 配置文件示例
```

## 构建说明

### 使用 CMake 构建

```bash
mkdir build
cd build
cmake ..
make
```

### 运行主程序

```bash
cd build/bin
./MPSBaseline
```

### 运行测试程序

```bash
cd build/bin
./test_particle    # 粒子类功能测试
./test_config      # 配置文件读取测试
```

## 配置文件使用

项目支持通过 INI 格式的配置文件进行参数设置。配置文件 `config.ini` 包含以下部分：

- **[File]**: 文件路径配置
  - `InputDir`: 输入目录
  - `OutputDir`: 输出目录
  - `FluidParticleFile`: 流体粒子文件名
  - `SolidParticleFile`: 固体粒子文件名

- **[Simulation]**: 仿真参数配置
  - `TimeStep`: 时间步长
  - `TotalTime`: 总仿真时间
  - `MaxIterations`: 最大迭代次数
  - `Density`: 流体密度
  - `Viscosity`: 动力粘度
  - `GravityX/Y/Z`: 重力加速度分量

- **[Particle]**: 粒子参数配置
  - `ParticleRadius`: 粒子半径
  - `SmoothingRadius`: 平滑半径
  - `ParticleCount`: 粒子数量

### 使用配置文件

```cpp
#include "config/MPSConfig.hpp"
#include "ini/SimpleIni.h"

using namespace mps;

// 直接使用 SimpleIni
SimpleIni ini;
if (ini.LoadFile("config.ini")) {
    std::string value = ini.GetValue("File", "InputDir", "data");
    int count = ini.GetIntValue("Particle", "ParticleCount", 1000);
    double time_step = ini.GetDoubleValue("Simulation", "TimeStep", 0.001);
}

// 使用参数类（推荐）
MPSConfig mps_config;
SimpleIni ini;
ini.LoadFile("config.ini");
mps_config.LoadFromConfig(ini);
const auto& sim_config = mps_config.GetSimulationConfig();
```

## CMake 结构说明

项目使用模块化的 CMake 配置：

- **CMakeLists.txt**: 主配置文件，包含项目基本信息和模块引用
- **cmake/Common.cmake**: 通用配置，包括：
  - C++ 标准设置
  - 编译选项
  - 核心源文件定义
  - 通用函数（如数据目录复制）
- **cmake/MainProgram.cmake**: 主程序配置
- **cmake/TestProgram.cmake**: 测试程序配置

这种结构使得 CMake 文件更加清晰和易于维护。

## 编程规范

- 遵循 Google C++ 编程规范
- 注重内存安全，防止内存泄漏
- 代码设计考虑未来扩展到 CPU/GPU 并行架构
- 尽量简化代码实现和第三方库调用，轻量化程序

