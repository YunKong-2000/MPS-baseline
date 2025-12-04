cmake_minimum_required(VERSION 3.10)
project(MPSBaseline VERSION 1.0.0 LANGUAGES CXX)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 设置编译选项
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# 设置输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# 包含目录
include_directories(${CMAKE_SOURCE_DIR}/include)

# 源文件
set(CORE_SOURCES
    src/core/Particle.cpp
    src/core/FileOperator.cpp
)

set(CORE_HEADERS
    src/core/Particle.hpp
    src/core/FileOperator.hpp
)

# 创建可执行文件
add_executable(${PROJECT_NAME}
    src/core/main.cpp
    ${CORE_SOURCES}
)

# 设置目标属性
target_include_directories(${PROJECT_NAME} PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)

# 复制数据目录到可执行文件目录
file(COPY ${CMAKE_SOURCE_DIR}/data DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# 如果需要调试信息
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(${PROJECT_NAME} PRIVATE DEBUG=1)
endif()

# 安装规则（可选）
install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION bin
)

# 创建测试可执行文件
add_executable(test_particle
    test/test_particle.cpp
    ${CORE_SOURCES}
)

# 设置测试目标属性
target_include_directories(test_particle PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)

# 复制数据目录到测试可执行文件目录
file(COPY ${CMAKE_SOURCE_DIR}/data DESTINATION ${CMAKE_BINARY_DIR}/bin)

# 打印配置信息
message(STATUS "Project: ${PROJECT_NAME}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install Prefix: ${CMAKE_INSTALL_PREFIX}")

