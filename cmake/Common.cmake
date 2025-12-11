# 通用配置
# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 设置编译选项
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# 包含目录
include_directories(${CMAKE_SOURCE_DIR}/include)
# third_party 目录用于包含第三方库（如 third_party/ini/SimpleIni.h）
include_directories(${CMAKE_SOURCE_DIR}/third_party)

# 核心源文件
set(CORE_SOURCES
    src/core/Particle.cpp
    src/core/FileOperator.cpp
)

# 配置模块源文件
set(CONFIG_SOURCES
    src/config/MPSConfig.cpp
)

# 邻居列表模块源文件
set(NEIGHBOR_LIST_SOURCES
    src/neighbour_list/NeighborListSearcher.cpp
)

# 自由面判定模块源文件
set(SURFACE_DETECTION_SOURCES
    src/surface_detection/SurfaceDetector.cpp
)

# 复制数据目录到可执行文件目录
function(copy_data_directory target_name)
    file(COPY ${CMAKE_SOURCE_DIR}/data DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endfunction()

