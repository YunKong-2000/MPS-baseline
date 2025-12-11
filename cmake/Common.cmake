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

# 配置fmt库（用于错误处理和格式化输出）
# 优先查找系统已安装的fmt库
find_package(fmt QUIET)
if(NOT fmt_FOUND)
    # 如果系统没有安装，则使用FetchContent下载
    include(FetchContent)
    # 使用commit hash而不是tag，更稳定可靠
    FetchContent_Declare(
      fmt
      GIT_REPOSITORY https://github.com/fmtlib/fmt.git
      GIT_TAG        e69e5f977d458f2650bb346dadf2ad30c5320281  # fmt 10.2.1的commit hash
      GIT_SHALLOW    TRUE    # 只克隆最新提交，加快下载速度
    )
    # 显示下载进度
    set(FETCHCONTENT_QUIET OFF)
    FetchContent_MakeAvailable(fmt)
endif()

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

