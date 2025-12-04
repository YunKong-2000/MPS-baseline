# 主程序配置
add_executable(${PROJECT_NAME}
    src/main.cpp
    ${CORE_SOURCES} # 核心模块
    ${CONFIG_SOURCES} # 配置模块
    ${NEIGHBOR_LIST_SOURCES} # 邻居列表模块
)

# 设置目标属性
target_include_directories(${PROJECT_NAME} PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)

# 复制数据目录和配置文件
copy_data_directory(${PROJECT_NAME})
file(COPY ${CMAKE_SOURCE_DIR}/config.ini DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# 如果需要调试信息
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(${PROJECT_NAME} PRIVATE DEBUG=1)
endif()

# 安装规则
install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION bin
)

