# 测试程序配置
add_executable(test_particle
    test/test_particle.cpp
    ${CORE_SOURCES}
    ${CONFIG_SOURCES}
    ${NEIGHBOR_LIST_SOURCES}
)

add_executable(test_config
    test/test_config.cpp
    ${CORE_SOURCES}
    ${CONFIG_SOURCES}
    ${NEIGHBOR_LIST_SOURCES}
)

add_executable(test_neighbor_list
    test/test_neighbor_list.cpp
    ${CORE_SOURCES}
    ${CONFIG_SOURCES}
    ${NEIGHBOR_LIST_SOURCES}
)

# 设置测试目标属性
target_include_directories(test_particle PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)

target_include_directories(test_config PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)

target_include_directories(test_neighbor_list PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)

# 复制数据目录和配置文件
copy_data_directory(test_particle)
copy_data_directory(test_config)
copy_data_directory(test_neighbor_list)
file(COPY ${CMAKE_SOURCE_DIR}/config.ini DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

