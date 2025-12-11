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

add_executable(test_surface_detection
    test/test_surface_detection.cpp
    ${CORE_SOURCES}
    ${CONFIG_SOURCES}
    ${NEIGHBOR_LIST_SOURCES}
    ${SURFACE_DETECTION_SOURCES}
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

target_include_directories(test_surface_detection PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)

# 链接fmt库到所有测试目标
target_link_libraries(test_particle PRIVATE fmt::fmt)
target_link_libraries(test_config PRIVATE fmt::fmt)
target_link_libraries(test_neighbor_list PRIVATE fmt::fmt)
target_link_libraries(test_surface_detection PRIVATE fmt::fmt)

# 复制数据目录和配置文件
copy_data_directory(test_particle)
copy_data_directory(test_config)
copy_data_directory(test_neighbor_list)
copy_data_directory(test_surface_detection)
file(COPY ${CMAKE_SOURCE_DIR}/config.ini DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

