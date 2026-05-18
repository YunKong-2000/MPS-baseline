cmake_minimum_required(VERSION 3.10)
project(MPS2DTools LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(pressure_probe_postprocess
    pressure_probe_postprocess.cpp
)

add_executable(dambreak_front_postprocess
    dambreak_front_postprocess.cpp
)

add_executable(liquid_height_postprocess
    liquid_height_postprocess.cpp
)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(pressure_probe_postprocess PRIVATE DEBUG=1)
    target_compile_definitions(dambreak_front_postprocess PRIVATE DEBUG=1)
    target_compile_definitions(liquid_height_postprocess PRIVATE DEBUG=1)
endif()

