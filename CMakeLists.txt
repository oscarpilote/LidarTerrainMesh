cmake_minimum_required(VERSION 3.13)
cmake_policy(VERSION 3.13)

project(swiss_lidar VERSION 0.1)

set( CMAKE_EXPORT_COMPILE_COMMANDS ON )
set( CMAKE_VERBOSE_MAKEFILE OFF )


add_subdirectory(src)
add_subdirectory(extern/miniply)
add_subdirectory(extern/meshoptimizer/src)
add_subdirectory(extern/lazperf)




