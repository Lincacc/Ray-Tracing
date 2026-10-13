cmake_minimum_required(VERSION 3.15)
project(RayTracing CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 1. 把 third_party 加入包含路径
include_directories(third_party)
include_directories(third_party/hosek)

# 2. 递归收集 src 目录下所有的 .cpp 和 .h 文件
file(GLOB_RECURSE SOURCES "src/*.cpp")
file(GLOB_RECURSE HEADERS "src/*.h")

# 关闭 c4819 警告
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4819")

# 3. 递归把 src 及其所有子目录都设为 include 路径
#    这样你在代码里既可以用 #include "core/vec3.h"，也可以用相对路径
include_directories(src)

# 4. 构建可执行文件
add_executable(RayTracing ${SOURCES} ${HEADERS}
    third_party/hosek/ArHosekSkyModel.cc)