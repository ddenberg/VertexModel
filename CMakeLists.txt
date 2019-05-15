cmake_minimum_required(VERSION 3.5.1)
project(VertexModel)

set(CMAKE_CXX_FLAGS "-std=c++11 -fopenmp")

include_directories(include)

file(GLOB_RECURSE SOURCES "src/*.cpp")

# Add the executable
add_executable(VertexModel ${SOURCES})
