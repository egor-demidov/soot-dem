cmake_minimum_required(VERSION 3.27)
project(soot_afm)

set(CMAKE_CXX_STANDARD 17)

include_directories(deps/libgran/include)
include_directories(deps/libtimestep/include)
include_directories(deps/eigen)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -flto=auto -O3 -march=native -fopenmp")

add_compile_definitions(LIBGRAN_USE_OMP)

add_executable(soot_afm main.cpp writer.cpp)
add_executable(aggregate_deposition aggregate_deposition.cpp writer.cpp)
