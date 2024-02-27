cmake_minimum_required(VERSION 3.27)
project(soot_afm)

set(CMAKE_CXX_STANDARD 17)

include_directories(deps/libgran/include)
include_directories(deps/libtimestep/include)
include_directories(deps/eigen)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -flto=auto -O3 -march=native -fopenmp")

add_compile_definitions(LIBGRAN_USE_OMP)

set(CMAKE_UTILITY_SOURCES energy.cpp writer.cpp reader.cpp)

add_executable(soot_afm main.cpp ${CMAKE_UTILITY_SOURCES})
add_executable(aggregate_deposition aggregate_deposition.cpp ${CMAKE_UTILITY_SOURCES})
add_executable(restructuring restructuring.cpp ${CMAKE_UTILITY_SOURCES})
