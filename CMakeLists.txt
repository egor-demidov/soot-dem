cmake_minimum_required(VERSION 3.20)
project(soot_afm)

set(CMAKE_CXX_STANDARD 17)

include_directories(SYSTEM deps/libgran/include)
include_directories(SYSTEM deps/libtimestep/include)
include_directories(SYSTEM deps/eigen)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -flto=auto -O3 -march=native -fopenmp")

add_compile_definitions(LIBGRAN_USE_OMP)

set(CMAKE_UTIL_SOURCES writer.cpp energy.cpp reader.cpp)

add_executable(soot_afm main.cpp ${CMAKE_UTIL_SOURCES})
add_executable(aggregate_deposition aggregate_deposition.cpp ${CMAKE_UTIL_SOURCES})
add_executable(restructuring restructuring.cpp ${CMAKE_UTIL_SOURCES})
