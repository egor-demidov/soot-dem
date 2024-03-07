cmake_minimum_required(VERSION 3.20)
project(soot_afm)

set(CMAKE_CXX_STANDARD 17)

include_directories(SYSTEM deps/libgran/include)
include_directories(SYSTEM deps/libtimestep/include)
include_directories(SYSTEM deps/eigen)
include_directories(SYSTEM deps)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -flto=auto -O3 -march=native -fopenmp")

include(FetchContent)
FetchContent_Declare(cmake_git_version_tracking
        GIT_REPOSITORY https://github.com/andrew-hardin/cmake-git-version-tracking.git
        GIT_TAG 904dbda1336ba4b9a1415a68d5f203f576b696bb
)
FetchContent_MakeAvailable(cmake_git_version_tracking)

add_compile_definitions(LIBGRAN_USE_OMP)

set(CMAKE_UTIL_SOURCES writer.cpp energy.cpp reader.cpp break_neck.cpp remove_overlap.cpp aggregate_stats.cpp
        parameter_loader.cpp io_common.cpp deps/tinyxml2/tinyxml2.cpp)

add_executable(soot_afm main.cpp ${CMAKE_UTIL_SOURCES})
add_executable(aggregate_deposition aggregate_deposition.cpp ${CMAKE_UTIL_SOURCES})
add_executable(restructuring restructuring.cpp ${CMAKE_UTIL_SOURCES})
add_executable(aggregate_deposition_no_necks aggregate_deposition_no_necks.cpp ${CMAKE_UTIL_SOURCES})
add_executable(afm_no_necks afm_no_necks.cpp ${CMAKE_UTIL_SOURCES})
add_executable(afm_necking_fraction afm_necking_fraction.cpp ${CMAKE_UTIL_SOURCES})

target_link_libraries(aggregate_deposition cmake_git_version_tracking)
