cmake_minimum_required(VERSION 3.20)
project(soot_afm)

set(CMAKE_CXX_STANDARD 20)

include_directories(SYSTEM deps/libgran/include)
include_directories(SYSTEM deps/libtimestep/include)
include_directories(SYSTEM deps/eigen)
include_directories(SYSTEM deps)

add_compile_definitions(_USE_MATH_DEFINES)

if (${MSVC})
    set(CMAKE_CXX_FLAGS "/O2 /EHsc /GL /fp:except /Wall")
endif ()

if (${CMAKE_COMPILER_IS_GNUCXX})
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3")
    if (NOT ${MINGW})
        set(CMAKE_CXX_FLAGS "-flto=auto ${CMAKE_CXX_FLAGS}")
    endif ()
    if (${MINGW})
        set(CMAKE_CXX_FLAGS "-static ${CMAKE_CXX_FLAGS}")
    endif ()
endif ()

find_package(OpenMP REQUIRED)

if (${CMAKE_CXX_COMPILER_ID} MATCHES "AppleClang")
    message(FATAL_ERROR "AppleClang compiler is not supported")
endif ()

if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -flto=auto -O3 -fopenmp")
endif ()

include(FetchContent)
FetchContent_Declare(cmake_git_version_tracking
        GIT_REPOSITORY https://github.com/andrew-hardin/cmake-git-version-tracking.git
        GIT_TAG 6c0cb87edd029ddfb403a8e24577c144a03605a6
)
FetchContent_MakeAvailable(cmake_git_version_tracking)

set(LIBRARY_LIST cmake_git_version_tracking OpenMP::OpenMP_CXX ${LIBRARY_LIST})

add_compile_definitions(LIBGRAN_USE_OMP)

add_compile_definitions("COMPILER_STRING=\"${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}\"")
add_compile_definitions("SYSTEM_STRING=\"${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_VERSION}\"")

set(CMAKE_UTIL_SOURCES src/random_engine.cpp src/writer.cpp src/energy.cpp src/reader.cpp
        src/break_neck.cpp src/aggregate_stats.cpp
        src/parameter_loader.cpp src/io_common.cpp deps/tinyxml2/tinyxml2.cpp)

add_executable(aggregate_deposition src/aggregate_deposition.cpp ${CMAKE_UTIL_SOURCES})
add_executable(restructuring src/restructuring.cpp ${CMAKE_UTIL_SOURCES})
add_executable(afm_necking_fraction src/afm_necking_fraction.cpp ${CMAKE_UTIL_SOURCES})
add_executable(afm_breaking_necks src/afm_breaking_necks.cpp ${CMAKE_UTIL_SOURCES})
add_executable(aggregation src/aggregation.cpp ${CMAKE_UTIL_SOURCES})
add_executable(anchored_restructuring src/anchored_restructuring.cpp ${CMAKE_UTIL_SOURCES})
add_executable(restructuring_breaking_necks src/restructuring_breaking_necks.cpp ${CMAKE_UTIL_SOURCES})
add_executable(restructuring_breaking_necks_sphere_force src/restructuring_breaking_necks_sphere_force.cpp ${CMAKE_UTIL_SOURCES})
add_executable(anchored_restructuring_breaking_necks_sphere_force src/anchored_restructuring_breaking_necks_sphere_force.cpp ${CMAKE_UTIL_SOURCES})
add_executable(mechanics_testing src/mechanics_testing.cpp ${CMAKE_UTIL_SOURCES})

# VALIDATION EXECUTABLES
add_executable(01_particles_colliding validation/01_particles_colliding.cpp ${CMAKE_UTIL_SOURCES})
add_executable(02_particles_colliding validation/02_particles_colliding.cpp ${CMAKE_UTIL_SOURCES})
add_executable(03_particles_colliding validation/03_particles_colliding.cpp ${CMAKE_UTIL_SOURCES})
add_executable(04_particles_colliding validation/04_particles_colliding.cpp ${CMAKE_UTIL_SOURCES})
add_executable(05_hamaker validation/05_hamaker.cpp ${CMAKE_UTIL_SOURCES})
add_executable(06_sintering validation/06_sintering.cpp ${CMAKE_UTIL_SOURCES})
add_executable(07_sintering validation/07_sintering.cpp ${CMAKE_UTIL_SOURCES})
add_executable(08_sintering validation/08_sintering.cpp ${CMAKE_UTIL_SOURCES})
add_executable(09_sintering validation/09_sintering.cpp ${CMAKE_UTIL_SOURCES})

# NECK STRENGTH PARAMETRIZATION EXECUTABLE
add_executable(tensile_deformation neck_strength_parametrization/tensile_deformation.cpp
        neck_strength_parametrization/force_writer.cpp ${CMAKE_UTIL_SOURCES})

target_link_libraries(aggregate_deposition ${LIBRARY_LIST})
target_link_libraries(restructuring ${LIBRARY_LIST})
target_link_libraries(afm_necking_fraction ${LIBRARY_LIST})
target_link_libraries(afm_breaking_necks ${LIBRARY_LIST})
target_link_libraries(aggregation ${LIBRARY_LIST})
target_link_libraries(anchored_restructuring ${LIBRARY_LIST})
target_link_libraries(restructuring_breaking_necks ${LIBRARY_LIST})
target_link_libraries(restructuring_breaking_necks_sphere_force ${LIBRARY_LIST})
target_link_libraries(anchored_restructuring_breaking_necks_sphere_force ${LIBRARY_LIST})
target_link_libraries(mechanics_testing ${LIBRARY_LIST})

target_link_libraries(01_particles_colliding ${LIBRARY_LIST})
target_link_libraries(02_particles_colliding ${LIBRARY_LIST})
target_link_libraries(03_particles_colliding ${LIBRARY_LIST})
target_link_libraries(04_particles_colliding ${LIBRARY_LIST})
target_link_libraries(05_hamaker ${LIBRARY_LIST})
target_link_libraries(06_sintering ${LIBRARY_LIST})
target_link_libraries(07_sintering ${LIBRARY_LIST})
target_link_libraries(08_sintering ${LIBRARY_LIST})
target_link_libraries(09_sintering ${LIBRARY_LIST})

target_link_libraries(tensile_deformation ${LIBRARY_LIST})
