cmake_minimum_required(VERSION 3.20)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(sphinxsys_capillary_validation LANGUAGES CXX)
endif()
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# The capillary pair law and its checks are usable without SPHinXsys.
add_library(pf_capillary INTERFACE)
target_include_directories(pf_capillary INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/src")
add_executable(pf_capillary_test tests/test_pf_capillary.cpp)
target_link_libraries(pf_capillary_test PRIVATE pf_capillary)
enable_testing()
add_test(NAME pf_capillary_pair_regression COMMAND pf_capillary_test)
add_executable(pf_liquid_test tests/test_pf_liquid.cpp)
target_link_libraries(pf_liquid_test PRIVATE pf_capillary)
add_test(NAME pf_liquid_regression COMMAND pf_liquid_test)
add_executable(pf_restart_test tests/test_pf_restart.cpp)
target_link_libraries(pf_restart_test PRIVATE pf_capillary)
add_test(NAME pf_restart_regression COMMAND pf_restart_test)
add_executable(pf_snapshot_test tests/test_pf_snapshot.cpp)
target_link_libraries(pf_snapshot_test PRIVATE pf_capillary)
add_test(NAME pf_snapshot_regression COMMAND pf_snapshot_test)
add_executable(pf_rigid_test tests/test_pf_rigid.cpp)
target_link_libraries(pf_rigid_test PRIVATE pf_capillary)
add_test(NAME pf_rigid_regression COMMAND pf_rigid_test)
add_executable(pf_bridge_test tests/test_pf_bridge.cpp)
target_link_libraries(pf_bridge_test PRIVATE pf_capillary)
add_test(NAME pf_bridge_regression COMMAND pf_bridge_test)
add_executable(pf_neighbors_test tests/test_pf_neighbors.cpp)
target_link_libraries(pf_neighbors_test PRIVATE pf_capillary)
add_test(NAME pf_neighbors_regression COMMAND pf_neighbors_test)
# General N-body DEM driver: same PF-WCSPH core, independent of SPHinXsys I/O.
add_executable(pf_dem src/pf_dem.cpp)
target_link_libraries(pf_dem PRIVATE pf_capillary)
find_package(TBB QUIET COMPONENTS tbb)
if(TARGET TBB::tbb)
  target_compile_definitions(pf_dem PRIVATE PF_DEM_USE_TBB=1)
  target_link_libraries(pf_dem PRIVATE TBB::tbb)
endif()
add_executable(pf_dem_test tests/test_pf_dem.cpp)
target_link_libraries(pf_dem_test PRIVATE pf_capillary)
add_test(NAME pf_dem_regression COMMAND pf_dem_test)
option(PF_WITH_SOOT_DEM "Build actual soot-dem contact/neck backend" OFF)
set(SOOT_DEM_SOURCE_DIR "" CACHE PATH "Existing soot-dem checkout with libgran and Eigen initialized")
if(PF_WITH_SOOT_DEM)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  if(NOT EXISTS "${SOOT_DEM_SOURCE_DIR}/src/aggregate.h")
    message(FATAL_ERROR "Set SOOT_DEM_SOURCE_DIR to your soot-dem checkout")
  endif()
  set(PF_SOOT_OVERLAY "${CMAKE_CURRENT_BINARY_DIR}/soot-backend")
  execute_process(COMMAND ${Python3_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/tools/prepare_soot_backend.py"
    --source "${SOOT_DEM_SOURCE_DIR}" --output "${PF_SOOT_OVERLAY}"
    COMMAND_ERROR_IS_FATAL ANY)
  file(GLOB_RECURSE PF_SOOT_HEADERS CONFIGURE_DEPENDS "${SOOT_DEM_SOURCE_DIR}/deps/libgran/include/*.h")
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    ${PF_SOOT_HEADERS}
    "${SOOT_DEM_SOURCE_DIR}/src/aggregate.h"
    "${SOOT_DEM_SOURCE_DIR}/src/break_neck.h"
    "${SOOT_DEM_SOURCE_DIR}/deps/libgran/include/libgran/contact_force/contact_force.h"
    "${SOOT_DEM_SOURCE_DIR}/deps/libgran/include/libgran/sinter_bridge/alt_sinter_bridge.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/tools/prepare_soot_backend.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/pf_soot_dem.h")
  add_library(pf_soot_backend INTERFACE)
  target_compile_features(pf_soot_backend INTERFACE cxx_std_20)
  target_include_directories(pf_soot_backend SYSTEM INTERFACE
    "${PF_SOOT_OVERLAY}" "${SOOT_DEM_SOURCE_DIR}/deps/libgran/include" "${SOOT_DEM_SOURCE_DIR}/deps/eigen")
  target_link_libraries(pf_soot_backend INTERFACE pf_capillary)
  add_executable(pf_soot_dem src/pf_dem.cpp)
  target_compile_definitions(pf_soot_dem PRIVATE PF_WITH_SOOT_DEM=1)
  target_link_libraries(pf_soot_dem PRIVATE pf_soot_backend)
  target_include_directories(pf_soot_dem SYSTEM BEFORE PRIVATE "${PF_SOOT_OVERLAY}")
  if(TARGET TBB::tbb)
    target_compile_definitions(pf_soot_dem PRIVATE PF_DEM_USE_TBB=1)
    target_link_libraries(pf_soot_dem PRIVATE TBB::tbb)
  endif()
  add_executable(pf_soot_dem_test tests/test_pf_soot_dem.cpp)
  target_link_libraries(pf_soot_dem_test PRIVATE pf_soot_backend)
  target_include_directories(pf_soot_dem_test SYSTEM BEFORE PRIVATE "${PF_SOOT_OVERLAY}")
  add_test(NAME pf_soot_dem_regression COMMAND pf_soot_dem_test)
endif()
add_executable(capillary_clock_test tests/test_time_window.cpp)
add_test(NAME capillary_clock_regression COMMAND capillary_clock_test)

add_executable(capillary_viscous_test tests/test_viscous_stability.cpp)
add_test(NAME capillary_viscous_regression COMMAND capillary_viscous_test)
# SPHinXsys is only used by the optional historical validation adapters.
# Default PF-SPH/DEM targets need neither this checkout nor a bootstrap step.
add_executable(pf_density_diagnostics_test tests/test_pf_density_diagnostics.cpp)
target_link_libraries(pf_density_diagnostics_test PRIVATE pf_capillary)
add_test(NAME pf_density_diagnostics_regression COMMAND pf_density_diagnostics_test)

option(PF_BUILD_SPHINXSYS_VALIDATION "Build optional SPHinXsys validation adapters" OFF)
option(CAPILLARY_PAIR_TESTS_ONLY "Build dependency-free PF/DEM drivers and checks only" OFF)
if(CAPILLARY_PAIR_TESTS_ONLY OR NOT PF_BUILD_SPHINXSYS_VALIDATION)
  return()
endif()
set(SPHINXSYS_PIN "319520e02f631847f85160c0e8b55648c395eafa")
set(SPHINXSYS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/SPHinXsys" CACHE PATH "Pinned SPHinXsys checkout")
if(NOT EXISTS "${SPHINXSYS_SOURCE_DIR}/CMakeLists.txt")
  message(FATAL_ERROR
    "PF_BUILD_SPHINXSYS_VALIDATION=ON requires an existing SPHinXsys checkout. "
    "Set SPHINXSYS_SOURCE_DIR, or disable PF_BUILD_SPHINXSYS_VALIDATION to build PF-SPH/DEM only.")
endif()
find_package(Git REQUIRED)
execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD WORKING_DIRECTORY ${SPHINXSYS_SOURCE_DIR} OUTPUT_VARIABLE actual_pin OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
if(NOT actual_pin STREQUAL SPHINXSYS_PIN)
  message(FATAL_ERROR "SPHinXsys checkout does not match pinned version ${SPHINXSYS_PIN}")
endif()
set(SPHINXSYS_2D OFF CACHE BOOL "" FORCE)
set(SPHINXSYS_3D ON CACHE BOOL "" FORCE)
set(SPHINXSYS_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(SPHINXSYS_DEVELOPER_MODE OFF CACHE BOOL "" FORCE)
add_subdirectory(${SPHINXSYS_SOURCE_DIR} _deps/sphinxsys EXCLUDE_FROM_ALL)
add_executable(pf_capillary_validation src/pf_validation.cpp)
target_link_libraries(pf_capillary_validation PRIVATE sphinxsys_3d pf_capillary)
add_test(NAME pf_liquid_adapter_smoke COMMAND pf_capillary_validation --dx 0.3 --end-time 0.0001 --steps 2 --dump-every 0)
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pf-smoke")
set_tests_properties(pf_liquid_adapter_smoke PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pf-smoke" TIMEOUT 60)
add_test(NAME pf_sessile_adapter_smoke COMMAND pf_capillary_validation --case sessile --dx 0.3 --end-time 0.0001 --steps 2 --dump-every 0)
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pf-sessile-smoke")
set_tests_properties(pf_sessile_adapter_smoke PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pf-sessile-smoke" TIMEOUT 60)
# Experimental 6DOF DEM/SPH adapter: native build required to execute this smoke.
add_test(NAME pf_immersed_adapter_smoke COMMAND pf_capillary_validation
  --case immersed --dx 0.12 --radius 0.8 --sphere-radius 0.24
  --sphere-density 2 --sphere-vx 0.05 --sphere-omega-z 0.5
  --end-time 0.002 --steps 10 --dump-every 0 --record-every 1
  --checkpoint-interval 0.002)
# Experimental native two-fixed-sphere bridge; measurement is not yet a
# calibrated Young-Laplace capillary force.
add_test(NAME pf_bridge_adapter_smoke COMMAND pf_capillary_validation
  --case bridge --dx 0.1 --radius 0.35 --sphere-radius 0.3
  --bridge-gap 0.2 --fixed-spheres 1 --sff 0.001 --sfs 0.002
  --end-time 0.001 --steps 4 --dump-every 0 --record-every 1
  --checkpoint-interval 0.001)
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pf-bridge-smoke")
set_tests_properties(pf_bridge_adapter_smoke PROPERTIES
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pf-bridge-smoke" TIMEOUT 60)
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pf-immersed-smoke")
set_tests_properties(pf_immersed_adapter_smoke PROPERTIES
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pf-immersed-smoke" TIMEOUT 60)
add_executable(capillary_validation src/validation.cpp)
target_link_libraries(capillary_validation PRIVATE sphinxsys_3d)
enable_testing()
add_test(NAME capillary_smoke COMMAND capillary_validation --case laplace --dx 0.3 --steps 2 --end-time 0.01 --dump-every 0)
set_tests_properties(capillary_smoke PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} TIMEOUT 60)


