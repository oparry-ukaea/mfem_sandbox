cmake_minimum_required(VERSION 3.9...3.12)

# Choose library name
set(MFEMSB_LIB_NAME mfemsb)

# List of directories under ./solvers to build
set(SOLVERS ex9p)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
  cmake_policy(VERSION 3.12)
endif()

project(mfem_sandbox LANGUAGES CXX)

# Find CMake scripts
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Get the Git revision
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC MFEMSB_REVISION)
if(MFEMSB_REVISION STREQUAL "GITDIR-NOTFOUND")
  set(MFEMSB_REVISION "Unknown")
endif()
message(STATUS "Git revision: ${MFEMSB_REVISION}")

# ================================ Dependencies ===============================
find_package(MFEM REQUIRED)
find_package(MPI REQUIRED)

# ====================== MFEM-sandbox library (./common) ======================
set(MFEMSB_SOURCES common/tmp.cpp common/tmp.hpp)

add_library(${MFEMSB_LIB_NAME} ${MFEMSB_SOURCES})
target_link_libraries(${MFEMSB_LIB_NAME} PUBLIC mfem::mfem MPI::MPI_CXX)
target_compile_features(${MFEMSB_LIB_NAME} PUBLIC cxx_std_17)
set_target_properties(${MFEMSB_LIB_NAME} PROPERTIES CXX_EXTENSIONS OFF)

target_include_directories(
  ${MFEMSB_LIB_NAME}
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/common
         $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>/common
         $<INSTALL_INTERFACE:include>)

install(TARGETS ${MFEMSB_LIB_NAME})

# ================================== Solvers ==================================
# Build and install one exec per solver
include(Solvers)
foreach(SOLVER IN LISTS SOLVERS)
  add_mfemsb_solver(${SOLVER})
endforeach()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
