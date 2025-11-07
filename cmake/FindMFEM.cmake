# FindMFEM.cmake
#
# Workaround the fact that the mfem spack package uses make, rather than cmake,
# so doesn't produce a .config.cmake file
#
# Sets up an interface target mfem::mfem.
#
# The directories containing include/mfem.hpp and lib/libmfem.so or
# lib64/libmfem.so must be in CMAKE_PREFIX_PATH for this to work.

find_path(
  MFEM_INC
  NAMES "mfem.hpp"
  HINTS ${CMAKE_PREFIX_PATH}
  PATH_SUFFIXES "include" REQUIRED)

find_library(
  MFEM_LIB
  NAMES "mfem"
  HINTS ${CMAKE_PREFIX_PATH}
  PATH_SUFFIXES "lib" "lib64" REQUIRED)

add_library(mfem::mfem INTERFACE IMPORTED)
target_link_libraries(mfem::mfem INTERFACE ${MFEM_LIB} -lmfem -lHYPRE)
target_include_directories(mfem::mfem INTERFACE ${MFEM_INC} ${MFEM_INC}/mfem)
