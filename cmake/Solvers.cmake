function(add_mfemsb_solver SOLVER_NAME)

  add_executable(
    ${SOLVER_NAME}
    ${CMAKE_CURRENT_SOURCE_DIR}/solvers/${SOLVER_NAME}/${SOLVER_NAME}.cpp)

  target_include_directories(
    ${SOLVER_NAME}
    PRIVATE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/solvers/${SOLVER_NAME})
  target_link_libraries(${SOLVER_NAME} PRIVATE ${MFEMSB_LIB_NAME} mfem::mfem)

  install(TARGETS ${SOLVER_NAME})
endfunction()
