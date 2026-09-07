# partrac_add_app(<target> SOURCE <file> [NEEDS_DOLFIN])
function(partrac_add_app name)
  cmake_parse_arguments(APP "NEEDS_DOLFIN" "SOURCE" "" ${ARGN})
  if (APP_NEEDS_DOLFIN AND NOT PARTRAC_ENABLE_DOLFIN)
    message(STATUS "Skipping app ${name}: needs dolfin.")
    return()
  endif()

  add_executable(${name} ${APP_SOURCE})
  target_link_libraries(${name} PRIVATE ${PROJECT_NAME}_core)

  # in the build tree, so two configurations can coexist
  set_target_properties(${name} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin/")
  install(TARGETS ${name} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endfunction()
