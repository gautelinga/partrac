cmake_minimum_required(VERSION 3.14)

# Testing library
FetchContent_Declare(
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v2.13.7)
FetchContent_MakeAvailable(catch)
# Adds Catch2::Catch2

# Listed, not globbed: a glob is only re-evaluated when cmake runs again
add_executable(unittest
  test_main.cpp
  test_io.cpp
  test_h5direct.cpp
  test_params.cpp
  test_rng.cpp
  test_expressions.cpp
  test_schemas.cpp
  test_mesh_core.cpp
  test_mesh_tables.cpp
  test_param_files.cpp
  test_interpol_core.cpp
  test_cell_tree.cpp
  test_simplex_load.cpp
  test_stamp_buffer.cpp
  test_xdmf_load.cpp
  test_stamped_formats.cpp
  test_freq_load.cpp
  test_split_eval.cpp
  test_split_interpol.cpp
)

target_link_libraries(unittest
	PRIVATE
	Catch2::Catch2
	partrac_base
	partrac_io
	partrac_particles)
target_include_directories(unittest PRIVATE ../apps/)
target_compile_features(unittest PRIVATE cxx_std_17)
target_compile_definitions(unittest PRIVATE PARTRAC_SOURCE_DIR="${PROJECT_SOURCE_DIR}")

add_test(NAME unittesttest COMMAND unittest)

# The python tests drive the built apps, and skip themselves if bin/ is empty.
# Note find_package prefers an active virtualenv, so this is not always the
# system python; register the test only if that interpreter can run it.
find_package(Python3 COMPONENTS Interpreter)
if (Python3_Interpreter_FOUND)
  execute_process(COMMAND ${Python3_EXECUTABLE} -c "import pytest"
    RESULT_VARIABLE PARTRAC_NO_PYTEST OUTPUT_QUIET ERROR_QUIET)
endif()
if (Python3_Interpreter_FOUND AND NOT PARTRAC_NO_PYTEST)
  add_test(NAME pytest
    COMMAND ${Python3_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR} -q)
  # the apps live in the build tree, which the tests cannot guess
  set_tests_properties(pytest PROPERTIES
    ENVIRONMENT "PARTRAC_BIN=${PROJECT_BINARY_DIR}/bin")
elseif (Python3_Interpreter_FOUND)
  message(STATUS "No pytest for ${Python3_EXECUTABLE}, not running the python tests.")
endif()
