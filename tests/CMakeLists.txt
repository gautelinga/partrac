cmake_minimum_required(VERSION 3.14)

# Testing library
FetchContent_Declare(
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v2.13.7)
FetchContent_MakeAvailable(catch)
# Adds Catch2::Catch2

set (WILDCARD_SOURCE test_*.cpp)

file (GLOB_RECURSE TEST_SOURCES ${WILDCARD_SOURCE})

add_executable(unittest
  ${TEST_SOURCES}
)

if (PARTRAC_ENABLE_DOLFIN)
  target_link_libraries(${PROJECT_NAME} dolfin)
endif()
target_link_libraries(${PROJECT_NAME} ${HDF5_LIBRARIES}
  ${HDF5_HL_LIBRARIES} Eigen3::Eigen)
target_link_libraries(unittest 
	PRIVATE 
	Catch2::Catch2
	partrac_io)
target_include_directories(unittest PRIVATE ../src/)
target_compile_features(unittest PRIVATE cxx_std_17)

add_test(NAME unittesttest COMMAND unittest)
