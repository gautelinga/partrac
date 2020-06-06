cmake_minimum_required(VERSION 3.0)
project(partrac)

set(CMAKE_CXX_COMPILER g++-9)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-pedantic-errors -Wall -Wextra -Werror -O3 -ffast-math -fconcepts -lstdc++")

find_package(HDF5 REQUIRED COMPONENTS CXX HL)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(${HDF5_INCLUDE_DIR})

aux_source_directory(src/ SRC_LIST)

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME} ${HDF5_LIBRARIES}
  ${HDF5_HL_LIBRARIES} Eigen3::Eigen)
