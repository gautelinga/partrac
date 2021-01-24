cmake_minimum_required(VERSION 3.5)

set(PROJECT_NAME partrac)
project(partrac)

# Set CMake behavior
cmake_policy(SET CMP0004 NEW)

set(CMAKE_BUILD_TYPE "Developer" CACHE STRING "Choose build type, options are: Developer" FORCE)

set(CMAKE_CXX_COMPILER g++-10)
set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_FLAGS "-pedantic-errors -Wall -Wextra -Werror -O3 -ffast-math -fconcepts -lstdc++")
#set(CMAKE_CXX_FLAGS "-O3 -std=c++17 -ffast-math -L/usr/lib -Iinclude/ -lstdc++ -lm ${CMAKE_CXX_FLAGS}")
#set(CMAKE_CXX_FLAGS "-Wextra -O3 -fconcepts -lstdc++ -lm -std=c++17")
#set(CMAKE_CXX_FLAGS "-pedantic-errors -Wextra -O3 -ggdb3 -fconcepts -lstdc++ -lm -std=c++17 -march=native")
set(CMAKE_CXX_FLAGS "-lm -ggdb3 -fsanitize=address")

#set(CMAKE_EXE_LINKER_FLAGS "/usr/lib/x86_64-linux-gnu/libprofiler.so.0")
set(CMAKE_EXE_LINKER_FLAGS "-fsanitize=address") #-static-libasan

#find_package(coz-profiler)
add_definitions(-DDEBUG)



#set(ENABLE_DOLFIN "true" CACHE STRING "Set to true to enable Dolfin import.")
include(FeatureSummary)

set(OPTIONAL_PACKAGES "")
list(APPEND OPTIONAL_PACKAGES "Dolfin")
foreach (OPTIONAL_PACKAGE ${OPTIONAL_PACKAGES})
  string (TOUPPER "PARTRAC_ENABLE_${OPTIONAL_PACKAGE}" OPTION_NAME)
  option(${OPTION_NAME} "Compile with support for ${OPTIONAL_PACKAGE}." ON)
  add_feature_info(${OPTION_NAME} ${OPTION_NAME} "Compile with support for ${OPTIONAL_PACKAGE}.")
endforeach()

if(PARTRAC_ENABLE_DOLFIN)
  message(STATUS "*** with dolfin support ***")
  add_definitions(-DUSE_DOLFIN)
  #include_directories(${DOLFIN_SOURCE_DIR} ${DOLFIN_SOURCE_DIR}/dolfin ${DOLFIN_BINARY_DIR})
  find_package(DOLFIN)
  include(${DOLFIN_USE_FILE})
endif()

find_package(HDF5 REQUIRED COMPONENTS CXX HL)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(${HDF5_INCLUDE_DIR})

aux_source_directory(src/ SRC_LIST)

add_executable(${PROJECT_NAME} ${SRC_LIST})

if (PARTRAC_ENABLE_DOLFIN)
  target_link_libraries(${PROJECT_NAME} dolfin)
endif()
target_link_libraries(${PROJECT_NAME} ${HDF5_LIBRARIES}
  ${HDF5_HL_LIBRARIES} Eigen3::Eigen)
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)

install(TARGETS ${PROJECT_NAME} DESTINATION /usr/bin/)
