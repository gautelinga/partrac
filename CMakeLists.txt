cmake_minimum_required(VERSION 3.14)

set(PROJECT_NAME partrac)
project(${PROJECT_NAME}
  VERSION 0.1
  DESCRIPTION "Particle tracker"
  LANGUAGES C CXX
  )

# Set CMake behavior
cmake_policy(SET CMP0004 NEW)

#set(CMAKE_CXX_COMPILER g++-10)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

include(CTest)
find_package(Doxygen)
if (Doxygen_FOUND)
  add_subdirectory(docs)
else()
  message(STATUS "Doxygen not found, not building docs.")
endif()

#set(ENABLE_DOLFIN "true" CACHE STRING "Set to true to enable Dolfin import.")
include(FeatureSummary)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)
endif()

set(OPTIONAL_PACKAGES "")
list(APPEND OPTIONAL_PACKAGES "Dolfin" "mpi")
foreach (OPTIONAL_PACKAGE ${OPTIONAL_PACKAGES})
  string (TOUPPER "PARTRAC_ENABLE_${OPTIONAL_PACKAGE}" OPTION_NAME)
  option(${OPTION_NAME} "Compile with support for ${OPTIONAL_PACKAGE}." ON)
  add_feature_info(${OPTION_NAME} ${OPTION_NAME} "Compile with support for ${OPTIONAL_PACKAGE}.")
endforeach()

if(PARTRAC_ENABLE_DOLFIN)
  message(STATUS "*** with dolfin support ***")
  add_definitions(-DUSE_DOLFIN)
  find_package(DOLFIN)
  include(${DOLFIN_USE_FILE})
endif()

if (PARTRAC_ENABLE_MPI)
  message(STATUS "*** with MPI support ***")
  find_package(MPI)
  include_directories(${MPI_INCLUDE_PATH})
endif()

message(STATUS "****** ${CMAKE_BUILD_TYPE}")
message(STATUS "dolfin cxx flags: ${DOLFIN_CXX_FLAGS}")
message(STATUS "cxx flags: ${CMAKE_CXX_FLAGS}")

if( CMAKE_BUILD_TYPE MATCHES "Release")
  # Optimized build
  message(STATUS "Release build")
  set(CMAKE_CXX_FLAGS "${DOLFIN_CXX_FLAGS} ${CMAKE_CXX_FLAGS} -DNDEBUG -O3 -march=native")
else() # (CMAKE_BUILD_TYPE STREQUAL "Debug")
  # # Debug build
  message(STATUS "Debug build")
  set(CMAKE_CXX_FLAGS "${DOLFIN_CXX_FLAGS} ${CMAKE_CXX_FLAGS} -DDEBUG -O2 -march=native -ggdb3 -pedantic-errors -Wextra -pg")
endif()

message(STATUS "cxx flags: ${CMAKE_CXX_FLAGS}")

# Profiling: Use as LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libprofiler.so.0 CPUPROFILE=prof.txt ./partrac ...
# pprof-symbolize --ps ./partrac prof.txt > prof.ps
#list(APPEND DOLFIN_LIBRARIES "/usr/lib/x86_64-linux-gnu/libprofiler.so.0")

# Profiling build
#set(CMAKE_CXX_FLAGS "${DOLFIN_CXX_FLAGS} ${CMAKE_CXX_FLAGS} -DNDEBUG -O3 -march=native -ggdb3 -pg -pedantic-errors -Wextra")

find_package(HDF5 REQUIRED COMPONENTS CXX HL)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(${HDF5_INCLUDE_DIR})

add_subdirectory(src)

set(PROJECT_TARGETS ${PROJECT_NAME})
#aux_source_directory(src SRC_LIST)
add_subdirectory(apps)

message(STATUS ${PROJECT_TARGETS})

separate_arguments(PROJECT_TARGETS)
install(TARGETS ${PROJECT_TARGETS} DESTINATION /usr/bin/)

# FetchContent added in CMake 3.11, downloads during the configure step

if (BUILD_TESTING)
  include(FetchContent)
  enable_testing()
  add_subdirectory(tests)
endif()


# # # Debug build
# # set(CMAKE_CXX_FLAGS "${DOLFIN_CXX_FLAGS} ${CMAKE_CXX_FLAGS} -DDEBUG -O3 -march=native -ggdb3 -pedantic-errors -Wextra")
# # set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)

# # Profiling: Use as LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libprofiler.so.0 CPUPROFILE=prof.txt ./partrac ...
# # pprof-symbolize --ps ./partrac prof.txt > prof.ps
# #list(APPEND DOLFIN_LIBRARIES "/usr/lib/x86_64-linux-gnu/libprofiler.so.0")

# # Optimized build
# set(CMAKE_CXX_FLAGS "${DOLFIN_CXX_FLAGS} ${CMAKE_CXX_FLAGS} -DNDEBUG -O3 -march=native -ggdb3 -pedantic-errors -Wextra")
# set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)

message(STATUS "cxx flags: ${CMAKE_CXX_FLAGS}")
