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
include(GNUInstallDirs)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(PartracApp)
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
list(APPEND OPTIONAL_PACKAGES "Dolfin")
foreach (OPTIONAL_PACKAGE ${OPTIONAL_PACKAGES})
  string (TOUPPER "PARTRAC_ENABLE_${OPTIONAL_PACKAGE}" OPTION_NAME)
  option(${OPTION_NAME} "Compile with support for ${OPTIONAL_PACKAGE}." ON)
  add_feature_info(${OPTION_NAME} ${OPTION_NAME} "Compile with support for ${OPTIONAL_PACKAGE}.")
endforeach()

if(PARTRAC_ENABLE_DOLFIN)
  message(STATUS "*** with dolfin support ***")
  add_definitions(-DUSE_DOLFIN)
  find_package(DOLFIN REQUIRED)
  include(${DOLFIN_USE_FILE})
endif()

message(STATUS "****** ${CMAKE_BUILD_TYPE}")
message(STATUS "dolfin cxx flags: ${DOLFIN_CXX_FLAGS}")
message(STATUS "cxx flags: ${CMAKE_CXX_FLAGS}")

# Tuning for the build machine. Off by default: such a binary crashes with
# SIGILL on a machine with an older instruction set.
option(PARTRAC_NATIVE_ARCH "Optimize for the CPU of the build machine." OFF)
if (PARTRAC_NATIVE_ARCH)
  set(PARTRAC_ARCH_FLAGS "-march=native")
endif()

# Warnings belong to every build type: a bug only Release can see is still a
# bug. Unused parameters are off: the interfaces take arguments some
# implementations ignore.
# -fno-permissive comes last on purpose: DOLFIN_CXX_FLAGS carries -fpermissive,
# which would otherwise downgrade errors across the whole tree.
set(PARTRAC_WARNING_FLAGS "-Wall -Wextra -Wno-unused-parameter -fno-permissive")

if( CMAKE_BUILD_TYPE MATCHES "Release")
  # Optimized build
  message(STATUS "Release build")
  set(CMAKE_CXX_FLAGS "${DOLFIN_CXX_FLAGS} ${CMAKE_CXX_FLAGS} ${PARTRAC_WARNING_FLAGS} -DNDEBUG -O3 ${PARTRAC_ARCH_FLAGS}")
else() # (CMAKE_BUILD_TYPE STREQUAL "Debug")
  # # Debug build
  message(STATUS "Debug build")
  set(CMAKE_CXX_FLAGS "${DOLFIN_CXX_FLAGS} ${CMAKE_CXX_FLAGS} ${PARTRAC_WARNING_FLAGS} -DDEBUG -O2 ${PARTRAC_ARCH_FLAGS} -ggdb3 -pedantic-errors -pg")
endif()

message(STATUS "cxx flags: ${CMAKE_CXX_FLAGS}")

# Profiling: Use as LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libprofiler.so.0 CPUPROFILE=prof.txt ./partrac ...
# pprof-symbolize --ps ./partrac prof.txt > prof.ps
#list(APPEND DOLFIN_LIBRARIES "/usr/lib/x86_64-linux-gnu/libprofiler.so.0")

# Profiling build
#set(CMAKE_CXX_FLAGS "${DOLFIN_CXX_FLAGS} ${CMAKE_CXX_FLAGS} -DNDEBUG -O3 -march=native -ggdb3 -pg -pedantic-errors -Wextra")

find_package(HDF5 REQUIRED COMPONENTS CXX HL)
# A parallel HDF5 includes mpi.h from its public headers, so it has to be on the
# include path even though partrac itself is serial. This is the only HDF5 some
# distributions package next to dolfin.
if (HDF5_IS_PARALLEL)
  message(STATUS "HDF5 is parallel, adding MPI for its headers")
  find_package(MPI REQUIRED COMPONENTS CXX)
endif()
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)
# Headers only, and used whether or not dolfin is enabled
find_package(Boost REQUIRED)

include_directories(${HDF5_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_DIRS})

add_subdirectory(src)

# Each app installs itself, see partrac_add_app
add_subdirectory(apps)
# Build features for the tests
if(PARTRAC_ENABLE_DOLFIN)
  file(WRITE "${PROJECT_BINARY_DIR}/bin/build_features.txt" "dolfin=on\n")
else()
  file(WRITE "${PROJECT_BINARY_DIR}/bin/build_features.txt" "dolfin=off\n")
endif()

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
