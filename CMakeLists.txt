cmake_minimum_required(VERSION 3.5)

set(PROJECT_NAME partrac)
project(${PROJECT_NAME})

# Set CMake behavior
cmake_policy(SET CMP0004 NEW)

set(CMAKE_CXX_COMPILER g++-10)
set(CMAKE_CXX_STANDARD 17)

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
  find_package(DOLFIN)
  include(${DOLFIN_USE_FILE})
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)
endif()

# # Debug build
# set(CMAKE_CXX_FLAGS "${DOLFIN_CXX_FLAGS} ${CMAKE_CXX_FLAGS} -DDEBUG -O3 -march=native -ggdb3 -pedantic-errors -Wextra")
# set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)

# Optimized build
set(CMAKE_CXX_FLAGS "${DOLFIN_CXX_FLAGS} ${CMAKE_CXX_FLAGS} -DNDEBUG -O3 -march=native -ggdb3 -pedantic-errors -Wextra")

# Profiling: Use as LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libprofiler.so.0 CPUPROFILE=prof.txt ./partrac ...
# pprof-symbolize --ps ./partrac prof.txt > prof.ps
#list(APPEND DOLFIN_LIBRARIES "/usr/lib/x86_64-linux-gnu/libprofiler.so.0")


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





# # Get DOLFIN configuration data (DOLFINConfig.cmake must be in
# # DOLFIN_CMAKE_CONFIG_PATH)
# find_package(DOLFIN REQUIRED)
# add_definitions(-DUSE_DOLFIN)

# include(${DOLFIN_USE_FILE})
# set(CMAKE_CXX_STANDARD 17)

# # # Debug build
# # set(CMAKE_CXX_FLAGS "${DOLFIN_CXX_FLAGS} ${CMAKE_CXX_FLAGS} -DDEBUG -O3 -march=native -ggdb3 -pedantic-errors -Wextra")
# # set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)

# # Profiling: Use as LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libprofiler.so.0 CPUPROFILE=prof.txt ./partrac ...
# # pprof-symbolize --ps ./partrac prof.txt > prof.ps
# #list(APPEND DOLFIN_LIBRARIES "/usr/lib/x86_64-linux-gnu/libprofiler.so.0")

# # Optimized build
# set(CMAKE_CXX_FLAGS "${DOLFIN_CXX_FLAGS} ${CMAKE_CXX_FLAGS} -DNDEBUG -O3 -march=native -ggdb3 -pedantic-errors -Wextra")
# set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)

# find_package(HDF5 REQUIRED COMPONENTS CXX HL)
# find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# include_directories(${HDF5_INCLUDE_DIR})

# aux_source_directory(src/ SRC_LIST)
# add_executable(${PROJECT_NAME} ${SRC_LIST})

# target_link_libraries(${PROJECT_NAME} ${HDF5_LIBRARIES}
#   ${HDF5_HL_LIBRARIES} Eigen3::Eigen ${DOLFIN_LIBRARIES})
