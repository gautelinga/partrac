# One library per directory; each exports only its own headers, so a header
# can include only what its library depends on.

# Parameter schema layer. No dependencies, so the unit tests link it on its own.
add_library(${PROJECT_NAME}_params
  params/Params.cpp
  params/param_print.cpp
  )
target_include_directories(${PROJECT_NAME}_params PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/params)
target_compile_features(${PROJECT_NAME}_params PUBLIC cxx_std_17)
target_compile_options(${PROJECT_NAME}_params PRIVATE -Wall -Wextra)

# Types, geometry, strings, stamps, files, generators, run folders
add_library(${PROJECT_NAME}_base
  base/Timestamps.cpp
  base/files.cpp
  )
target_include_directories(${PROJECT_NAME}_base PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/base)
target_link_libraries(${PROJECT_NAME}_base PUBLIC ${PROJECT_NAME}_params
  Eigen3::Eigen OpenMP::OpenMP_CXX stdc++fs)
target_compile_features(${PROJECT_NAME}_base PUBLIC cxx_std_17)

# HDF5 and text input and output
add_library(${PROJECT_NAME}_io
  io/io.cpp
  )
target_include_directories(${PROJECT_NAME}_io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/io
  ${HDF5_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME}_io PUBLIC ${PROJECT_NAME}_base
  ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

# Interpolators; the mesh ones need dolfin
set(PARTRAC_INTERPOL_SOURCES interpol/interpol_factory.cpp)
if (PARTRAC_ENABLE_DOLFIN)
  list(APPEND PARTRAC_INTERPOL_SOURCES
    interpol/DolfInterpol.cpp
    interpol/TetInterpol.cpp
    interpol/TriangleInterpol.cpp
    interpol/TriangleFreqInterpol.cpp
    interpol/XDMFTriangleInterpol.cpp
    interpol/XDMFTetInterpol.cpp
    interpol/Tet.cpp
    interpol/Triangle.cpp
    interpol/dolfin_helpers.cpp
    interpol/xdmf_helpers.cpp
    )
endif()
add_library(${PROJECT_NAME}_interpol ${PARTRAC_INTERPOL_SOURCES})
target_include_directories(${PROJECT_NAME}_interpol PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/interpol)
target_link_libraries(${PROJECT_NAME}_interpol PUBLIC ${PROJECT_NAME}_io ${PROJECT_NAME}_base
  OpenMP::OpenMP_CXX)
if (PARTRAC_ENABLE_DOLFIN)
  target_link_libraries(${PROJECT_NAME}_interpol PUBLIC dolfin)
  # dolfin.h is most of a loader's compile time
  if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.16)
    target_precompile_headers(${PROJECT_NAME}_interpol PRIVATE <dolfin.h> <Eigen/Dense>)
  endif()
endif()

# The step loops, compiled once per interpolator; the same list as with_concrete
# (interpol_dispatch.hpp) and set_interpolate_mode (interpol_factory.cpp)
set(PARTRAC_STEP_INTERPOLATORS AnalyticInterpol StructuredInterpol)
if (PARTRAC_ENABLE_DOLFIN)
  list(APPEND PARTRAC_STEP_INTERPOLATORS DolfInterpol TetInterpol TriangleInterpol
    TriangleFreqInterpol XDMFTriangleInterpol XDMFTetInterpol)
endif()
set(PARTRAC_STEP_SOURCES)
foreach(INTERP ${PARTRAC_STEP_INTERPOLATORS})
  configure_file(particles/steps.cpp.in ${CMAKE_CURRENT_BINARY_DIR}/steps_${INTERP}.cpp @ONLY)
  list(APPEND PARTRAC_STEP_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/steps_${INTERP}.cpp)
endforeach()

# Particles, meshes, initializers, steps and the run loop
add_library(${PROJECT_NAME}_particles
  particles/mesh.cpp
  particles/Topology.cpp
  particles/Initializer.cpp
  particles/stepping.cpp
  ${PARTRAC_STEP_SOURCES}
  )
target_include_directories(${PROJECT_NAME}_particles PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/particles)
target_link_libraries(${PROJECT_NAME}_particles PUBLIC ${PROJECT_NAME}_interpol)
if (PARTRAC_ENABLE_DOLFIN)
  # the step sources include the mesh loaders
  if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.16)
    target_precompile_headers(${PROJECT_NAME}_particles PRIVATE <dolfin.h> <Eigen/Dense>)
  endif()
endif()

# What every app links against
add_library(${PROJECT_NAME}_core INTERFACE)
target_link_libraries(${PROJECT_NAME}_core INTERFACE ${PROJECT_NAME}_particles)
if (PARTRAC_ENABLE_MPI)
  target_link_libraries(${PROJECT_NAME}_core INTERFACE ${MPI_LIBRARIES})
endif()
