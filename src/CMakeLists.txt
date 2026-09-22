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

# Interpolators; only the reference loader (mode fenics) needs dolfin
set(PARTRAC_INTERPOL_SOURCES interpol/interpol_factory.cpp interpol/mesh_tables.cpp)
# Shewchuk's exact predicates, vendored verbatim: C, none of our warnings, and
# no floating-point contraction, which their exactness rests on. An object
# library of its own so that a C source takes none of the C++ options.
add_library(${PROJECT_NAME}_predicates OBJECT interpol/predicates/predicates.c)
target_compile_options(${PROJECT_NAME}_predicates PRIVATE -w -ffp-contract=off)
set_property(TARGET ${PROJECT_NAME}_predicates PROPERTY POSITION_INDEPENDENT_CODE ON)

# The cell tree and its Morton order need no dolfin, nor do the mesh core and
# the loaders that read a dolfin HDF5 checkpoint or an XDMF case from its arrays
list(APPEND PARTRAC_INTERPOL_SOURCES
  interpol/cell_tree.cpp
  interpol/MeshCore.cpp
  interpol/SimplexInterpol.cpp
  interpol/SimplexInterpol_eval.cpp
  interpol/SplitInterpol.cpp
  interpol/SplitInterpol_eval.cpp
  interpol/XDMFInterpol_eval.cpp
  interpol/simplex_load.cpp
  interpol/SimplexFreqInterpol.cpp
  interpol/XDMFInterpol.cpp
  interpol/near_wall.cpp
  interpol/xdmf_helpers.cpp
  interpol/phase_timing.cpp
  $<TARGET_OBJECTS:${PROJECT_NAME}_predicates>
  )

# The reference loader (mode fenics) and the dolfin::Cell constructors it
# builds its cells with: an object library of its own, the only one that
# includes <dolfin.h>, so the dolfin-free loaders never see it
if (PARTRAC_ENABLE_DOLFIN)
  add_library(${PROJECT_NAME}_dolfin OBJECT
    interpol/DolfInterpol.cpp
    interpol/dolfin_ref.cpp
    interpol/Tet.cpp
    interpol/Triangle.cpp
    )
  target_include_directories(${PROJECT_NAME}_dolfin PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/interpol)
  target_link_libraries(${PROJECT_NAME}_dolfin PUBLIC ${PROJECT_NAME}_io ${PROJECT_NAME}_base
    OpenMP::OpenMP_CXX dolfin)
  set_property(TARGET ${PROJECT_NAME}_dolfin PROPERTY POSITION_INDEPENDENT_CODE ON)
  list(APPEND PARTRAC_INTERPOL_SOURCES $<TARGET_OBJECTS:${PROJECT_NAME}_dolfin>)
endif()
add_library(${PROJECT_NAME}_interpol ${PARTRAC_INTERPOL_SOURCES})
target_include_directories(${PROJECT_NAME}_interpol PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/interpol)
target_link_libraries(${PROJECT_NAME}_interpol PUBLIC ${PROJECT_NAME}_io ${PROJECT_NAME}_base
  OpenMP::OpenMP_CXX)
if (PARTRAC_ENABLE_DOLFIN)
  target_link_libraries(${PROJECT_NAME}_interpol PUBLIC dolfin)
endif()

# The step loops, compiled once per interpolator; the same list as with_concrete
# (interpol_dispatch.hpp) and set_interpolate_mode (interpol_factory.cpp)
set(PARTRAC_STEP_INTERPOLATORS AnalyticInterpol StructuredInterpol StructuredConstInterpol
  TetInterpol TriangleInterpol TriangleFreqInterpol TetFreqInterpol
  XDMFTriangleInterpol XDMFTetInterpol SplitTriangleInterpol SplitTetInterpol)
set(PARTRAC_DOLFIN_STEP_INTERPOLATORS)
if (PARTRAC_ENABLE_DOLFIN)
  set(PARTRAC_DOLFIN_STEP_INTERPOLATORS DolfTriangleInterpol DolfTetInterpol)
endif()
set(PARTRAC_STEP_SOURCES)
set(PARTRAC_DOLFIN_STEP_SOURCES)
foreach(INTERP ${PARTRAC_STEP_INTERPOLATORS} ${PARTRAC_DOLFIN_STEP_INTERPOLATORS})
  set(STEP_LIST PARTRAC_STEP_SOURCES)
  if (INTERP IN_LIST PARTRAC_DOLFIN_STEP_INTERPOLATORS)
    set(STEP_LIST PARTRAC_DOLFIN_STEP_SOURCES)
  endif()
  configure_file(particles/steps.cpp.in ${CMAKE_CURRENT_BINARY_DIR}/steps_${INTERP}.cpp @ONLY)
  list(APPEND ${STEP_LIST} ${CMAKE_CURRENT_BINARY_DIR}/steps_${INTERP}.cpp)
  foreach(ELEMENT Point Vector Tensor)
    configure_file(particles/steps_rk4.cpp.in ${CMAKE_CURRENT_BINARY_DIR}/steps_rk4_${ELEMENT}_${INTERP}.cpp @ONLY)
    list(APPEND ${STEP_LIST} ${CMAKE_CURRENT_BINARY_DIR}/steps_rk4_${ELEMENT}_${INTERP}.cpp)
  endforeach()
endforeach()

# The step loops of the reference loader, for the same reason
if (PARTRAC_ENABLE_DOLFIN)
  add_library(${PROJECT_NAME}_dolfin_steps OBJECT ${PARTRAC_DOLFIN_STEP_SOURCES})
  target_include_directories(${PROJECT_NAME}_dolfin_steps PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/particles)
  target_link_libraries(${PROJECT_NAME}_dolfin_steps PUBLIC ${PROJECT_NAME}_interpol)
  set_property(TARGET ${PROJECT_NAME}_dolfin_steps PROPERTY POSITION_INDEPENDENT_CODE ON)
  list(APPEND PARTRAC_STEP_SOURCES $<TARGET_OBJECTS:${PROJECT_NAME}_dolfin_steps>)
endif()

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

# What every app links against
add_library(${PROJECT_NAME}_core INTERFACE)
target_link_libraries(${PROJECT_NAME}_core INTERFACE ${PROJECT_NAME}_particles)
if (PARTRAC_ENABLE_MPI)
  target_link_libraries(${PROJECT_NAME}_core INTERFACE ${MPI_LIBRARIES})
endif()
