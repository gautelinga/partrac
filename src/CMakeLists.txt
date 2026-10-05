# One library per directory; each exports only its own headers, so a header
# can include only what its library depends on.

# Parameter schema layer. No dependencies, so the unit tests link it on its own.
add_library(${PROJECT_NAME}_params
  params/Params.cpp
  params/param_print.cpp
  )
target_include_directories(${PROJECT_NAME}_params PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/params)
target_compile_features(${PROJECT_NAME}_params PUBLIC cxx_std_17)

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
# A coverage build's one OpenMP unit: it splits its work by omp_get_max_threads
if (PARTRAC_COVERAGE)
  separate_arguments(PARTRAC_OPENMP_FLAGS UNIX_COMMAND "${OpenMP_CXX_FLAGS}")
  set_source_files_properties(interpol/mesh_tables.cpp PROPERTIES COMPILE_OPTIONS "${PARTRAC_OPENMP_FLAGS}")
endif()
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
  interpol/SimplexInterpol_held.cpp
  interpol/SplitInterpol.cpp
  interpol/SplitInterpol_eval.cpp
  interpol/SplitInterpol_held.cpp
  interpol/XDMFInterpol_eval.cpp
  interpol/XDMFInterpol_held.cpp
  interpol/simplex_load.cpp
  interpol/SimplexFreqInterpol.cpp
  interpol/XDMFInterpol.cpp
  interpol/near_wall.cpp
  interpol/xdmf_helpers.cpp
  interpol/phase_timing.cpp
  interpol/openfoam_host.cpp
  interpol/openfoam_split.cpp
  interpol/openfoam_nodes.cpp
  interpol/openfoam_phase.cpp
  interpol/OpenFoamInterpol.cpp
  interpol/OpenFoamInterpol_eval.cpp
  interpol/OpenFoamInterpol_held.cpp
  )
set(PARTRAC_INTERPOL_OBJECTS $<TARGET_OBJECTS:${PROJECT_NAME}_predicates>)

# The OpenFOAM reader, the only unit that includes an OpenFOAM header: its
# include tree as system headers and its flags, linked with DT_RPATH to the
# dummy Pstream, since OpenFOAM's libraries carry no RUNPATH of their own
if (PARTRAC_ENABLE_OPENFOAM)
  add_library(${PROJECT_NAME}_openfoam SHARED interpol/openfoam_load.cpp)
  target_include_directories(${PROJECT_NAME}_openfoam PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/interpol)
  target_include_directories(${PROJECT_NAME}_openfoam SYSTEM PRIVATE
    ${FOAM_SRC}/OpenFOAM/lnInclude ${FOAM_SRC}/OSspecific/POSIX/lnInclude
    ${FOAM_SRC}/finiteVolume/lnInclude ${FOAM_SRC}/meshTools/lnInclude)
  target_compile_definitions(${PROJECT_NAME}_openfoam PRIVATE ${FOAM_ARCH}
    WM_ARCH_OPTION=${FOAM_ARCH_OPTION} WM_${FOAM_PRECISION} WM_LABEL_SIZE=${FOAM_LABEL_SIZE} NoRepository)
  target_compile_features(${PROJECT_NAME}_openfoam PRIVATE cxx_std_17)
  target_compile_options(${PROJECT_NAME}_openfoam PRIVATE -ftemplate-depth-100)
  target_link_libraries(${PROJECT_NAME}_openfoam PRIVATE ${FOAM_LIBBIN}/libfiniteVolume.so
    ${FOAM_LIBBIN}/libOpenFOAM.so ${FOAM_LIBBIN}/dummy/libPstream.so)
  target_link_options(${PROJECT_NAME}_openfoam PRIVATE -Wl,--disable-new-dtags -Wl,--no-as-needed)
  set_target_properties(${PROJECT_NAME}_openfoam PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
    BUILD_RPATH "${FOAM_LIBBIN}/dummy;${FOAM_LIBBIN}"
    INSTALL_RPATH "${FOAM_LIBBIN}/dummy;${FOAM_LIBBIN}")
  install(TARGETS ${PROJECT_NAME}_openfoam LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
  # The build tree's plugin, and the installed one's place from an installed app
  file(RELATIVE_PATH PARTRAC_OPENFOAM_LIBDIR_REL ${CMAKE_INSTALL_FULL_BINDIR} ${CMAKE_INSTALL_FULL_LIBDIR})
  set_source_files_properties(interpol/openfoam_host.cpp PROPERTIES
    COMPILE_DEFINITIONS "PARTRAC_OPENFOAM_PLUGIN=\"$<TARGET_FILE:${PROJECT_NAME}_openfoam>\";PARTRAC_OPENFOAM_LIBDIR_REL=\"${PARTRAC_OPENFOAM_LIBDIR_REL}\";PARTRAC_WM_PROJECT_DIR=\"${FOAM_PROJECT_DIR}\"")
endif()

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
  list(APPEND PARTRAC_INTERPOL_OBJECTS $<TARGET_OBJECTS:${PROJECT_NAME}_dolfin>)
endif()
# The units an object library, which compiles while the other object libraries
# do; the archive takes all of them
add_library(${PROJECT_NAME}_interpol_units OBJECT ${PARTRAC_INTERPOL_SOURCES})
add_library(${PROJECT_NAME}_interpol $<TARGET_OBJECTS:${PROJECT_NAME}_interpol_units>
  ${PARTRAC_INTERPOL_OBJECTS})
foreach(t ${PROJECT_NAME}_interpol_units ${PROJECT_NAME}_interpol)
  target_include_directories(${t} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/interpol
    ${CMAKE_CURRENT_BINARY_DIR}/generated)
  target_link_libraries(${t} PUBLIC ${PROJECT_NAME}_io ${PROJECT_NAME}_base OpenMP::OpenMP_CXX)
  if (PARTRAC_ENABLE_DOLFIN)
    target_link_libraries(${t} PUBLIC dolfin)
  endif()
  if (PARTRAC_ENABLE_OPENFOAM)
    target_link_libraries(${t} PUBLIC ${CMAKE_DL_LIBS})
  endif()
endforeach()
if (PARTRAC_ENABLE_OPENFOAM)
  add_dependencies(${PROJECT_NAME}_interpol ${PROJECT_NAME}_openfoam)
endif()

# The interpolators, name|header|type, one row each, a loader needing dolfin in
# the second table: the step loops are compiled per row, and interpolators.hpp,
# written here, includes the headers, names the types and lists them for
# with_concrete; the order is the cast chain's and the step objects' in the
# archive. A new loader is a row here, its mode in interpol_modes() and
# set_interpolate_mode, and, without regions, an entry in the RK4cells list below
set(PARTRAC_INTERPOLATOR_TABLE
  "AnalyticInterpol|AnalyticInterpol.hpp|AnalyticInterpol"
  "StructuredInterpol|StructuredInterpol.hpp|StructuredInterpol"
  "StructuredConstInterpol|StructuredInterpol.hpp|StructuredConstInterpol"
  "TriangleInterpol|SimplexInterpol.hpp|SimplexInterpol<Triangle>"
  "TetInterpol|SimplexInterpol.hpp|SimplexInterpol<Tet>"
  "TriangleFreqInterpol|SimplexFreqInterpol.hpp|SimplexFreqInterpol<Triangle>"
  "TetFreqInterpol|SimplexFreqInterpol.hpp|SimplexFreqInterpol<Tet>"
  "XDMFTriangleInterpol|XDMFInterpol.hpp|XDMFInterpol<Triangle>"
  "XDMFTetInterpol|XDMFInterpol.hpp|XDMFInterpol<Tet>"
  "SplitTriangleInterpol|SplitInterpol.hpp|SplitInterpol<Triangle>"
  "SplitTetInterpol|SplitInterpol.hpp|SplitInterpol<Tet>"
  "OpenFoamTriangleInterpol|OpenFoamInterpol.hpp|OpenFoamInterpol<Triangle>"
  "OpenFoamTetInterpol|OpenFoamInterpol.hpp|OpenFoamInterpol<Tet>")
set(PARTRAC_DOLFIN_INTERPOLATOR_TABLE)
if (PARTRAC_ENABLE_DOLFIN)
  set(PARTRAC_DOLFIN_INTERPOLATOR_TABLE
    "DolfTriangleInterpol|DolfInterpol.hpp|DolfInterpol<Triangle>"
    "DolfTetInterpol|DolfInterpol.hpp|DolfInterpol<Tet>")
endif()
set(PARTRAC_STEP_INTERPOLATORS)
set(PARTRAC_DOLFIN_STEP_INTERPOLATORS)
set(INTERP_INCLUDES "")
set(INTERP_ALIASES "")
set(INTERP_LIST "")
foreach(ROW ${PARTRAC_INTERPOLATOR_TABLE} ${PARTRAC_DOLFIN_INTERPOLATOR_TABLE})
  string(REPLACE "|" ";" FIELDS "${ROW}")
  list(GET FIELDS 0 NAME)
  list(GET FIELDS 1 HEADER)
  list(GET FIELDS 2 TYPE)
  if (ROW IN_LIST PARTRAC_DOLFIN_INTERPOLATOR_TABLE)
    list(APPEND PARTRAC_DOLFIN_STEP_INTERPOLATORS ${NAME})
  else()
    list(APPEND PARTRAC_STEP_INTERPOLATORS ${NAME})
  endif()
  set(INTERP_HEADER_${NAME} ${HEADER})
  set(INTERP_ALIAS_${NAME} "")
  if (NOT TYPE STREQUAL NAME)
    set(INTERP_ALIAS_${NAME} "using ${NAME} = ${TYPE};")
    string(APPEND INTERP_ALIASES "using ${NAME} = ${TYPE};\n")
  endif()
  if (NOT INTERP_INCLUDES MATCHES "\"${HEADER}\"")
    string(APPEND INTERP_INCLUDES "#include \"${HEADER}\"\n")
  endif()
  string(APPEND INTERP_LIST " X(${NAME})")
endforeach()
configure_file(interpol/interpolators.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/generated/interpolators.hpp @ONLY)
# The interpolators with regions (has_regions), which RK4cells steps in cells
set(PARTRAC_CELLS_INTERPOLATORS ${PARTRAC_STEP_INTERPOLATORS})
list(REMOVE_ITEM PARTRAC_CELLS_INTERPOLATORS AnalyticInterpol StructuredConstInterpol)
set(PARTRAC_STEP_SOURCES)
set(PARTRAC_CELLS_STEP_SOURCES)
set(PARTRAC_DOLFIN_STEP_SOURCES)
# The step objects in the order of this loop, where CMake writes them for the
# object libraries below
set(PARTRAC_STEP_OBJECTS)
set(PARTRAC_OBJDIR ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles)
foreach(INTERP ${PARTRAC_STEP_INTERPOLATORS} ${PARTRAC_DOLFIN_STEP_INTERPOLATORS})
  set(STEP_LIST PARTRAC_STEP_SOURCES)
  set(STEP_DIR ${PARTRAC_OBJDIR}/${PROJECT_NAME}_steps_units.dir)
  if (INTERP IN_LIST PARTRAC_DOLFIN_STEP_INTERPOLATORS)
    set(STEP_LIST PARTRAC_DOLFIN_STEP_SOURCES)
    set(STEP_DIR ${PARTRAC_OBJDIR}/${PROJECT_NAME}_dolfin_steps.dir)
  endif()
  set(INTERP_HEADER "${INTERP_HEADER_${INTERP}}")
  set(INTERP_ALIAS "${INTERP_ALIAS_${INTERP}}")
  configure_file(particles/steps.cpp.in ${CMAKE_CURRENT_BINARY_DIR}/steps_${INTERP}.cpp @ONLY)
  list(APPEND ${STEP_LIST} ${CMAKE_CURRENT_BINARY_DIR}/steps_${INTERP}.cpp)
  list(APPEND PARTRAC_STEP_OBJECTS ${STEP_DIR}/steps_${INTERP}.cpp${CMAKE_CXX_OUTPUT_EXTENSION})
  foreach(ELEMENT Point Vector Tensor)
    configure_file(particles/steps_rk4.cpp.in ${CMAKE_CURRENT_BINARY_DIR}/steps_rk4_${ELEMENT}_${INTERP}.cpp @ONLY)
    list(APPEND ${STEP_LIST} ${CMAKE_CURRENT_BINARY_DIR}/steps_rk4_${ELEMENT}_${INTERP}.cpp)
    list(APPEND PARTRAC_STEP_OBJECTS ${STEP_DIR}/steps_rk4_${ELEMENT}_${INTERP}.cpp${CMAKE_CXX_OUTPUT_EXTENSION})
    if (INTERP IN_LIST PARTRAC_CELLS_INTERPOLATORS)
      configure_file(particles/steps_cells.cpp.in ${CMAKE_CURRENT_BINARY_DIR}/steps_cells_${ELEMENT}_${INTERP}.cpp @ONLY)
      list(APPEND PARTRAC_CELLS_STEP_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/steps_cells_${ELEMENT}_${INTERP}.cpp)
      list(APPEND PARTRAC_STEP_OBJECTS
        ${PARTRAC_OBJDIR}/${PROJECT_NAME}_cells_units.dir/steps_cells_${ELEMENT}_${INTERP}.cpp${CMAKE_CXX_OUTPUT_EXTENSION})
    endif()
  endforeach()
endforeach()

# Precompiled step headers; no prologue, which would silence their warnings
set(CMAKE_PCH_PROLOGUE "")
add_library(${PROJECT_NAME}_steps_units OBJECT ${PARTRAC_STEP_SOURCES})
target_precompile_headers(${PROJECT_NAME}_steps_units PRIVATE particles/steps_impl.hpp)
add_library(${PROJECT_NAME}_cells_units OBJECT ${PARTRAC_CELLS_STEP_SOURCES})
target_precompile_headers(${PROJECT_NAME}_cells_units PRIVATE particles/steps_cells_impl.hpp)
set(PARTRAC_STEP_TARGETS ${PROJECT_NAME}_steps_units ${PROJECT_NAME}_cells_units)

# The step loops of the reference loader; no precompiled header: too large
if (PARTRAC_ENABLE_DOLFIN)
  add_library(${PROJECT_NAME}_dolfin_steps OBJECT ${PARTRAC_DOLFIN_STEP_SOURCES})
  set_property(TARGET ${PROJECT_NAME}_dolfin_steps PROPERTY POSITION_INDEPENDENT_CODE ON)
  list(APPEND PARTRAC_STEP_TARGETS ${PROJECT_NAME}_dolfin_steps)
endif()

# Particles, meshes, initializers, steps and the run loop
add_library(${PROJECT_NAME}_particles_units OBJECT
  particles/ParticleSet.cpp
  particles/mesh.cpp
  particles/Topology.cpp
  particles/Initializer.cpp
  particles/stepping.cpp
  )
# Step objects in the loop's order
set_source_files_properties(${PARTRAC_STEP_OBJECTS} PROPERTIES EXTERNAL_OBJECT ON GENERATED ON)
add_library(${PROJECT_NAME}_particles $<TARGET_OBJECTS:${PROJECT_NAME}_particles_units>
  ${PARTRAC_STEP_OBJECTS})
add_dependencies(${PROJECT_NAME}_particles ${PARTRAC_STEP_TARGETS})
foreach(t ${PROJECT_NAME}_particles_units ${PARTRAC_STEP_TARGETS} ${PROJECT_NAME}_particles)
  target_include_directories(${t} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/particles)
  target_link_libraries(${t} PUBLIC ${PROJECT_NAME}_interpol)
endforeach()

# What every app links against
add_library(${PROJECT_NAME}_core INTERFACE)
target_link_libraries(${PROJECT_NAME}_core INTERFACE ${PROJECT_NAME}_particles)
