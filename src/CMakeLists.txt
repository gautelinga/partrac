# Parameter schema layer. No dependencies, so the unit tests link it on its own.
add_library(${PROJECT_NAME}_params
  Params.cpp Params.hpp
  )
target_include_directories(${PROJECT_NAME}_params PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(${PROJECT_NAME}_params PUBLIC cxx_std_17)
target_compile_options(${PROJECT_NAME}_params PRIVATE -Wall -Wextra)

add_library(${PROJECT_NAME}_utils
  Timestamps.cpp Timestamps.hpp
  FreqStamps.hpp
  )
target_link_libraries(${PROJECT_NAME}_utils Eigen3::Eigen)
target_compile_features(${PROJECT_NAME}_utils PRIVATE cxx_std_17)

add_library(${PROJECT_NAME}_io
  io.hpp io.cpp
  )
target_link_libraries(${PROJECT_NAME}_io stdc++fs ${HDF5_LIBRARIES}
  ${HDF5_HL_LIBRARIES})
target_compile_features(${PROJECT_NAME}_io PRIVATE cxx_std_17)

# Only the analytic and structured interpolators are header-only; every source
# file here needs dolfin, so the library exists only when dolfin is enabled.
if (PARTRAC_ENABLE_DOLFIN)
add_library(${PROJECT_NAME}_interpol
  Interpol.hpp
  AnalyticInterpol.hpp
  StructuredInterpol.hpp
  StructuredConstInterpol.hpp
  DolfInterpol.cpp DolfInterpol.hpp
  TetInterpol.cpp TetInterpol.hpp
  TriangleInterpol.cpp TriangleInterpol.hpp
  TriangleFreqInterpol.cpp TriangleFreqInterpol.hpp
  XDMFTriangleInterpol.cpp XDMFTriangleInterpol.hpp
  XDMFTetInterpol.cpp XDMFTetInterpol.hpp
  Tet.cpp Tet.hpp
  Triangle.cpp Triangle.hpp
  dolfin_helpers.cpp dolfin_helpers.hpp
  xdmf_helpers.cpp xdmf_helpers.hpp
)


target_link_libraries(${PROJECT_NAME}_interpol PUBLIC stdc++fs ${HDF5_LIBRARIES}
  ${HDF5_HL_LIBRARIES} Eigen3::Eigen OpenMP::OpenMP_CXX dolfin)
target_compile_features(${PROJECT_NAME}_interpol PRIVATE cxx_std_17)
endif()

# What every app links against
add_library(${PROJECT_NAME}_core INTERFACE)
target_link_libraries(${PROJECT_NAME}_core INTERFACE
  ${PROJECT_NAME}_io ${PROJECT_NAME}_utils ${PROJECT_NAME}_params
  stdc++fs ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES}
  Eigen3::Eigen OpenMP::OpenMP_CXX)
target_include_directories(${PROJECT_NAME}_core INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR} ${HDF5_INCLUDE_DIRS})
target_compile_features(${PROJECT_NAME}_core INTERFACE cxx_std_17)
if (PARTRAC_ENABLE_DOLFIN)
  target_link_libraries(${PROJECT_NAME}_core INTERFACE ${PROJECT_NAME}_interpol)
endif()
if (PARTRAC_ENABLE_MPI)
  target_link_libraries(${PROJECT_NAME}_core INTERFACE ${MPI_LIBRARIES})
endif()