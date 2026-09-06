# Parameter schema layer.  Deliberately has no dependencies -- no HDF5, no
# DOLFIN, no MPI, no boost -- so the unit tests link it on its own.  New code,
# so it is held to stricter warnings than the legacy tree (which is built with
# -fpermissive inherited from DOLFIN_CXX_FLAGS).
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

#add_subdirectory(interpol)
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


find_package(OpenMP REQUIRED)
#target_link_libraries(${APP_NAME})

target_link_libraries(${PROJECT_NAME}_interpol stdc++fs ${HDF5_LIBRARIES}
  ${HDF5_HL_LIBRARIES} Eigen3::Eigen OpenMP::OpenMP_CXX)
if (PARTRAC_ENABLE_DOLFIN)
  target_link_libraries(${PROJECT_NAME}_interpol dolfin)
endif()
target_compile_features(${PROJECT_NAME}_interpol PRIVATE cxx_std_17)