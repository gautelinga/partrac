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
  Tet.cpp Tet.hpp
  Triangle.cpp Triangle.hpp
  dolfin_helpers.cpp dolfin_helpers.hpp
)


find_package(OpenMP REQUIRED)
#target_link_libraries(${APP_NAME})

target_link_libraries(${PROJECT_NAME}_interpol stdc++fs ${HDF5_LIBRARIES}
  ${HDF5_HL_LIBRARIES} Eigen3::Eigen OpenMP::OpenMP_CXX)
if (PARTRAC_ENABLE_DOLFIN)
  target_link_libraries(${PROJECT_NAME}_interpol dolfin)
endif()
target_compile_features(${PROJECT_NAME}_interpol PRIVATE cxx_std_17)