#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

add_library(${PROJECT_NAME}_utils
  Timestamps.cpp
  )
target_link_libraries(${PROJECT_NAME}_utils Eigen3::Eigen)
target_compile_features(${PROJECT_NAME}_utils PRIVATE cxx_std_17)

add_library(${PROJECT_NAME}_io
  io.hpp io.cpp
  )
target_link_libraries(${PROJECT_NAME}_io stdc++fs ${HDF5_LIBRARIES}
  ${HDF5_HL_LIBRARIES})
target_compile_features(${PROJECT_NAME}_io PRIVATE cxx_std_17)

add_library(${PROJECT_NAME}_interpol
  Interpol.hpp
  AnalyticInterpol.hpp
  StructuredInterpol.hpp
  DolfInterpol.cpp DolfInterpol.hpp
  TetInterpol.cpp TetInterpol.hpp
  TriangleInterpol.cpp TriangleInterpol.hpp
  Tet.cpp Tet.hpp
  Triangle.cpp Triangle.hpp
)


target_link_libraries(${PROJECT_NAME}_interpol stdc++fs ${HDF5_LIBRARIES}
  ${HDF5_HL_LIBRARIES} Eigen3::Eigen)
if (PARTRAC_ENABLE_DOLFIN)
  target_link_libraries(${PROJECT_NAME}_interpol dolfin)
endif()
target_compile_features(${PROJECT_NAME}_interpol PRIVATE cxx_std_17)

#add_executable(${PROJECT_NAME} traj.cpp)

#if (PARTRAC_ENABLE_DOLFIN)
#  target_link_libraries(${PROJECT_NAME} dolfin)
#endif()
#target_link_libraries(${PROJECT_NAME} ${HDF5_LIBRARIES}
#  ${HDF5_HL_LIBRARIES} Eigen3::Eigen partrac_lib)
#target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)

#install(TARGETS ${PROJECT_NAME} DESTINATION /usr/bin/)
