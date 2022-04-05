include_directories(
  ${PROJECT_SOURCE_DIR}/src
  ${HDF5_INCLUDE_DIRS}
)

add_executable(${PROJECT_NAME} partrac.cpp)

if (PARTRAC_ENABLE_DOLFIN)
  target_link_libraries(${PROJECT_NAME} dolfin)
endif()

if (PARTRAC_ENABLE_MPI)
  target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})
endif()

if(MPI_COMPILE_FLAGS)
  set_target_properties(${PROJECT_NAME}  PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${PROJECT_NAME}  PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

target_link_libraries(${PROJECT_NAME} ${HDF5_LIBRARIES}
  ${HDF5_HL_LIBRARIES} Eigen3::Eigen ${PROJECT_NAME}_io ${PROJECT_NAME}_utils ${PROJECT_NAME}_interpol)
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)

set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/")