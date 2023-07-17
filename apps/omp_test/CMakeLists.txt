if (PARTRAC_ENABLE_DOLFIN)

include_directories(
  ${PROJECT_SOURCE_DIR}/src
)

set(APP_NAME omp_test)
set(PROJECT_TARGETS "${PROJECT_TARGETS} ${APP_NAME}" PARENT_SCOPE)
#message(STATUS ${PROJECT_TARGETS})

find_package(OpenMP REQUIRED)

add_executable(${APP_NAME} omp_test.cpp)

target_link_libraries(${APP_NAME} dolfin)
target_link_libraries(${APP_NAME} OpenMP::OpenMP_CXX)

if (PARTRAC_ENABLE_MPI)
  target_link_libraries(${APP_NAME} ${MPI_LIBRARIES})
endif()

if(MPI_COMPILE_FLAGS)
  set_target_properties(${APP_NAME}  PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${APP_NAME}  PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

target_link_libraries(${APP_NAME} ${HDF5_LIBRARIES}
  ${HDF5_HL_LIBRARIES} Eigen3::Eigen ${PROJECT_NAME}_io ${PROJECT_NAME}_utils ${PROJECT_NAME}_interpol)
target_compile_features(${APP_NAME} PRIVATE cxx_std_17)

set_target_properties(${APP_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/")

endif()