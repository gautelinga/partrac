option(BUILD_APP_INTERPOL "Build interpol app" TRUE)
option(BUILD_APP_FILAMENTS "Build filaments app" TRUE)
option(BUILD_APP_FILAMENTS_TRIANGLERK4 "Build filaments_triangleRK4 app" TRUE)
option(BUILD_APP_FILAMENTS_FELBMRK4 "Build filaments_felbmRK4 app" TRUE)
option(BUILD_APP_TRACERS_TRIANGLERK4 "Build tracers_triangleRK app" TRUE)
option(BUILD_APP_TRACERVECTORS_TRIANGLERK4 "Build tracervectors_triangleRK4 app" TRUE)
option(BUILD_APP_TRACERTENSORS_TRIANGLERK4 "Build tracertensors_triangleRK4 app" TRUE)
option(BUILD_APP_TRACERVECTORS_TRIANGLEFREQRK4 "Build tracervectors_trianglefreqRK4 app" TRUE)
option(BUILD_APP_WEIGHTED_WALKERS "Build weighted_walkers app" TRUE)
option(BUILD_APP_SPATIAL "Build spatial app" TRUE)
option(BUILD_APP_TESTS "Build test apps" TRUE)

partrac_add_app(partrac SOURCE partrac/partrac.cpp)
if (BUILD_APP_INTERPOL)
    partrac_add_app(interpol SOURCE interpol/interpol.cpp)
endif()
if (BUILD_APP_FILAMENTS)
    partrac_add_app(filaments SOURCE filaments/filaments.cpp)
endif()
if (BUILD_APP_FILAMENTS_TRIANGLERK4)
    partrac_add_app(filaments_triangleRK4
      SOURCE filaments_triangleRK4/filaments_triangleRK4.cpp NEEDS_DOLFIN)
endif()
if (BUILD_APP_FILAMENTS_FELBMRK4)
    partrac_add_app(filaments_felbmRK4
      SOURCE filaments_felbmRK4/filaments_felbmRK4.cpp)
endif()
if (BUILD_APP_TRACERS_TRIANGLERK4)
    partrac_add_app(tracers_triangleRK4
      SOURCE tracers_triangleRK4/tracers_triangleRK4.cpp NEEDS_DOLFIN)
endif()
if (BUILD_APP_TRACERVECTORS_TRIANGLERK4)
    partrac_add_app(tracervectors_triangleRK4
      SOURCE tracervectors_triangleRK4/tracervectors_triangleRK4.cpp NEEDS_DOLFIN)
endif()
if (BUILD_APP_TRACERVECTORS_TRIANGLEFREQRK4)
    partrac_add_app(tracervectors_trianglefreqRK4
      SOURCE tracervectors_trianglefreqRK4/tracervectors_trianglefreqRK4.cpp NEEDS_DOLFIN)
    partrac_add_app(tracervectors_analyticRK4
      SOURCE tracervectors_analyticRK4/tracervectors_analyticRK4.cpp)
endif()
if (BUILD_APP_TRACERTENSORS_TRIANGLERK4)
    partrac_add_app(tracertensors_triangleRK4
      SOURCE tracertensors_triangleRK4/tracertensors_triangleRK4.cpp NEEDS_DOLFIN)
endif()
if (BUILD_APP_SPATIAL)
    partrac_add_app(static_space_stepper
      SOURCE spatial/static_space_stepper.cpp)
    partrac_add_app(tracervectors_triangle_spatial
      SOURCE tracervectors_triangle_spatial/tracervectors_triangle_spatial.cpp NEEDS_DOLFIN)
endif()
if (BUILD_APP_WEIGHTED_WALKERS)
    partrac_add_app(weighted_walkers
      SOURCE weighted_walkers/weighted_walkers.cpp NEEDS_DOLFIN)
endif()
if (BUILD_APP_TESTS)
    partrac_add_app(omp_test SOURCE omp_test/omp_test.cpp NEEDS_DOLFIN)
endif()
