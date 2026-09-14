option(BUILD_APP_INTERPOL "Build interpol app" TRUE)
option(BUILD_APP_FILAMENTS "Build filaments app" TRUE)
option(BUILD_APP_WEIGHTED_WALKERS "Build weighted_walkers app" TRUE)
option(BUILD_APP_SPATIAL "Build spatial app" TRUE)

partrac_add_app(partrac SOURCE partrac/partrac.cpp)
if (BUILD_APP_INTERPOL)
    partrac_add_app(interpol SOURCE interpol/interpol.cpp)
endif()
if (BUILD_APP_FILAMENTS)
    partrac_add_app(filaments SOURCE filaments/filaments.cpp)
endif()
# Retired filaments names
partrac_add_wrapper(filaments_triangleRK4 TARGET filaments
  ARGS mode=triangle scheme=RK4 resize=doublings resize_target=ds_init outside=reinject)
partrac_add_wrapper(filaments_felbmRK4 TARGET filaments
  ARGS mode=felbm scheme=RK4 resize=doublings resize_target=ds_init outside=reinject)
# Tracer apps and retired names
partrac_add_app(tracers SOURCE tracers/tracers.cpp)
partrac_add_wrapper(tracers_triangleRK4 TARGET tracers
  ARGS mode=triangle scheme=RK4 outside=reinject)
partrac_add_app(tracervectors SOURCE tracervectors/tracervectors.cpp)
partrac_add_wrapper(tracervectors_analyticRK4 TARGET tracervectors
  ARGS mode=analytic scheme=explicit outside=ignore output_J=true output_phi=false output_cell_type=false
  DROP print_found)
partrac_add_wrapper(tracervectors_trianglefreqRK4 TARGET tracervectors
  ARGS mode=trianglefreq scheme=explicit outside=reinject output_phi=false output_cell_type=false
  DROP print_found)
partrac_add_wrapper(tracervectors_triangleRK4 TARGET tracervectors
  ARGS mode=xdmftriangle
  DROP print_found)
partrac_add_app(tracertensors SOURCE tracertensors/tracertensors.cpp)
partrac_add_wrapper(tracertensors_triangleRK4 TARGET tracertensors
  ARGS mode=tet scheme=RK4 outside=reinject)
if (BUILD_APP_SPATIAL)
    # Marching apps and retired name
    partrac_add_app(static_space_stepper
      SOURCE spatial/static_space_stepper.cpp)
    partrac_add_app(tracervectors_spatial
      SOURCE tracervectors_spatial/tracervectors_spatial.cpp)
    partrac_add_wrapper(tracervectors_triangle_spatial TARGET tracervectors_spatial
      ARGS mode=triangle
      RENAME Lt=Ln COPY dt=dxn DROP T ds_max)
endif()
if (BUILD_APP_WEIGHTED_WALKERS)
    partrac_add_app(weighted_walkers
      SOURCE weighted_walkers/weighted_walkers.cpp)
endif()
