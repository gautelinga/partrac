# partrac
**partrac** is a **par**ticle **trac**ker that can advect passive (possibly diffusive) particles in time-dependent velocity fields. On top of this, stretching of lines and sheets with automatic refinement/coarsening is possible, which implements the diffusive strip/sheet method (post-processing in the sibling repository `partractools`). It is written in C++.

## Input modes
* `structured/lbm/felbm`: Trilinear interpolation of cubically ordered velocity data, typically the output of the Lattice Boltzmann code [FELBM](https://bitbucket.org/ulbm/felbm/).
* `fenics`: Fenics/Dolfin HDF5 meshes and fields. Needs dolfin.
* `triangle/tet`: As `fenics`, faster, without dolfin.
* `trianglefreq/tetfreq`: As `triangle/tet`, a time series given as frequency components.
* `xdmftriangle/xdmftet`: P1 fields in XDMF files.
* `openfoam`: OpenFOAM cases, read in place, with `partrac_params.dat` in the case folder.
* `analytic`: Analytic input.

## Conversion from FELBM output
`python3 python/parse_xdmf.py FELBM_OUTPUT/output.xdmf` writes `timestamps.dat` next to `output.xdmf`, the input for **partrac**.

## Compilation
```
cmake -S . -B build
make -C build -j
```
The executables end up in `build/bin/`. Options:
* `-DPARTRAC_ENABLE_DOLFIN=OFF`: no `mode=fenics` (on by default).
* `-DPARTRAC_ENABLE_OPENFOAM=ON`: `mode=openfoam`, against an installed OpenFOAM.
* `-DCMAKE_BUILD_TYPE=Debug`: the default is `Release`.
* `-DPARTRAC_BUILD_<APP>=OFF`: leaves out `INTERPOL`, `FILAMENTS`, `WEIGHTED_WALKERS` or `SPATIAL`.
* `-DPARTRAC_SANITIZE=ON`, `-DPARTRAC_COVERAGE=ON`, `-DPARTRAC_BUILD_DOCS=ON`.

Tests:
```
ctest --test-dir build --output-on-failure
```
`python3 -m pytest tests --mpi` (or `PARTRAC_MPI=1`) adds the multi-rank divfree tests.

## Running
Passive tracers example:
`./build/bin/partrac data_example/plane_poiseuille/expr_params.dat mode=analytic init_mode=uniform_x Nrw=100 Nrw_max=10000 ds_max=0.4 ds_min=0.1 Dm=0 dt=0.01 T=1.0 int_order=1 dump_intv=0.1 stat_intv=0.1`
This creates the folder `data_example/plane_poiseuille/RandomWalkers/Dm0..../` and puts the simulation data into it.

`num_threads` sets the OpenMP threads, which wait passively unless `OMP_WAIT_POLICY` is set.

### Restarting
Give the folder that holds `Checkpoints/`:
```
./build/bin/partrac data_example/plane_poiseuille/expr_params.dat restart_folder=data_example/plane_poiseuille/RandomWalkers/Dm0..../0 T=2.0
```

## Apps
| App                     | What it follows                                                           |
|-------------------------|---------------------------------------------------------------------------|
| `partrac`               | Points, strips and sheets, with refinement, coarsening and injection      |
| `tracers`               | A cloud of points                                                         |
| `tracervectors`         | Points carrying a material line element                                   |
| `tracertensors`         | Points carrying the deformation gradient                                  |
| `static_space_stepper`  | Points, strips and sheets marched in path length along steady streamlines |
| `tracervectors_spatial` | Line elements marched in path length                                      |
| `weighted_walkers`      | Diffusive walkers past an exit plane, resampled by weight                 |
| `filaments`             | Pairs of points and their stretching                                      |
| `interpol`              | Probes the fields at random points                                        |

The older names (`tracervectors_triangleRK4`, `filaments_felbmRK4`, ...) still run, as their app with the old choices fixed.

### Time schemes
Set with `scheme`.
* `explicit`: Euler, or second order with `int_order=2`. Diffusive when `Dm > 0`.
* `RK4`: Fourth-order Runge-Kutta. Ignores `Dm`.
* `RK4cells`: RK4 with the steps cut at cell facets, `Dm=0`.

## Mesh examples
The mesh examples in `data_example` ship a `generate_up.py` (needs dolfin) rather than the mesh:
```
cd data_example/ppf_triangle_p2 && python3 generate_up.py -dim 1
```

## Visualization
Plotting the position:
`python3 python/plot_pos.py FOLDER/RandomWalkers/Dm0..../`

## Parameters
Parameters are given as `key=value` after the input file. `--help` lists an app's parameters, with types and defaults:
```
./build/bin/partrac --help
```
`wall_upgrade=true` raises the order of the velocity next to walls at rest (P1 to quadratic, P2 by quartic terms); it is on by default for `xdmftet`, `xdmftriangle` and `openfoam`.

### Initialization modes
Set with `init_mode`. The trailing axes select the direction(s) involved.

| Mode                              | Description                                              | Also reads              |
|-----------------------------------|----------------------------------------------------------|-------------------------|
| `point`                           | Nrw particles at a single point                          | x0, y0, z0              |
| `uniform_[xyz]`                   | Nrw particles spread uniformly along an axis             |                         |
| `strip_[xyz]`                     | Nrw particles along a strip of length La                 | La                      |
| `sheet_[xy,xz,yz]`                | A refined sheet spanning La by Lb                        | La, Lb, ds_init         |
| `ellipsoid_[xy,xz,yz]`            | A refined ellipsoid with semi-axes La and Lb             | La, Lb, x0, y0, z0      |
| `pair_[xyz]`                      | Two particles separated by ds_init, randomly oriented    | ds_init, x0, y0, z0     |
| `pairs_[xyz]_[xyz]`               | Nrw/2 pairs separated by ds_init, randomly oriented      | ds_init                 |
| `points_[xyz]`                    | Nrw particles at random positions, weighted by `init_weight` (`uniform`, `u`, `ux`, `uy`, `uz`) | init_weight; ds_init    |
| `randomgaussianstrip_[xyz]_[xyz]` | A strip of length La with gaussian spread Lb             | La, Lb                  |
| `randomgaussiancircle_[xyz]`      | A circle of diameter La with gaussian spread Lb          | La, Lb                  |
| `file:<file.h5>`                  | Positions read from a file                               | x0, y0, z0, t0          |

## Divergence-free velocity fields
`python/divfree/divfree_clean.py` makes a dolfin HDF5 case divergence-free in every cell, which keeps tracers from stopping at no-slip walls; the result is read with `divfree=true`:
```
python3 python/divfree/divfree_clean.py CASE/dolfin_params.dat --out CLEANED
mpirun -n 8 python3 python/divfree/divfree_clean.py CASE/dolfin_params.dat --out CLEANED
```
It needs `h5py`, `scipy`, `petsc4py` and `mpi4py`; run it under `mpirun` with `OMP_NUM_THREADS=1`.
