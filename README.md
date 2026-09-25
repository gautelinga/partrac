# partrac
**partrac** is a **par**ticle **trac**ker that can advect passive (possibly diffusive) particles in time-dependent velocity fields. On top of this, stretching of lines and sheets with automatic refinement/coarsening is possible, which implements the diffusive strip/sheet method (post-processing in the sibling repository `partractools`). It is written in C++.

## Input modes
* `structured/lbm/felbm`: Uses trilinear interpolation of cubically ordered velocity data. The typical input for **partrac** is the output of the Lattice Boltzmann code FELBM, but any other data can be used as input by a suitable conversion.
* `fenics`: Uses unstructured meshes, in the form of Fenics/Dolfin HDF5 files. Needs dolfin.
* `triangle/tet`: Better and faster implementation of `fenics`, without dolfin.
* `trianglefreq/tetfreq`: As `triangle/tet`, a time series given as frequency components.
* `xdmftriangle/xdmftet`: P1 fields in XDMF files.
* `openfoam`: OpenFOAM cases, read in place, with `partrac_params.dat` beside `constant/`.
* `analytic`: Analytic input

## Conversion from FELBM output
`python/parse_xdmf.py` translates from FELBM output XDMF file to **partrac** input. This is run by:
`python3 python/parse_xdmf.py FELBM_OUTPUT/output.xdmf`
This creates a file `timestamps.dat` in the same folder as `output.xdmf` that is used as input for **partrac**.

## Compilation
```
cmake -S . -B build
make -C build -j
```
The executables end up in `build/bin/`. `-DPARTRAC_ENABLE_OPENFOAM=ON` adds
`mode=openfoam`, built against an installed OpenFOAM. To run the tests:
```
ctest --test-dir build --output-on-failure
```

## Running
Passive tracers example:
`./build/bin/partrac data_example/plane_poiseuille/expr_params.dat mode=analytic init_mode=uniform_x Nrw=100 Nrw_max=10000 ds_max=0.4 ds_min=0.1 Dm=0 dt=0.01 T=1.0 int_order=1 dump_intv=0.1 stat_intv=0.1`
This creates the folder `data_example/plane_poiseuille/RandomWalkers/Dm0..../` and puts the simulation data into it.

## Apps
One app per kind of thing followed, all on the same run loop. Every app reads
its field through any interpolator (`mode=`); the ones that step in time take
`scheme=explicit` (with noise when `Dm > 0`) or `scheme=RK4`.

| App | What it follows |
|-----|-----------------|
| `partrac` | Points, strips and sheets, with refinement, coarsening and injection |
| `filaments` | Pairs of points and their stretching; `resize=` and `outside=reinject` for edges stuck in an underresolved field |
| `tracers` | A cloud of points |
| `tracervectors` | Points carrying a material line element |
| `tracertensors` | Points carrying the deformation gradient |
| `static_space_stepper` | Points, strips and sheets marched in path length along the streamlines, the fields frozen |
| `tracervectors_spatial` | Line elements marched in path length |
| `weighted_walkers` | Diffusive walkers past an exit plane, resampled by weight |
| `interpol` | Probes the fields at random points |

The older per-interpolator names (`tracervectors_triangleRK4`,
`filaments_felbmRK4`, ...) still run: each is its app with the interpolator
and the choices it used to fix pinned (`apps/CMakeLists.txt`).

## Mesh examples
The `data_example` folders for the mesh modes (`ppf_triangle_p2`,
`test_triangle_p2`, `test_tet_p1`, `test_tet_p2`, `sine_trianglefreq_p2`,
`sine_tetfreq_p2`) ship a
`generate_up.py` rather than the mesh itself. Run it inside the folder to write
`mesh.h5` and `up_0.h5`:
```
cd data_example/ppf_triangle_p2 && python3 generate_up.py -dim 1
```
It needs FEniCS/dolfin.

## Divergence-free velocity fields
`python/divfree/divfree_clean.py` prepares a dolfin HDF5 case so that the
velocity is divergence-free in every cell, which keeps tracers from stopping at
no-slip walls. The output is a case of its own whose parameter file carries
`divfree=true`, the key the loaders read it back with.
```
python3 python/divfree/divfree_clean.py CASE/dolfin_params.dat --out CLEANED
mpirun -n 8 python3 python/divfree/divfree_clean.py CASE/dolfin_params.dat --out CLEANED
python3 python/divfree/divfree_clean.py CLEANED/dolfin_params.dat --check
```
It needs `h5py`, `scipy`, `petsc4py` and `mpi4py`, not dolfin; `--help` lists the
options. Under `mpirun` each rank holds its part of the mesh; run it with
`OMP_NUM_THREADS=1`.

## Visualization
Plotting the position:
`python3 python/plot_pos.py FOLDER/RandomWalkers/Dm0..../`

## Parameters
Each app declares the parameters it accepts, so the set differs between them and
an unrecognised parameter is an error rather than being silently ignored. Run an
app with `--help` for its own list, with types, defaults and which parameters are
required:
```
./build/bin/partrac --help
```
Parameters are given as `key=value` after the input file. A parameter is either
required, optional with a default, or computed by the program (`folder`, `t`,
`Lx`, `Ly`, `Lz`) and then only read back when restarting. Some are required only
in certain configurations, for instance `La` when `init_mode` is a strip, sheet or
ellipsoid.

### Initialization modes
Set with `init_mode`. The trailing axes select the direction(s) involved.

| Mode                         | Description                                             | Also reads              |
|------------------------------|---------------------------------------------------------|-------------------------|
| `point`                      | Nrw particles at a single point                          | x0, y0, z0              |
| `uniform_[xyz]`              | Nrw particles spread uniformly along an axis             |                         |
| `strip_[xyz]`                | Nrw particles along a strip of length La                 | La                      |
| `sheet_[xy,xz,yz]`           | A refined sheet spanning La by Lb                        | La, Lb, ds_init         |
| `ellipsoid_[xy,xz,yz]`       | A refined ellipsoid with semi-axes La and Lb             | La, Lb, x0, y0, z0      |
| `pair_[xyz]`                 | Two particles separated by ds_init, randomly oriented    | ds_init, x0, y0, z0     |
| `pairs_[xyz]_[xyz]`          | Nrw/2 pairs separated by ds_init, randomly oriented      | ds_init                 |
| `points_[xyz]`               | Nrw particles at random positions                        | init_weight, ds_init    |
| `randomgaussianstrip_[xyz]_[xyz]` | A strip of length La with gaussian spread Lb        | La, Lb                  |
| `randomgaussiancircle_[xyz]` | A circle of diameter La with gaussian spread Lb          | La, Lb                  |
| `from_file:<file.h5>`        | Positions read from a file                               | x0, y0, z0, t0          |

`init_weight` selects how `points_*` samples positions: `none` (uniform), `u`, or
one velocity component `ux`, `uy`, `uz`.
