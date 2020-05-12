# partrac
**partrac** is a **par**ticle **trac**ker that uses trilinear interpolation of cubically ordered velocity data to advect passive (possibly diffusive) particles. The typical input for **partrac** is the output of the Lattice Boltzmann code FELBM, but any other data can be used as input by a suitable conversion. It is written in C++.

## Conversion from FELBM output
`parse_xdmf.py` translates from FELBM output XDMF file to **partrac** input. This is run by:
`python3 parse_xdmf.py data_example/L64x256_a1/felbm_output/output.xdmf`
This creates a file `timestamps.dat` in the same folder as `output.xdmf` that is used as input for **partrac**.

## Compilation
Straightforward:
`make clean && make`

## Running
Passive tracers example:
`./trace data_example/L64x256_a1/felbm_output/timestamps.dat Dm=0 T=591000 dt=1.0 Nrw=10000 dump_intv=100.0 stat_intv=100.0 checkpoint_intv=1000 verbose=true x0=0 y0=0 z0=0 int_order=1 init_mode=line_x init_weight=uy write_mode=hdf5 interpolation_test=0 dump_chunk_size=10 refine=false refine_intv=10 hist_chunk_size=0 ds_max=0.05 Nrw_max=10000`
This creates the folder `data_example/L64x256_a1/felbm_output/RandomWalkers/Dm0..../` and puts the simulation data into it.

## Visualization
Plotting the position:
`python3 plot_pos.py data_example/L64x256_a1/felbm_output/RandomWalkers/Dm0..../`

## Parameters
| Parameter          |  Default value | Description                                            |
|--------------------|----------------|--------------------------------------------------------|
| folder             | ""             | Folder to store files in                               |
| restart_folder     | ""             | Folder to restart from                                 |
| Dm                 | 1.0            | Diffusion constant                                     |
| t0                 | 0.0            | Initial simulation time                                |
| T                  | 10000000.0     | Total simulation time                                  |
| Nrw                | 100            | Initial number of particles                            |
| dump_intv          | 100.0          | Time interval between dumping particle data            |
| stat_intv          | 100.0          | Time interval between dumpting ensemble statistics     |
| checkpoint_intv    | 1000.0         | Time interval between each checkpoint                  |
| dump_chunk_size    | 50             | Number of time stamps in a single hdf5 file            |
| verbose            | false          | Verbose output                                         |
| U0                 | 1.0            | Rescale velocity (not active)                          |
| x0                 | 0.0            | x parameter for initial distribution                   |
| y0                 | 0.0            | y parameter for initial distribution                   |
| z0                 | 0.0            | z parameter for initial distribution                   |
| dt                 | 1.0            | Time step size                                         |
| int_order          | 1              | Explicit integration order (1 or 2)                    |
| init_mode          | "line_x"       | Initialization mode                                    |
| init_weight        | "none"         | Initialization weighting                               |
| write_mode         | "hdf5"         | Dump mode ("hdf5" or "txt")                            |
| interpolation_test | 0              | Number of particles to test the interpolation with     |
| refine             | false          | Refine edges                                           |
| refine_intv        | 100.0          | Time interval between when to refine                   |
| hist_chunk_size    | 10             | Number of stats_intv between when to output histograms |
| ds_max             | 1.0            | Maximum accepted edge length                           |
| Nrw_max            | -1             | Maximum number of particles                            |

Parameters that are also stored in the event of a restart:
`t, Lx, Ly, Lz, nx, ny, nz, n_accepted, n_declined`
