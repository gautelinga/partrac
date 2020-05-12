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
