# partrac
Particle tracker

`python3 parse_xdmf.py data_example/L64x256_a1/felbm_output/output.xdmf`

`timestamps.dat`

`make clean && make`

`./trace data_example/L64x256_a1/felbm_output/timestamps.dat Dm=0 T=591000 dt=1.0 Nrw=10000 dump_intv=100.0 stat_intv=100.0 checkpoint_intv=1000 verbose=true x0=0 y0=0 z0=0 int_order=1 init_mode=line_x init_weight=uy write_mode=hdf5 interpolation_test=0 dump_chunk_size=10 refine=false refine_intv=10 hist_chunk_size=0 ds_max=0.05 Nrw_max=10000`

`data_example/L64x256_a1/felbm_output/RandomWalkers/Dm0..../tdata.dat`

`python3 plot_pos.py data_example/L64x256_a1/felbm_output/RandomWalkers/Dm0..../`
