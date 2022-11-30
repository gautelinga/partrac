import sys
import os
import h5py
import numpy as np
basedir = os.path.join(os.path.dirname(__file__), "..")
sys.path.append(basedir)
from utils import Params, read_timestamps
import mpi4py.MPI as MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def mpi_print(*args):
    if rank == 0:
        print(*args)

def select_timestamps(felbm_folder, t0, t1):
    timestamps = read_timestamps(os.path.join(felbm_folder, "timestamps.dat"))
    timestamps_out = []
    for timestamp in timestamps:
        if (t0 is None or timestamp[0] >= t0) and (t1 is None or timestamp[0] <= t1):
            timestamps_out.append(timestamp)
    #timestamps = timestamps_out
    mpi_print(timestamps_out)
    return timestamps_out

def get_fluid_domain(felbm_folder):
    with h5py.File(os.path.join(felbm_folder, "output_is_solid.h5"), "r") as h5f:
        is_solid = np.array(h5f["is_solid"])
        nx, ny, nz = is_solid.shape
        mpi_print(nx, ny, nz)

    is_solid_xy = is_solid.reshape((nz, ny, nx))[nz // 2, :, :]
    is_fluid_xy = np.logical_not(is_solid_xy)
    return is_fluid_xy, (nx, ny, nz)

def make_folder_safe(analysisfolder):
    if rank == 0 and not os.path.exists(analysisfolder):
        os.makedirs(analysisfolder)