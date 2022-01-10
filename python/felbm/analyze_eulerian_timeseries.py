import argparse
import sys
import os
basedir = os.path.join(os.path.dirname(__file__), "..")
sys.path.append(basedir)
from utils import Params, read_timestamps
import numpy as np
import os
import h5py
from scipy import signal
import mpi4py.MPI as MPI
from tqdm import tqdm

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def mpi_print(*args):
    if rank == 0:
        print(*args)


def parse_args():
    parser = argparse.ArgumentParser(description="Plot pos")
    parser.add_argument("folder", type=str, help="Folder")
    parser.add_argument("-t0", type=float, default=None, help="")
    parser.add_argument("-t1", type=float, default=None, help="")
    parser.add_argument("-axis", type=int, default=0, help="Axis")
    parser.add_argument("--show", action="store_true", help="Show plot")
    parser.add_argument("--phase_by_density", action="store_true", help="Show plot")
    parser.add_argument("--blocks", type=int, default=1, help="Number of blocks to process")
    args = parser.parse_args()
    return args


def select_timestamps(felbm_folder, t0, t1):
    timestamps = read_timestamps(os.path.join(felbm_folder, "timestamps.dat"))
    timestamps_out = []
    for timestamp in timestamps:
        if (t0 is None or timestamp[0] >= t0) and (t1 is None or timestamp[0] <= t1):
            timestamps_out.append(timestamp)
    #timestamps = timestamps_out
    mpi_print(timestamps_out)
    return timestamps_out


if __name__ == "__main__":
    args = parse_args()
    #params = Params(args.folder)
    #t0 = params.get_tmin()
    #Lx = params.get("Lx", t0)

    felbm_folder = args.folder
    analysisfolder = os.path.join(felbm_folder, "Analysis")
    if rank == 0 and not os.path.exists(analysisfolder):
        os.makedirs(analysisfolder)

    timestamps = select_timestamps(felbm_folder, args.t0, args.t1)

    with h5py.File(os.path.join(felbm_folder, "output_is_solid.h5"), "r") as h5f:
        is_solid = np.array(h5f["is_solid"])
        nx, ny, nz = is_solid.shape
        mpi_print(nx, ny, nz)
        is_solid = is_solid.reshape((nz, ny, nx))

    #is_solid = is_solid[nz // 2, :, :].flatten()
    #is_fluid = np.logical_not(is_solid)
    is_fluid = np.logical_not(is_solid)
    is_plane = np.zeros_like(is_solid, dtype=bool)
    is_plane[nz // 2, :, :] = True
    is_domain = np.logical_and(is_plane, is_fluid)
    is_domain_flat = is_domain.flatten()
    is_fluid_xy = is_fluid[nz//2, :, :]
    fluid2xy = np.argwhere(is_fluid_xy.flatten()).flatten()
    mpi_print(len(fluid2xy), fluid2xy[-1])

    #print(np.sum(is_domain_flat), np.sum(is_fluid_xy_flat))
    #p = np.zeros_like(is_solid, dtype=float)
    rho = np.zeros_like(is_fluid_xy, dtype=float).flatten()
    rho_loc = np.zeros_like(rho)
    #c = np.zeros_like(is_solid, dtype=float)
    #ux = np.zeros_like(is_solid, dtype=float)
    #uy = np.zeros_like(is_solid, dtype=float)
    ux = np.zeros_like(is_fluid_xy, dtype=float).flatten()
    ux_loc = np.zeros_like(ux)
    du = np.zeros_like(rho)
    du_loc = np.zeros_like(du)
    # x, y = np.meshgrid(np.arange(nx) + .5, np.arange(ny) + 0.5)
    freq = np.zeros_like(ux)
    freq_loc = np.zeros_like(freq)

    Nt = len(timestamps)

    #print(is_domain.shape)

    Ndofs = int(np.sum(is_domain.flatten()))
    mpi_print("porosity={}".format(Ndofs/(nx*ny)))
    mpi_print("Is fluid: {}".format(Ndofs))
    block_size = int(np.ceil(Ndofs / args.blocks))

    blocks = list(range(args.blocks))
    bstart = args.blocks * rank // size
    bstop = args.blocks * (rank + 1) // size
    blocks_proc = blocks[bstart:bstop]

    t_ = np.array([t for t, _ in timestamps])

    comm.Barrier()
    if rank == 0:
        pbar = tqdm(total=len(blocks_proc)*len(timestamps))

    for iblock in blocks_proc:
        istart = iblock * block_size
        istop = min((iblock + 1) * block_size, Ndofs)
        uxt_block = np.zeros((istop-istart, len(timestamps)))
        uyt_block = np.zeros_like(uxt_block)
        rho_block = np.zeros_like(uxt_block)
        #Ux = np.zeros((Nt, ny, nx))
        #Uy = np.zeros_like(Ux)
        for it, timestamp in enumerate(timestamps):
            t = timestamp[0]
            h5filename = timestamp[1]
            with h5py.File(os.path.join(felbm_folder, h5filename), "r") as h5f:
                #p[:, :] = np.array(h5f["pressure"]).reshape((nz, ny, nx))[nz // 2, :, :]
                #rho[:, :] = np.array(h5f["density"]).reshape((nz, ny, nx))[nz // 2, :, :]
                #ux[:, :] = np.array(h5f["u_x"]).reshape((nz, ny, nx))[nz // 2, :, :]
                #uy[:, :] = np.array(h5f["u_y"]).reshape((nz, ny, nx))[nz // 2, :, :]
                rho_block[:, it] = np.array(h5f["density"]).flatten()[is_domain_flat][istart:istop]
                uxt_block[:, it] = np.array(h5f["u_x"]).flatten()[is_domain_flat][istart:istop]
                uyt_block[:, it] = np.array(h5f["u_y"]).flatten()[is_domain_flat][istart:istop]
                #print(w)
                if rank == 0:
                    pbar.update(1)
        freq_block = np.zeros(uxt_block.shape[0])
        rhoavg_block = rho_block.mean(axis=1)
        for i in range(uxt_block.shape[0]):
            f, P2 = signal.welch(uxt_block[i, :], 1./(t_[1]-t_[0]), 'flattop', uxt_block.shape[1], scaling='spectrum')
            P = np.sqrt(P2)
            sumP = np.sum(P)
            freq_block[i] = np.sum(f*P)/sumP if sumP > 0. else 0.
            #plt.plot(t_, uxt_block[0, :])
            #plt.plot(f, Pxx_spec)
            #plt.show()
        ids = fluid2xy[istart:istop]
        ux_loc[ids] = np.sqrt(np.mean(uxt_block**2, axis=1))
        freq_loc[ids] = freq_block[:]
        rho_loc[ids] = rhoavg_block[:]
        du_loc[ids] = np.sqrt(np.mean((uxt_block - np.outer(uxt_block.mean(axis=1), np.ones(uxt_block.shape[1])))**2 
                                      + (uyt_block - np.outer(uyt_block.mean(axis=1), np.ones(uyt_block.shape[1])))**2, axis=1))

    comm.Reduce(ux_loc, ux, op=MPI.SUM, root=0)
    comm.Reduce(freq_loc, freq, op=MPI.SUM, root=0)
    comm.Reduce(rho_loc, rho, op=MPI.SUM, root=0)
    comm.Reduce(du_loc, du, op=MPI.SUM, root=0)

    comm.Barrier()
    if rank == 0:
        freq = freq.reshape((ny, nx))
        ux = ux.reshape((ny, nx))
        rho = rho.reshape((ny, nx))
        du = du.reshape((ny, nx))

        np.savetxt(os.path.join(analysisfolder, "uxnorm_avg.dat"), ux)
        np.savetxt(os.path.join(analysisfolder, "freq_avg.dat"), freq)
        np.savetxt(os.path.join(analysisfolder, "rho_avg.dat"), rho)
        np.savetxt(os.path.join(analysisfolder, "is_fluid_xy.dat"), is_fluid_xy)
        np.savetxt(os.path.join(analysisfolder, "du.dat"), du)
