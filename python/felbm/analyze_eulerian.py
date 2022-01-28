from analyze_eulerian_timeseries import *

from scipy.ndimage import measurements, filters
import matplotlib.pyplot as plt
from skimage import feature
from scipy import spatial

def parse_args():
    parser = argparse.ArgumentParser(description="Plot pos")
    parser.add_argument("folder", type=str, help="Folder")
    parser.add_argument("-t0", type=float, default=None, help="")
    parser.add_argument("-t1", type=float, default=None, help="")
    parser.add_argument("-axis", type=int, default=0, help="Axis")
    parser.add_argument("--show", action="store_true", help="Show plot")
    parser.add_argument("--phase_by_density", action="store_true", help="Show plot")
    parser.add_argument("-o", "--output_dir", type=str, default="", help="Output directory")
    args = parser.parse_args()
    return args


def get_cluster_sizes(phase):
    lw, num = measurements.label(phase)
    # periodic boundaries
    for y in range(lw.shape[0]):
        if lw[y, 0] > 0 and lw[y, -1] > 0:
            lw[lw == lw[y, -1]] = lw[y, 0]
    for x in range(lw.shape[1]):
        if lw[0, x] > 0 and lw[-1, x] > 0:
            lw[lw == lw[-1, x]] = lw[0, x]
    old2new = np.zeros(lw.max()+1, dtype=int)
    old = np.unique(lw)
    old2new[old] = np.arange(len(old))
    lw = old2new[lw]

    label_list = np.arange(lw.max())
    S = measurements.sum(phase, lw, label_list)
    return S, lw

def get_pores(is_solid_xy):
    if rank == 0:
        S, lw = get_cluster_sizes(is_solid_xy)
        rad_est = np.sqrt(S[1:].mean()/np.pi)
        print("estimated_radius = {}".format(rad_est))
        d = 2 * int(np.ceil(rad_est))

        n_buf = 5*d
        is_solid_xy_embed = np.zeros((ny + 2 * n_buf, nx + 2 * n_buf))
        is_solid_xy_embed[n_buf:-n_buf, n_buf:-n_buf] = is_solid_xy
        is_solid_xy_embed[:n_buf, :] = is_solid_xy_embed[ny:-n_buf, :]
        is_solid_xy_embed[-n_buf:, :] = is_solid_xy_embed[n_buf:2*n_buf, :]
        is_solid_xy_embed[:, :n_buf] = is_solid_xy_embed[:, nx:-n_buf]
        is_solid_xy_embed[:, -n_buf:] = is_solid_xy_embed[:, n_buf:2*n_buf]
        
        lw, num = measurements.label(is_solid_xy_embed)
        label_list = np.arange(lw.max())
        S = measurements.sum(is_solid_xy_embed, lw, label_list)
        x_cm = measurements.center_of_mass(is_solid_xy_embed, lw, label_list[S > 0])
        x_ = np.array(x_cm)[:, ::-1]
        #plt.imshow(is_solid_xy_embed)
        #plt.plot(x_[:, 1], x_[:, 0], '.')

        tri = spatial.Delaunay(x_)
        pts = tri.points
        edges = set()
        for ids in tri.simplices:
            id0, id1, id2 = sorted(ids)
            edges.add((id0, id1))
            edges.add((id0, id2))
            edges.add((id1, id2))
        edges = np.array(list(edges), dtype=int)
        x_mid = 0.5*(pts[edges[:, 0], :] + pts[edges[:, 1], :])
        dx = pts[edges[:, 1], :] - pts[edges[:, 0], :]
        lp = np.linalg.norm(dx, axis=1)-d
        n = np.zeros_like(dx)
        n[:, 0] = -dx[:, 1]
        n[:, 1] = dx[:, 0]
        nnorm = np.linalg.norm(n, axis=1)
        #print(nnorm)
        n[:, 0] /= nnorm
        n[:, 1] /= nnorm

        x_mid[:, 0] -= n_buf
        x_mid[:, 1] -= n_buf
        ids = np.logical_and(np.logical_and(x_mid[:, 0] < nx, x_mid[:, 1] < ny),
                            np.logical_and(x_mid[:, 0] >= 0, x_mid[:, 1] >= 0))
        x_mid = x_mid[ids, :]
        n = n[ids, :]
        lp = lp[ids]
    else:
        x_mid = None
        n = None
        lp = None
        d = None
    x_mid = comm.bcast(x_mid, root=0)
    n = comm.bcast(n, root=0)
    lp = comm.bcast(lp, root=0)
    d = comm.bcast(d, root=0)
    return x_mid, n, lp, d


if __name__ == "__main__":
    args = parse_args()

    felbm_folder = args.folder
    if args.output_dir != "":
        output_folder = args.output_dir
    else:
        output_folder = felbm_folder
    analysisfolder = os.path.join(output_folder, "Analysis")
    if rank == 0 and not os.path.exists(analysisfolder):
        os.makedirs(analysisfolder)

    timestamps = select_timestamps(felbm_folder, args.t0, args.t1)

    with h5py.File(os.path.join(felbm_folder, "output_is_solid.h5"), "r") as h5f:
        is_solid = np.array(h5f["is_solid"])
        nx, ny, nz = is_solid.shape
        mpi_print(nx, ny, nz)
    is_solid_xy = is_solid.reshape((nz, ny, nx))[nz // 2, :, :]
    is_fluid_xy = np.logical_not(is_solid_xy)
    
    x_p, n_p, l_p, d = get_pores(is_solid_xy)
    ij_p = x_p[:, ::-1].round().astype(int)
    ij_p[:, 0] %= ny
    ij_p[:, 1] %= nx

    rho = np.zeros_like(is_fluid_xy, dtype=float)
    ux = np.zeros_like(rho)
    uy = np.zeros_like(rho)
    p = np.zeros_like(rho)

    t_ = np.array([t for t, _ in timestamps])
    uxt = np.zeros_like(t_)
    uyt = np.zeros_like(t_)
    uxt_loc = np.zeros_like(uxt)
    uyt_loc = np.zeros_like(uyt)
    S1t_loc = np.zeros_like(t_)
    S2t_loc = np.zeros_like(t_)
    S1t = np.zeros_like(t_)
    S2t = np.zeros_like(t_)

    istart = (rank * len(timestamps)) // size
    istop = ((rank + 1) * len(timestamps)) // size
    timestamps_block = list(enumerate(timestamps))[istart:istop]

    comm.Barrier()
    if rank == 0:
        pbar = tqdm(total=len(timestamps_block))

    un_p = np.zeros((len(ij_p), len(timestamps)))
    un_p_loc = np.zeros_like(un_p)

    for it, timestamp in timestamps_block:
        t = timestamp[0]
        h5filename = timestamp[1]
        with h5py.File(os.path.join(felbm_folder, h5filename), "r") as h5f:
            #p[:, :] = np.array(h5f["pressure"]).reshape((nz, ny, nx))[nz // 2, :, :]
            rho[:, :] = np.array(h5f["density"]).reshape((nz, ny, nx))[nz // 2, :, :]
            ux[:, :] = np.array(h5f["u_x"]).reshape((nz, ny, nx))[nz // 2, :, :]
            uy[:, :] = np.array(h5f["u_y"]).reshape((nz, ny, nx))[nz // 2, :, :]

        # TODO: only works when there is a mass density difference between phases!!
        rho_mean = rho[is_fluid_xy].mean()
        phase1 = np.logical_and(rho < rho_mean, is_fluid_xy)
        phase2 = np.logical_and(np.logical_not(phase1), is_fluid_xy)
        #rho1_mean = rho[phase1].mean()
        #rho2_mean = rho[phase2].mean()
        #eps = (rho2_mean - rho1_mean) / 10
        #interface = np.logical_and(np.logical_and(rho > rho1_mean + eps, rho < rho2_mean - eps), is_fluid_xy)
        #interface = np.logical_and(filters.gaussian_filter(interface.astype(float), sigma=.3, mode='wrap').astype(bool), is_fluid_xy)

        #dp = np.ma.masked_where(np.logical_not(interface), np.zeros_like(p))
        #dp[interface] = p[interface]
        #pp = np.ma.masked_where(is_solid_xy, p)
        #plt.imshow(dp)
        #plt.show()

        S1, lw1 = get_cluster_sizes(phase1)
        S2, lw2 = get_cluster_sizes(phase2)

        S1t_loc[it] = S1.mean()
        S2t_loc[it] = S2.mean()
        uxt_loc[it] = ux[is_fluid_xy].mean()
        uyt_loc[it] = uy[is_fluid_xy].mean()

        if rank == 0:
            pbar.update(1)

        #for k in range(len(ij_p)):
        i, j = ij_p[:, 0], ij_p[:, 1]
        un_p_loc[:, it] = n_p[:, 0] * ux[i, j] + n_p[:, 1] * uy[i, j]
        
    comm.Reduce(uxt_loc, uxt, op=MPI.SUM, root=0)
    comm.Reduce(uyt_loc, uyt, op=MPI.SUM, root=0)
    comm.Reduce(S1t_loc, S1t, op=MPI.SUM, root=0)
    comm.Reduce(S2t_loc, S2t, op=MPI.SUM, root=0)
    comm.Reduce(un_p_loc, un_p, op=MPI.SUM, root=0)

    if rank == 0:
        tdata = np.vstack((t_, uxt, uyt, S1t, S2t))
        np.savetxt(os.path.join(analysisfolder, "tdata.dat"), tdata)
        
        with h5py.File(os.path.join(analysisfolder, "pore_data.h5"), "w") as h5f:
            h5f.create_dataset("un", data=un_p)
            h5f.create_dataset("n", data=n_p)
            h5f.create_dataset("x", data=x_p)
            h5f.create_dataset("lp", data=l_p)
            h5f.create_dataset("t", data=t_)
            h5f.attrs.create("d", d)