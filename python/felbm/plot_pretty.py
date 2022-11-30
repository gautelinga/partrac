import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import h5py
from utils import Params, read_timestamps
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import Normalize, colorConverter, LinearSegmentedColormap


parser = argparse.ArgumentParser(description="Plot pos")
parser.add_argument("folder", type=str, help="Folder")
parser.add_argument("-t_min", type=float, default=0.0, help="t_min")
parser.add_argument("-t_max", type=float, default=np.inf, help="t_max")
parser.add_argument("-cmap", type=str, default="copper_r", help="colormap")
parser.add_argument("-axis", type=int, default=2, help="Projection axis")
parser.add_argument("--show", action="store_true", help="Show plot")
parser.add_argument("--export", action="store_true", help="Export")
parser.add_argument("--cbar", action="store_true", help="Colorbar")
parser.add_argument("-cmax", type=float, default=7, help="cmax")
parser.add_argument("-cmin", type=float, default=-2, help="cmin")
parser.add_argument("-size_x", type=float, default=5, help="Image size_y")
parser.add_argument("-size_y", type=float, default=5, help="Image size_x")
parser.add_argument("-pointsize", type=float, default=1.0, help="Point/dot size")
parser.add_argument("--singlephase", action="store_true", help="Single phase")
parser.add_argument("--hideobstacles", action="store_true", help="Hide obstacles")
parser.add_argument("--skip", type=int, default=1, help="Skip timesteps")
args = parser.parse_args()

params = Params(args.folder)
t0 = params.get_tmin()
params.get("Lx", t0)

felbm_folder = os.path.dirname(os.path.dirname(os.path.dirname(os.path.join(args.folder, ""))))

timestamps = read_timestamps(os.path.join(felbm_folder, "timestamps.dat"))

Lx = float(params.get("Lx", t0))
Ly = float(params.get("Ly", t0))
Lz = float(params.get("Lz", t0))
L = [Lx, Ly, Lz]
nx = int(params.get("nx", t0))
ny = int(params.get("ny", t0))
nz = int(params.get("nz", t0))

files = os.listdir(args.folder)

posf = dict()
for file in files:
    if file[:11] == "data_from_t" and file[-3:] == ".h5":
        t = float(file[11:-3])
        posft = os.path.join(args.folder, file)
        try:
            with h5py.File(posft, "r") as h5f:
                for grp in h5f:
                    posf[float(grp)] = (posft, grp)
        except:
            pass

imgfolder = os.path.join(args.folder, "Images")
if not os.path.exists(imgfolder):
    os.makedirs(imgfolder)

posfolder = os.path.join(args.folder, "Positions")
if not os.path.exists(posfolder):
    os.makedirs(posfolder)

ts = []
for t in list(sorted(posf.keys())):
    if t >= args.t_min and t <= args.t_max:
        ts.append(t)
ts = ts[::args.skip]

proj_axis = [[1, 2],
             [2, 0],
             [0, 1]]
pax = proj_axis[args.axis]

cmap = plt.cm.get_cmap(args.cmap)

with h5py.File(os.path.join(felbm_folder, "output_is_solid.h5"), "r") as h5f:
    is_solid = np.array(h5f["is_solid"]).reshape((nz, ny, nx))

is_solid = is_solid[nz//2, :, :]
x, y = np.meshgrid(np.arange(nx), np.arange(ny))

figsize = (args.size_x, args.size_y)

timestamp_entry = 0
while timestamps[timestamp_entry+1][0] < ts[0]:
    timestamp_entry += 1

t_prev = timestamps[timestamp_entry][0]
t_next = timestamps[timestamp_entry+1][0]

with h5py.File(os.path.join(felbm_folder, timestamps[timestamp_entry][1]), "r") as h5f:
    rho_prev = np.array(h5f["density"]).reshape((nz, ny, nx))[nz//2, :, :]
with h5py.File(os.path.join(felbm_folder, timestamps[timestamp_entry+1][1]), "r") as h5f:
    rho_next = np.array(h5f["density"]).reshape((nz, ny, nx))[nz//2, :, :]

c_tg = colorConverter.to_rgba("#e6e6e6", alpha=0.0)
c_g = colorConverter.to_rgba("#e6e6e6", alpha=1.0)
c_tr = colorConverter.to_rgba("blue", alpha=0.0)
c_r = colorConverter.to_rgba("blue", alpha=0.2)
cmap_tg = LinearSegmentedColormap.from_list("tg_cmap", [c_tg, c_g])
cmap_tr = LinearSegmentedColormap.from_list("tr_cmap", [c_tr, c_r])

# hardcoded this
rho_a = 1.05
rho_b = 0.9
rho_mid = 0.5*(rho_a+rho_b)
levels = [rho_b, rho_mid, rho_a]

for t in ts:
    if bool(len(timestamps) > timestamp_entry+2
            and timestamps[timestamp_entry+1][0] < t):
        while timestamps[timestamp_entry+1][0] < t:
            timestamp_entry += 1
        rho_prev[:, :] = rho_next[:, :]
        t_prev = t_next
        t_next = timestamps[timestamp_entry+1][0]
        filename = os.path.join(felbm_folder,
                                timestamps[timestamp_entry+1][1])
        # print(filename)
        with h5py.File(filename, "r") as h5f:
            rho_next[:, :] = np.array(h5f["density"]).reshape((nz, ny, nx))[nz//2, :, :]

    posft, grp = posf[t]
    with h5py.File(posft, "r") as h5f:
        pos = np.array(h5f[grp + "/points"])
        elong = np.array(h5f[grp + "/e"])

    fig, ax = plt.subplots(figsize=figsize)
    x1 = np.remainder(pos[:, pax[0]], L[pax[0]])
    x2 = np.remainder(pos[:, pax[1]], L[pax[1]])

    c = np.log(elong[:, 0])
    label = "$\mathrm{log}(\delta \ell/\delta \ell_0)$"

    alpha_t = (t-t_prev)/(t_next-t_prev)
    rho = alpha_t*rho_next + (1.0-alpha_t)*rho_prev

    print(t, t_prev, t_next, alpha_t)
    
    rho[rho > 1.9] = rho_mid

    if not args.singlephase:
        ax.contourf(x, y, rho, levels=levels, cmap=cmap_tr)
    if not args.hideobstacles:
        ax.contourf(x, y, is_solid, cmap=cmap_tg)
    #ax.contour(x, y, is_solid, [0.9], colors='grey', linewidths=0.5)
    p = ax.scatter(x1, x2,
                   c=c, vmin=args.cmin, vmax=args.cmax,
                   marker=',', lw=0, s=args.pointsize, cmap=cmap)
    if args.cbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = fig.colorbar(p, ax=ax, cax=cax)
        cbar.set_label(label, rotation=270)

    plt.tick_params(
        axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        left=False,
        right=False,
        labelleft=False,
        labelbottom=False)  # labels along the bottom edge are off
    ax.set_xlim(0, L[pax[0]])
    ax.set_ylim(0, L[pax[1]])
    ax.set_aspect('equal')
    plt.tight_layout()
    if args.show:
        plt.show()
    if args.export:
        order = np.argsort(x1)
        np.savetxt(os.path.join(
            posfolder,
            "pos_{:06}.pos".format(int(t))),
                   np.vstack((x1[order], x2[order], c[order])).T)
    plt.savefig(os.path.join(
        imgfolder,
        "pos_{:06d}.png".format(int(t))))
    plt.close()
