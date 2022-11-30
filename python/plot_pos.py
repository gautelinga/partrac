import argparse
import os

import h5py
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

from utils import Params

parser = argparse.ArgumentParser(description="Plot pos")
parser.add_argument("folder", type=str, help="Folder")
parser.add_argument("-t_min", type=float, default=0.0, help="t_min")
parser.add_argument("-t_max", type=float, default=np.inf, help="t_max")
parser.add_argument("-cmap", type=str, default="parula", help="colormap")
parser.add_argument("-axis", type=int, default=2, help="Projection axis")
parser.add_argument("--show", action="store_true", help="Show plot")
parser.add_argument("--export", action="store_true", help="Export")
parser.add_argument("--elong", action="store_true", help="Plot elong")
parser.add_argument("--nocol", action="store_true", help="Plot no color")
parser.add_argument("--cbar", action="store_true", help="Colorbar")
parser.add_argument("-cmax", type=float, default=None, help="cmax")
parser.add_argument("-cmin", type=float, default=None, help="cmin")
parser.add_argument("-size_x", type=float, default=5, help="Image size x")
parser.add_argument("-size_y", type=float, default=10, help="Image size y")
parser.add_argument("-x_shift", type=float, default=0., help="Shift x")
parser.add_argument("-pointsize",
                    type=float,
                    default=1.0,
                    help="Point/dot size")
args = parser.parse_args()

params = Params(args.folder)
t0 = params.get_tmin()
params.get("Lx", t0)

Lx = float(params.get("Lx", t0))
Ly = float(params.get("Ly", t0))
Lz = float(params.get("Lz", t0))
L = [Lx, Ly, Lz]

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

proj_axis = [[1, 2], [2, 0], [0, 1]]
pax = proj_axis[args.axis]

cmap = plt.cm.get_cmap(args.cmap)

for t in ts:
    posft, grp = posf[t]
    with h5py.File(posft, "r") as h5f:
        pos = np.array(h5f[grp + "/points"])
        if args.nocol:
            col = np.zeros_like(pos)
        elif args.elong:
            elong = np.array(h5f[grp + "/e"])
        else:
            col = np.array(h5f[grp + "/c"])

    eps = 0
    #if t == ts[0]:
    #    eps = 1e-2*np.random.rand(len(pos[:, 1]))

    fig, ax = plt.subplots(figsize=(args.size_x, args.size_y))
    x1 = np.remainder(pos[:, pax[0]] - args.x_shift, L[pax[0]])
    x2 = np.remainder(pos[:, pax[1]], L[pax[1]])
    if args.elong:
        c = np.log(elong[:, 0])
        label = "$\mathrm{log}(\delta \ell/\delta \ell_0)$"
    else:
        c = col[:, 0]
        label = "Color"

    if args.cmin is None:
        cmin = c.min()
    else:
        cmin = args.cmin
    if args.cmax is None:
        cmax = c.max()
    else:
        cmax = args.cmax

    p = ax.scatter(x1,
                   x2 + eps,
                   c=c,
                   vmin=cmin,
                   vmax=cmax,
                   marker=',',
                   lw=0,
                   s=args.pointsize,
                   cmap=cmap)
    if args.cbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.15)
        #v1 = np.linspace(cmin, cmax, 8, endpoint=True)
        cbar = fig.colorbar(p, ax=ax, cax=cax)  # , ticks=v1)
        cbar.set_label(label, rotation=270, labelpad=15.0)
        #cbar.ax.set_yticklabels([""] + ["{:4.2f}".format(i)
        #                                for i in v1[1:-1]] + [""],
        #                        fontsize='7')

    plt.tick_params(
        axis='both',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        left=False,
        right=False,
        labelleft=False,
        labelbottom=False)  # labels along the bottom edge are off
    ax.set_xlim(0, L[pax[0]])
    ax.set_ylim(0, L[pax[1]])
    ax.set_aspect('equal')
    #plt.tight_layout()
    if args.show:
        plt.show()
    if args.export:
        order = np.argsort(x1)
        np.savetxt(os.path.join(posfolder, "pos_{:06.6f}.pos".format(t)),
                   np.vstack((x1[order], x2[order], c[order])).T)
    plt.savefig(os.path.join(imgfolder, "pos_{:06.6f}.png".format(t)))
    plt.close()
