import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import h5py
from utils import Params


parser = argparse.ArgumentParser(description="Plot velocity")
parser.add_argument("folder", type=str, help="Folder")
parser.add_argument("-t", type=float, default=0.0, help="t")
parser.add_argument("-cmap", type=str, default="parula", help="colormap")
parser.add_argument("-pcomp", type=int, default=2, help="Spatial comp")
parser.add_argument("-ucomp", type=int, default=2, help="Velocity comp")
parser.add_argument("--show", action="store_true", help="Show plot")
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
                for cat in h5f:
                    posf[float(cat)] = (posft, cat)
        except:
            pass

imgfolder = os.path.join(args.folder, "Images")
statsfolder = os.path.join(args.folder, "Statistics")
if not os.path.exists(imgfolder):
    os.makedirs(imgfolder)

ts = list(sorted(posf.keys()))
t_dist = []
for t in ts:
    t_dist.append(abs(t-args.t))
tq = ts[np.argmin(t_dist)]

pcomp = args.pcomp
ucomp = args.ucomp

if args.cmap == "parula":
    cmap = plt.cm.viridis
elif args.cmap == "twilight":
    cmap = plt.cm.twilight

for t in [tq]:
    posft, cat = posf[t]
    with h5py.File(posft, "r") as h5f:
        data = np.array(h5f[cat])

    fig, ax = plt.subplots(figsize=(5, 10))
    x = np.remainder(data[:, pcomp], L[pcomp])
    u = data[:, 3+ucomp]
    ax.plot(x, u, marker=',', lw=1)
    ax.set_xlim(0, L[pcomp])
    # ax.set_ylim(0, L[pax[1]])
    # ax.set_aspect('equal')
    plt.tight_layout()
    if args.show:
        plt.show()
    np.savetxt(os.path.join(
        statsfolder, "x{}_u{}_{:06d}.dat".format(pcomp, ucomp, int(t))),
               np.vstack((x, u)).T)
    plt.savefig(os.path.join(
        imgfolder,
        "vel_{:06d}.png".format(int(t))))
    plt.close()
