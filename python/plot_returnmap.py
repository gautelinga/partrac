import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import h5py
from utils import Params


parser = argparse.ArgumentParser(description="Plot velocity")
parser.add_argument("folder", type=str, help="Folder")
parser.add_argument("-t0", type=float, default=0.0, help="t0")
parser.add_argument("--show", action="store_true", help="Show plot")
parser.add_argument("-tau", type=float, default=1.0, help="tau")
parser.add_argument("-n", type=int, default=1, help="number of panels")
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

#print(posf)

imgfolder = os.path.join(args.folder, "Images")
statsfolder = os.path.join(args.folder, "Statistics")
if not os.path.exists(imgfolder):
    os.makedirs(imgfolder)

ts = list(sorted(posf.keys()))
print(ts)
assert(args.t0 in ts)
tq = [args.t0]
for i in range(1, 1+args.n):
    assert(args.t0+args.tau in ts)
    tq.append(args.t0+i*args.tau)

posft, cat = posf[tq[0]]
with h5py.File(posft, "r") as h5f:
    data_0 = np.array(h5f[cat]["points"])
x0 = np.remainder(data_0[:, 0].flatten(), L[0])
y0 = np.remainder(data_0[:, 1].flatten(), L[1])

xtau = []
ytau = []
for i in range(1, 1+args.n):
    posft, cat = posf[tq[i]]
    with h5py.File(posft, "r") as h5f:
        data_tau = np.array(h5f[cat]["points"])

    xtau.append(np.remainder(data_tau[:, 0].flatten(), L[0]))
    ytau.append(np.remainder(data_tau[:, 1].flatten(), L[1]))

c = np.vstack(((x0-x0.min())/L[0],
               np.zeros_like(x0)+0.2,
               (y0-y0.min())/L[1])).T

fig, ax_ = plt.subplots(1, args.n+1, figsize=(3*(args.n+1), 6))
# fig, ax = plt.subplots(figsize=(5, 10))
ax_[0].scatter(x0, y0, c=c, marker='.', s=1)
ax_[0].set_xlim(0, L[0])
ax_[0].set_ylim(0, L[1])
ax_[0].set_aspect('equal')
ax_[0].axis("off")
mapstr = r"\mathbf{x}"
ax_[0].set_title(r"$" + mapstr + r"$")
for i in range(1, 1+args.n):
    ax_[i].scatter(xtau[i-1], ytau[i-1], c=c, marker='.', s=1)
    ax_[i].set_xlim(0, L[0])
    ax_[i].set_ylim(0, L[1])
    ax_[i].set_aspect('equal')
    ax_[i].axis("off")
    mapstr = r"\Phi_{\tau} (" + mapstr + r")"
    ax_[i].set_title(r"$" + mapstr + r"$")
plt.tight_layout()
if args.show:
    plt.show()

if False:
    np.savetxt(os.path.join(
        statsfolder, "x{}_u{}_{:06d}.dat".format(pcomp, ucomp, int(t))),
               np.vstack((x, u)).T)
    plt.savefig(os.path.join(
        imgfolder,
        "vel_{:06d}.png".format(int(t))))
    plt.close()
