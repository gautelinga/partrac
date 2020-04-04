import argparse
import os
import matplotlib.pyplot as plt
import numpy as np


parser = argparse.ArgumentParser(description="Plot pos")
parser.add_argument("folder", type=str, help="Folder")
args = parser.parse_args()


def read_params(folder):
    paramsfile = os.path.join(folder, "params.dat")
    params = dict()
    with open(paramsfile) as pf:
        line = pf.readline()
        cnt = 1
        while line:
            item = line.strip().split("=")
            key = item[0]
            val = item[1]
            params[key] = val
            line = pf.readline()
            cnt += 1
    return params


folder = os.path.join(args.folder, "Positions")
params = read_params(args.folder)

Lx = float(params["Lx"])
Ly = float(params["Ly"])

files = os.listdir(folder)

posf = dict()
for file in files:
    if file[:4] == "xy_t" and file[-4:] == ".pos":
        t = float(file[4:-4])
        posf[t] = os.path.join(folder, file)

imgfolder = os.path.join(args.folder, "Images")
if not os.path.exists(imgfolder):
    os.makedirs(imgfolder)

ts = list(sorted(posf.keys()))
for t in ts:
    data = np.loadtxt(posf[t])
    fig, ax = plt.subplots(figsize=(5, 10))
    ax.scatter(np.remainder(data[:, 1], Lx),
               np.remainder(data[:, 2], Ly),
               c=data[:, 0],
               marker=',', lw=0, s=0.5, cmap=plt.cm.twilight)
    plt.tick_params(
        axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        left=False,
        right=False,
        labelleft=False,
        labelbottom=False) # labels along the bottom edge are off
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(os.path.join(
        imgfolder,
        "pos_{:06d}.png".format(int(t))))
    plt.close()
