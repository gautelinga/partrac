import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import h5py
from utils import Params, get_timeseries


parser = argparse.ArgumentParser(description="Make elongation pdf from sheet or strip")
parser.add_argument("folder", type=str, help="Folder")
parser.add_argument("-t", type=float, default=0.0, help="t")
parser.add_argument("-bins", type=int, default=256, help="Number of bins")
parser.add_argument("--show", action="store_true", help="Show")
parser.add_argument("--terminal", action="store_true", help="Print to terminal")
args = parser.parse_args()

params = Params(args.folder)
t0 = params.get_tmin()

ts, posf = get_timeseries(args.folder)

possible_fields = [["u", "Vector", "Node"],
                   ["c", "Scalar", "Node"],
                   ["p", "Scalar", "Node"],
                   ["rho", "Scalar", "Node"],
                   ["H", "Scalar", "Node"],
                   ["n", "Vector", "Node"],
                   ["dA", "Scalar", "Face"],
                   ["dA0", "Scalar", "Face"],
                   ["dl", "Scalar", "Edge"],
                   ["dl0", "Scalar", "Edge"]]

it = np.argmin(abs(np.array(ts)-args.t))

t = ts[it]

print("Found", t)

posft, grp = posf[t]
fields = []
with h5py.File(posft, "r") as h5f:
    if grp + "/dA" in h5f and grp + "/dA0" in h5f:
        print("Found face areas")
        dA = np.array(h5f[grp + "/dA"])[:, 0]
        dA0 = np.array(h5f[grp + "/dA0"])[:, 0]
        elong = dA/dA0
        w = dA0
    elif grp + "/dl" in h5f and grp + "/dl0" in h5f:
        print("Found edge lengths")
        dl = np.array(h5f[grp + "/dl"])[:, 0]
        dl0 = np.array(h5f[grp + "/dl0"])[:, 0]
        elong = dl/dl0
        w = dl0
    else:
        print("Does not contain this.")
        exit()

hist, bin_edges = np.histogram(np.log(elong), weights=w, density=True, bins=args.bins)
x = 0.5*(bin_edges[1:]+bin_edges[:-1])

if args.show:
    plt.plot(x, hist)
    plt.show()

if args.terminal:
    for xi, hi in zip(x, hist):
        print("{}\t{}".format(xi, hi))

