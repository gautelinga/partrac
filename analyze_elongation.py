import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import h5py
from utils import Params, get_timeseries
from scipy.optimize import curve_fit


parser = argparse.ArgumentParser(description="Make elongation pdf from sheet or strip")
parser.add_argument("folder", type=str, help="Folder")
parser.add_argument("-t", type=float, default=0.0, help="t")
parser.add_argument("-bins", type=int, default=256, help="Number of bins")
parser.add_argument("--show", action="store_true", help="Show")
parser.add_argument("--terminal", action="store_true", help="Print to terminal")
parser.add_argument("-o", "--output", type=str, default="", help="Output")
parser.add_argument("-nstd", type=int, default=5, help="Number of stds to include in data")
parser.add_argument("--all", action="store_true", help="Do it on all")
parser.add_argument("--nolog", action="store_true", help="No logarithm")
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

if not args.all:
    it = np.argmin(abs(np.array(ts)-args.t))
    t = ts[it]
    print("Found", t)
    ts = [t]

filename = os.path.join(args.folder, "elongdata.dat")
elongdatafile = open(filename, "w")

print("{}\t{}\t{}\t{}".format("t", "mean", "var", "std"))
for t in ts:
    posft, grp = posf[t]
    fields = []
    with h5py.File(posft, "r") as h5f:
        if grp + "/dA" in h5f and grp + "/dA0" in h5f:
            if not args.all:
                print("Found face areas")
            dA = np.array(h5f[grp + "/dA"])[:, 0]
            dA0 = np.array(h5f[grp + "/dA0"])[:, 0]
            elong = dA/dA0
            w = dA0
        elif grp + "/dl" in h5f and grp + "/dl0" in h5f:
            if not args.all:
                print("Found edge lengths")
            dl = np.array(h5f[grp + "/dl"])[:, 0]
            dl0 = np.array(h5f[grp + "/dl0"])[:, 0]
            elong = dl/dl0
            w = dl0
        else:
            print("Does not contain this.")
            exit()

    ids = np.argsort(elong)
    elong = elong[ids]
    w = w[ids]
    w /= w.sum()
    wcum = np.cumsum(w)
    ids = np.logical_and(wcum > 0.01, wcum < 0.99)
    elong = elong[ids]
    w = w[ids]
    wcum = wcum[ids]

    if args.nolog:
        data = elong
    else:
        data = np.log(elong)
    data_mean = np.mean(data*w)/np.mean(w)
    data_var = np.mean((data-data_mean)**2*w)/np.mean(w)
    data_std = np.sqrt(data_var)
    string = "{}\t{}\t{}\t{}".format(t, data_mean, data_var, data_std)
    print(string)
    elongdatafile.write(string + "\n")

hist, bin_edges = np.histogram(data, weights=w, density=True,
                               bins=args.bins,
                               range=(data_mean-args.nstd*data_std,
                                      data_mean+args.nstd*data_std))
x = 0.5*(bin_edges[1:]+bin_edges[:-1])

if args.show:
    plt.plot(x, hist)
    plt.show()


if args.terminal:
    for xi, hi in zip(x, hist):
        print("{}\t{}".format(xi, hi))

if args.output != "":
    with open(args.output, "w") as ofile:
        for xi, hi in zip(x, hist):
            ofile.write("{}\t{}\n".format(xi, hi))
