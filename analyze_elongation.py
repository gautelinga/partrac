import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import h5py
from utils import Params, get_timeseries
from scipy.optimize import curve_fit


parser = argparse.ArgumentParser(description="Make elongation pdf from sheet or strip")
parser.add_argument("folder", type=str, help="Folder")
parser.add_argument("-t", type=str, default="0.0", help="t")
parser.add_argument("--range", type=str, default=None, help="t")
parser.add_argument("-bins", type=int, default=100, help="Number of bins")
parser.add_argument("--show", action="store_true", help="Show")
parser.add_argument("--terminal", action="store_true", help="Print to terminal")
parser.add_argument("-o", "--output", type=str, default="", help="Output")
parser.add_argument("-nstd", type=int, default=5, help="Number of stds to include in data")
parser.add_argument("--single", action="store_true", help="Do it on a single one")
parser.add_argument("--nolog", action="store_true", help="No logarithm")
parser.add_argument("--tol", type=float, default=0.0, help="Tolerance for removing outliers")
parser.add_argument("--weights", type=str, default="dl0", help="Weights")
parser.add_argument("--skip", type=int, default=1, help="Skip")
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

ts = np.array(ts)

inn = args.t.split(":")
if len(inn) == 1:
    t0 = inn[0]
    t1 = ts[-1]
elif len(inn) == 2:
    t0, t1 = inn
else:
    print("Wrong input")
    exit()
t0 = float(t0)
t1 = float(t1)

if args.single:
    it = np.argmin(abs(np.array(ts)-t_in))
    t = ts[it]
    print("Found", t)
    ts = [t]

ts = ts[ts >= t0]
ts = ts[ts <= t1]

ax = plt.axes()
ax.set_prop_cycle('color', [plt.cm.viridis(i) for i in np.linspace(0, 1, len(ts[::args.skip]))])

ht = np.zeros((len(ts[::args.skip]), args.bins))

if True:
    filename = os.path.join(args.folder, "elongdata.dat")
    elongdatafile = open(filename, "w")

    print("{}\t{}\t{}\t{}".format("t", "mean", "var", "std"))
    for it, t in enumerate(ts[::args.skip]):
        posft, grp = posf[t]
        fields = []
        with h5py.File(posft, "r") as h5f:
            if grp + "/dA" in h5f and grp + "/dA0" in h5f:
                if args.single:
                    print("Found face areas")
                dA = np.array(h5f[grp + "/dA"])[:, 0]
                dA0 = np.array(h5f[grp + "/dA0"])[:, 0]
                elong = dA/dA0
                w = dA0
            elif grp + "/dl" in h5f and grp + "/dl0" in h5f:
                if args.single:
                    print("Found edge lengths")
                dl = np.array(h5f[grp + "/dl"])[:, 0]
                dl0 = np.array(h5f[grp + "/dl0"])[:, 0]
                elong = dl/dl0
                if args.weights == "dl":
                    w = dl
                elif args.weights == "1":
                    w = np.ones_like(dl0)
                else:
                    w = dl0
            else:
                print("Does not contain this.")
                exit()

        ids = np.argsort(elong)
        elong = elong[ids]
        w = w[ids]
        w /= w.sum()
        wcum = np.cumsum(w)
        ids = np.logical_and(wcum > args.tol, wcum < 1-args.tol)
        elong = elong[ids]
        w = w[ids]
        wcum = wcum[ids]

        if args.nolog:
            data = elong
            var = "\\rho"
        else:
            data = np.log(elong)
            var = "\\log(\\rho)"
        data_mean = np.mean(data*w)/np.mean(w)
        data_var = np.mean((data-data_mean)**2*w)/np.mean(w)
        data_std = np.sqrt(data_var)
        string = "{}\t{}\t{}\t{}".format(t, data_mean, data_var, data_std)
        print(string)
        elongdatafile.write(string + "\n")

        if args.range is None:
            brange = (data_mean-args.nstd*data_std,
                      data_mean+args.nstd*data_std)
        else:
            x0, x1 = [float(a) for a in args.range.split(":")]
            brange = (x0, x1)

        hist, bin_edges = np.histogram(data, weights=w, density=True,
                                       bins=args.bins,
                                       range=brange)
        ht[it, :] = hist
        x = 0.5*(bin_edges[1:]+bin_edges[:-1])

        if args.show:
            if not args.nolog:
                plt.plot(x, hist, label="$t={}$".format(t))
                plt.xlabel("$" + var + "$")
                plt.ylabel("$P(" + var + ")$")
            else:
                xf = (x[hist > 0])
                ff = np.log(hist[hist > 0])
                plt.plot(xf, ff)
                plt.xlabel("$" + var + "$")
                plt.ylabel("$\\log P(" + var + ")$")

        if args.terminal:
            for xi, hi in zip(x, hist):
                print("{}\t{}".format(xi, hi))

        if args.output != "":
            with open(args.output, "w") as ofile:
                for xi, hi in zip(x, hist):
                    ofile.write("{}\t{}\n".format(xi, hi))

if args.show:
    plt.show()

    plt.imshow(ht.T)
    plt.show()
