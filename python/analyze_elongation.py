import argparse
import os

import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from utils import Params, get_timeseries

parser = argparse.ArgumentParser(
    description="Make elongation pdf from sheet or strip")
parser.add_argument("folder", type=str, help="Folder")
parser.add_argument("-t", type=str, default="0.0", help="t")
parser.add_argument("--range", type=str, default=None, help="t")
parser.add_argument("-bins", type=int, default=100, help="Number of bins")
parser.add_argument("--show", action="store_true", help="Show")
parser.add_argument("--terminal",
                    action="store_true",
                    help="Print to terminal")
parser.add_argument("--output", action="store_true", help="Output")
parser.add_argument("-nstd",
                    type=int,
                    default=5,
                    help="Number of stds to include in data")
parser.add_argument("--single",
                    action="store_true",
                    help="Do it on a single one")
parser.add_argument("--nolog", action="store_true", help="No logarithm")
parser.add_argument("--tol",
                    type=float,
                    default=0.0,
                    help="Tolerance for removing outliers")
parser.add_argument("--weights", type=str, default="dl0", help="Weights")
parser.add_argument("--skip", type=int, default=1, help="Skip")
args = parser.parse_args()


def calc_moments(data, w):
    w /= w.sum()
    assert all(w >= 0)
    data_mean = np.sum(data * w)  # since w.sum() == 1
    data_var = np.sum(
        (data - data_mean)**2 *
        w)  # OK for large N, how to calculate sample mean with weights?
    data_std = np.sqrt(data_var)
    return data_mean, data_var, data_std


def calc_hist(data, w, data_mean, data_std, brange, nstd, bins):
    if brange is None:
        brange = (data_mean - args.nstd * data_std,
                  data_mean + args.nstd * data_std)
    else:
        x0, x1 = [float(a) for a in args.range.split(":")]
        brange = (x0, x1)

    hist, bin_edges = np.histogram(data,
                                   weights=w,
                                   density=True,
                                   bins=bins,
                                   range=brange)
    x = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    return x, hist


if __name__ == "__main__":
    params = Params(args.folder)
    t0 = params.get_tmin()

    ts, posf = get_timeseries(args.folder)

    possible_fields = [["u", "Vector", "Node"], ["c", "Scalar", "Node"],
                       ["p", "Scalar", "Node"], ["rho", "Scalar", "Node"],
                       ["H", "Scalar", "Node"], ["n", "Vector", "Node"],
                       ["dA", "Scalar", "Face"], ["dA0", "Scalar", "Face"],
                       ["dl", "Scalar", "Edge"], ["dl0", "Scalar", "Edge"]]

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
        it = np.argmin(abs(np.array(ts) - t_in))
        t = ts[it]
        print("Found", t)
        ts = [t]

    ts = ts[ts >= t0]
    ts = ts[ts <= t1]

    ax = plt.axes()
    ax.set_prop_cycle(
        'color',
        [plt.cm.viridis(i) for i in np.linspace(0, 1, len(ts[::args.skip]))])

    ht = np.zeros((len(ts[::args.skip]), args.bins))

    analysis_folder = os.path.join(args.folder, "Analysis")
    if not os.path.exists(analysis_folder):
        os.makedirs(analysis_folder)

    filename = os.path.join(analysis_folder, "elongdata.dat")
    elongdatafile = open(filename, "w")

    histdata = []

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
                elong = dA / dA0
                w = dA0
            elif grp + "/dl" in h5f and grp + "/dl0" in h5f:
                if args.single:
                    print("Found edge lengths")
                dl = np.array(h5f[grp + "/dl"])[:, 0]
                dl0 = np.array(h5f[grp + "/dl0"])[:, 0]
                elong = dl / dl0
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
        #ids = np.logical_and(wcum > args.tol, wcum < 1 - args.tol)
        #elong = elong[ids]
        #w = w[ids]
        #wcum = wcum[ids]

        logelong = np.log(elong)
        logelong_mean, logelong_var, logelong_std = calc_moments(logelong, w)
        elong_mean, elong_var, elong_std = calc_moments(elong, w)

        string = "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(t, elong_mean, elong_var,
                                                     elong_std, logelong_mean,
                                                     logelong_var,
                                                     logelong_std)
        print(string)
        elongdatafile.write(string + "\n")

        x_elong, hist_elong = calc_hist(elong, w, elong_mean, elong_std,
                                        args.range, args.nstd, args.bins)
        x_logelong, hist_logelong = calc_hist(logelong, w, logelong_mean,
                                              logelong_std, args.range,
                                              args.nstd, args.bins)
        ht[it, :] = hist_logelong

        if args.show:
            if not args.nolog:
                var = "\\rho"
                plt.plot(x_elong, hist_elong, label="$t={}$".format(t))
                plt.xlabel("$" + var + "$")
                plt.ylabel("$P(" + var + ")$")
            else:
                var = "\\log(\\rho)"
                ids = hist_logelong > 0
                xf = (x_logelong[ids])
                ff = np.log(hist_logelong[ids])
                plt.plot(xf, ff)
                plt.xlabel("$" + var + "$")
                plt.ylabel("$\\log P(" + var + ")$")

        if args.output:
            histdata.append((t, np.array(list(zip(x_elong, hist_elong))),
                             np.array(list(zip(x_logelong, hist_logelong)))))
    elongdatafile.close()

    if args.output:
        with h5py.File(os.path.join(analysis_folder, "histograms.h5"),
                       "w") as h5f:
            for t, hist_elong, hist_logelong in histdata:
                dset_elong = h5f.create_dataset("{}/elong".format(t),
                                                data=hist_elong)
                dset_logelong = h5f.create_dataset("{}/logelong".format(t),
                                                   data=hist_logelong)

    if args.show:
        plt.show()

        plt.imshow(ht.T)
        plt.show()
