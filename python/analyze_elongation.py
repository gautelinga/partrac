import argparse
import os
from sunau import AUDIO_FILE_ENCODING_DOUBLE

import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy import signal
from statsmodels.nonparametric.smoothers_lowess import lowess

from utils import Params, find_params, get_timeseries, read_params, get_folders

parser = argparse.ArgumentParser(
    description="Make elongation pdf from sheet or strip")
parser.add_argument("folder", type=str, help="Folder")
parser.add_argument("-t", type=str, default="0.0", help="t")
parser.add_argument("--range", type=str, default=None, help="t")
parser.add_argument("-bins", type=int, default=100, help="Number of bins")
parser.add_argument("--show", action="store_true", help="Show")
parser.add_argument("--save", action="store_true", help="Save")
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


def gaussian(x, mu,sig):
    return 1./(np.sqrt(2.*np.pi)*sig) * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


def lognormal(x, mu, sig):
    return 1./(np.sqrt(2.*np.pi)*sig*x) * np.exp(-np.power(np.log(x) - mu, 2.) / (2 * np.power(sig, 2.)))


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
        brange = (data_mean - nstd * data_std,
                  data_mean + nstd * data_std)
    elif isinstance(brange, str):
        x0, x1 = [float(a) for a in brange.split(":")]
        brange = (x0, x1)

    hist, bin_edges = np.histogram(data,
                                   weights=w,
                                   density=True,
                                   bins=bins,
                                   range=brange)
    x = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    return x, hist


if __name__ == "__main__":
    possible_fields = [["u", "Vector", "Node"],  ["c", "Scalar", "Node"],
                       ["p", "Scalar", "Node"],  ["rho", "Scalar", "Node"],
                       ["H", "Scalar", "Node"],  ["n", "Vector", "Node"],
                       ["dA", "Scalar", "Face"], ["dA0", "Scalar", "Face"],
                       ["dl", "Scalar", "Edge"], ["dl0", "Scalar", "Edge"]]

    folders = get_folders(args.folder)
    nfld = len(folders)
    if nfld == 0:
      folders = [args.folder]

    analysis_folder = os.path.join(args.folder, "Analysis")
    if not os.path.exists(analysis_folder):
        os.makedirs(analysis_folder)

    images_folder = os.path.join(args.folder, "Images")
    if not os.path.exists(images_folder):
        os.makedirs(images_folder)

    elongfilename = os.path.join(analysis_folder, "elongdata.dat")
    #elongdatafile = open(elongfilename, "w")

    #params_ = []
    ts_ = []
    posf_ = []
    for ifolder, folder in enumerate(folders):
        params = Params(folder)
        t0 = params.get_tmin()

        ts, posf = get_timeseries(folder)

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

        #if args.single:
        #    it = np.argmin(abs(np.array(ts) - t_in))
        #    t = ts[it]
        #    print("Found", t)
        #    ts = [t]

        ts = ts[ts >= t0]
        ts = ts[ts <= t1]

        ts_.append(ts)
        posf_.append(posf)

        #ax = plt.axes()
        #ax.set_prop_cycle(
        #    'color',
        #    [plt.cm.viridis(i) for i in np.linspace(0, 1, len(ts[::args.skip]))])

    ht = np.zeros((len(ts_[0][::args.skip]), args.bins))

    histdata = []
    for t in ts_[0]:
        for posf in posf_:
            assert(t in posf)

    c0 = None

    #print("{}\t{}\t{}\t{}".format("t", "mean", "var", "std"))
    t_ = ts_[0][::args.skip]
    rho_mean_ = np.zeros_like(t_)
    #rho_mean2 = np.zeros_like(t_)
    rho_var_ = np.zeros_like(t_)
    logrho_mean_ = np.zeros_like(t_)
    logrho_var_ = np.zeros_like(t_)
    for it, t in enumerate(t_):
        if it % 10:
            print("t = {} \t\t({}/{})".format(t, it, len(t_)))

        c_ = []
        dl_ = []
        dl0_ = []

        fig, (ax1, ax2) = plt.subplots(1, 2)
        for posf in posf_:
            posft, grp = posf[t]

            with h5py.File(posft, "r") as h5f:
                c = np.array(h5f[grp + "/c"][:, 0])
                dl = np.array(h5f[grp + "/dl"][:, 0])
                dl0 = np.array(h5f[grp + "/dl0"][:, 0])
                edges = np.array(h5f[grp + "/edges"], dtype=int)
                x = np.array(h5f[grp + "/points"])

                c_.append(c)
                dl_.append(dl)
                dl0_.append(dl0)

        """
            print("num nodes:", len(c))
            print("num edges:", len(dl))

            next_v = [None for _ in range(len(c))]
            prev_v = [None for _ in range(len(c))]

            for el in edges:
                v1, v2 = el[0], el[1]
                if c[v1] < c[v2]:
                    next_v[v1] = v2
                    prev_v[v2] = v1
                else:
                    next_v[v2] = v1
                    prev_v[v1] = v2

            vs = np.argsort(c)
            splits = []
            for i, v in enumerate(vs):
                if next_v[v] is None:
                    splits.append(i+1)
            lines = np.split(vs, splits)[:-1]

            for line in lines[:3]:
                ds = np.linalg.norm(x[line[1:], :]-x[line[:-1], :], axis=1)
                dc = c[line[1:]]-c[line[:-1]]
                if it==0:
                    c0 = dc.sum()
                dc /= c0
                ss = np.zeros(len(dc)+1)
                ss[1:] = np.cumsum(ds)
                #cc = c[line[:]]-c[line[0]]
                cc = np.zeros_like(ss)
                cc[1:] = np.cumsum(dc)
                #plt.plot(x[line, 0], x[line, 1])
                ax1.plot(cc, ss) # np.log(ds/dc))
                wlen =  2*(len(ss)//2)-1
                #dsdc = np.gradient(ss, cc) #signal.savgol_filter(ss, window_length=min(5, wlen), polyorder=min(2, wlen-1), deriv=1)
                lowess_tight = lowess(ss, cc, frac=0.2)
                #print(lowess_tight)
                #ax2.plot(lowess_tight[:, 0], lowess_tight[:, 1])
                dsdc = np.gradient(lowess_tight[:, 0], lowess_tight[:, 1])
                ax2.plot(cc, dsdc)
        plt.show()
        """

        dl = np.hstack(dl_)
        dl0 = np.hstack(dl0_)
        rho = dl/dl0
        logrho = np.log(rho)
        w = dl0/dl0.sum()
        #print(dl.sum(), dl0.sum())
        rho_mean = dl.sum()/dl0.sum()
        rho_mean_[it] = rho_mean
        #rho_mean2[it] = np.sum(rho * w)
        rho_var_[it] = np.sum((rho - rho_mean)**2 * w)
        logrho_mean = np.sum(logrho * w)
        logrho_var = np.sum((logrho - logrho_mean)**2 * w)
        logrho_mean_[it] = logrho_mean
        logrho_var_[it] = logrho_var
        logrho_std = np.sqrt(logrho_var)

        x_logrho, hist_logrho = calc_hist(logrho, w, logrho_mean, logrho_std,
                                        (logrho_mean - args.nstd*logrho_std, logrho_mean + args.nstd*logrho_std), args.nstd, args.bins)
        x_rho, hist_rho = calc_hist(rho, w, rho_mean, None,
                                    (np.exp(logrho_mean - args.nstd*logrho_std),
                                    np.exp(2*logrho_mean + 0*args.nstd*logrho_std)), args.nstd, args.bins)

        if args.show or args.save:
            fig, (ax1, ax2) = plt.subplots(1, 2)

            var = "\\rho"
            ax1.plot(x_rho, hist_rho, label="$t={}$".format(t))
            ax1.set_xlabel("$" + var + "$")
            ax1.set_ylabel("$P(" + var + ")$")
            ax1.set_yscale("log")
            xx = np.linspace(0, x_rho[-1], 1000)[1:]
            ax1.plot(xx, lognormal(xx, logrho_mean, logrho_std))

            var = "\\log(\\rho)"
            ax2.plot(x_logrho, hist_logrho, 'o')
            ax2.set_xlabel("$" + var + "$")
            ax2.set_ylabel("$P(" + var + ")$")
            nn = args.nstd
            xx = np.linspace(logrho_mean-nn*logrho_std, logrho_mean+nn*logrho_std, 1000)
            ax2.plot(xx, gaussian(xx, logrho_mean, logrho_std))

            if args.save:
                plt.savefig(os.path.join(images_folder, "rho_pdfs_t{}.png".format(t)))
            if args.show:
                plt.show()
            plt.close()

    elongdata = np.vstack((t_, np.log(rho_mean_), np.log(rho_var_), np.log(rho_var_)/2, logrho_mean_, logrho_var_, np.sqrt(logrho_var_))).T
    print(elongdata.shape)
    np.savetxt(elongfilename, elongdata)

    if args.show or args.save:
        fig, ax = plt.subplots(1, 1)
        tau = t_[1:] - t_[0]
        ax.plot(tau, np.log(rho_mean_[1:]), label='log(<rho>)')
        ax.plot(tau, np.log(rho_var_[1:]), label='log(Var(rho))')
        ax.plot(tau, logrho_mean_[1:], label='<log(rho)>')
        ax.plot(tau, logrho_var_[1:], label='Var(log(rho))')
        plt.legend()
        #ax.set_yscale("log")
        if args.save:
            plt.savefig(os.path.join(images_folder, "elong_t.png".format(t)))
        if args.show:
            plt.show()

        #edges_set = set([tuple(el) for el in edges])
        #print(edges_set)

        #exit()

if False:
    if False:
        if False:
            with h5py.File(posft, "r") as h5f:
                if grp + "/dA" in h5f and grp + "/dA0" in h5f:
                    if args.single:
                        print("Found face areas")
                    dA = np.array(h5f[grp + "/dA"])[:, 0]
                    dA0 = np.array(h5f[grp + "/dA0"])[:, 0]
                    elong = dA / dA0
                    ids = elong != 1.0
                    elong = elong[ids]
                    w = dA0[ids]
                elif grp + "/dl" in h5f and grp + "/dl0" in h5f:
                    if args.single:
                        print("Found edge lengths")
                    dl = np.array(h5f[grp + "/dl"])[:, 0]
                    dl0 = np.array(h5f[grp + "/dl0"])[:, 0]
                    elong = dl / dl0
                    ids = (abs(elong) - 1.0) > 1e-9
                    if args.weights == "dl":
                        w = dl
                    elif args.weights == "1":
                        w = np.ones_like(dl0)
                    else:
                        w = dl0
                    elong = elong[ids]
                    w = w[ids]
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

        string = "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(t, 
                                                     np.log(elong_mean), np.log(elong_var), np.log(elong_std),
                                                     logelong_mean, logelong_var, logelong_std)
        print(string)
        elongdatafile.write(string + "\n")

        x_elong, hist_elong = calc_hist(elong, w, elong_mean, elong_std,
                                        (np.exp(logelong_mean - args.nstd*logelong_std), np.exp(logelong_mean + args.nstd*logelong_std)), args.nstd, args.bins)
        x_logelong, hist_logelong = calc_hist(logelong, w, logelong_mean,
                                              logelong_std, None,
                                              args.nstd, args.bins)
        ht[it, :] = hist_logelong

        if args.show:
            var = "\\rho"
            
            fig, (ax1, ax2) = plt.subplots(1, 2)
            ax1.plot(x_elong, hist_elong, label="$t={}$".format(t))
            ax1.set_xlabel("$" + var + "$")
            ax1.set_ylabel("$P(" + var + ")$")
            #ax1.set_yscale("log")
            xx = np.linspace(0, x_elong[-1], 1000)[1:]
            ax1.plot(xx, lognormal(xx, logelong_mean, logelong_std))

            var = "\\log(\\rho)"
            ids = hist_logelong > 0
            xf = x_logelong[ids]
            ff = hist_logelong[ids]
            ax2.plot(xf, ff)
            ax2.set_xlabel("$" + var + "$")
            ax2.set_ylabel("$P(" + var + ")$")
            nn = 5.
            xx = np.linspace(logelong_mean-nn*logelong_std, logelong_mean+nn*logelong_std, 1000)
            ax2.plot(xx, gaussian(xx, logelong_mean, logelong_std))
            plt.show()

        if args.output:
            histdata.append((t, np.array(list(zip(x_elong, hist_elong))), np.array(list(zip(x_logelong, hist_logelong)))))
    elongdatafile.close()

    if args.output:
        with h5py.File(os.path.join(analysis_folder, "histograms.h5"), "w") as h5f:
            for t, hist_elong, hist_logelong in histdata:
                dset_elong = h5f.create_dataset("{}/elong".format(t),
                                                data=hist_elong)
                dset_logelong = h5f.create_dataset("{}/logelong".format(t),
                                                   data=hist_logelong)

    if args.show:
        plt.show()

        plt.imshow(ht.T[:, 1:])
        plt.show()
