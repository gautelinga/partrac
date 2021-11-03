import argparse
import os

import h5py
import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description="Try to collapse pdfs")
    parser.add_argument("folder", type=str, help="Folder")
    parser.add_argument("--show", action="store_true", help="Show")
    parser.add_argument("--lam",
                        type=float,
                        default=1.0,
                        help="Guessed Lyapunov exponent")
    parser.add_argument("--skip", type=int, default=1, help="Skip plots")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()

    analysis_folder = os.path.join(args.folder, "Analysis")

    with h5py.File(os.path.join(analysis_folder, "histograms.h5"), "r") as h5f:
        ts = sorted([float(a) for a in h5f])[1::args.skip]

        plt.rcParams["axes.prop_cycle"] = plt.cycler(
            "color", plt.cm.viridis(np.linspace(0, 1, len(ts))))

        for t in ts:
            d = h5f["{}".format(t)]["logelong"]
            mu = args.lam * t
            x = d[:, 0]
            h = d[:, 1]
            plt.plot(x / mu, h * mu)

    plt.show()
