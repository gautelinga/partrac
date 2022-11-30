import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import h5py
from utils import Params


parser = argparse.ArgumentParser(description="Plot pos")
parser.add_argument("folder", type=str, help="Folder")
parser.add_argument("--skip", default=1, type=int, help="Skip")
args = parser.parse_args()


params = Params(args.folder)
t0 = params.get_tmin()
params.get("Lx", t0)

Lx = float(params.get("Lx", t0))
Ly = float(params.get("Ly", t0))

histfolder = os.path.join(args.folder, "Histograms")
files = os.listdir(histfolder)

histf = dict()
for file in files:
    if file[-5:] == ".hist":
        if file[:9] == "logdens_t":
            t = float(file[9:-5])
        elif file[:10] == "logelong_t":
            t = float(file[10:-5])
        f = os.path.join(histfolder, file)
        histf[t] = f

ts = list(sorted(histf.keys()))
print(ts)

for t in ts[::args.skip]:
    data = np.loadtxt(histf[t])
    rho = data[:, 0]
    w = data[:, 4]
    fig, ax = plt.subplots(figsize=(5, 10))

    plt.title("t={}".format(t))
    plt.hist(rho, weights=w, bins=100)
    plt.xlabel("log(rho)")
    plt.ylabel("P(log(rho))")
    plt.show()
