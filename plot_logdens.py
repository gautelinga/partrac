import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import h5py
from utils import Params


parser = argparse.ArgumentParser(description="Plot pos")
parser.add_argument("folder", type=str, help="Folder")
args = parser.parse_args()


params = Params(args.folder)
t0 = params.get_tmin()
params.get("Lx", t0)

Lx = float(params.get("Lx", t0))
Ly = float(params.get("Ly", t0))

histfolder = os.path.join(args.folder, "Histograms")
files = os.listdir(histfolder)

print(files)

histf = dict()
for file in files:
    if file[:9] == "logdens_t" and file[-5:] == ".hist":
        t = float(file[9:-5])
        f = os.path.join(histfolder, file)
        histf[t] = f

ts = list(sorted(histf.keys()))
for t in ts:
    rho = np.loadtxt(histf[t])
    print(rho)
    fig, ax = plt.subplots(figsize=(5, 10))

    plt.hist(rho, bins=256)
    plt.show()
