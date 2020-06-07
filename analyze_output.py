import argparse
from utils import Params
import numpy as np
import os
import h5py


parser = argparse.ArgumentParser(description="Analyze timestamps")
parser.add_argument("infile", type=str, help="Timestamps file")
args = parser.parse_args()

timestamps = []
with open(args.infile, "r") as tf:
    for line in tf:
        line = line.strip()
        tstr, fname = line.split("\t")
        timestamps.append((float(tstr), fname))

folder = os.path.dirname(args.infile)

with h5py.File(os.path.join(folder, "output_is_solid.h5"), "r") as h5f:
    is_solid = np.array(h5f["is_solid"], dtype=bool)
    is_fluid_flat = np.logical_not(is_solid).flatten()
    print("porosity =", is_solid.mean())

u_mean = np.zeros((len(timestamps), 3))
#su_mean = np.zeros_like(u_mean)
for it, (t, h5fname) in enumerate(timestamps):
    filename = os.path.join(folder, h5fname)
    with h5py.File(filename, "r") as h5f:
        u_x = np.array(h5f["u_x"]).flatten()
        u_y = np.array(h5f["u_y"]).flatten()
        u_z = np.array(h5f["u_z"]).flatten()

    u_mean[it, 0] = u_x[is_fluid_flat].mean()
    u_mean[it, 1] = u_y[is_fluid_flat].mean()
    u_mean[it, 2] = u_z[is_fluid_flat].mean()

print("mean velocity =", *u_mean.mean(0))
