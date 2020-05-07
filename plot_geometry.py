import argparse
from utils import Params
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import h5py


parser = argparse.ArgumentParser(description="Plot pos")
parser.add_argument("folder", type=str, help="Folder")
parser.add_argument("-axis", type=int, default=0, help="Axis")
args = parser.parse_args()

with h5py.File(os.path.join(args.folder, "output_is_solid.h5"), "r") as h5f:
    data = np.array(h5f["is_solid"])
    nx, ny, nz = data.shape
    print(nx, ny, nz)
    data = np.array(data.flatten()).reshape(nz, ny, nx)

fig, ax = plt.subplots()
ax.imshow(data.mean(args.axis), cmap="Greys", origin="lower")
ax.set_aspect("equal")
plt.show()
