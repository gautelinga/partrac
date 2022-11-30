import argparse
from utils import Params
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


def normalize(arr):
    arr_min = np.min(arr)
    return (arr-arr_min)/(np.max(arr)-arr_min)


def explode(data):
    shape_arr = np.array(data.shape)
    size = shape_arr[:3]*2 - 1
    exploded = np.zeros(np.concatenate([size, shape_arr[3:]]), dtype=data.dtype)
    exploded[::2, ::2, ::2] = data
    return exploded


def expand_coordinates(indices):
    x, y, z = indices
    x[1::2, :, :] += 1
    y[:, 1::2, :] += 1
    z[:, :, 1::2] += 1
    return x, y, z


def plot_cube(cube, angle=320):
    cube = normalize(cube)

    facecolors = cm.viridis(cube)
    facecolors[:, :, :, -1] = cube
    facecolors = explode(facecolors)

    filled = facecolors[:, :, :, -1] != 0
    x, y, z = expand_coordinates(np.indices(np.array(filled.shape) + 1))

    fig = plt.figure(figsize=(30/2.54, 30/2.54))
    ax = fig.gca(projection='3d')
    ax.view_init(30, angle)
    #ax.set_xlim(right=IMG_DIM*2)
    #ax.set_ylim(top=IMG_DIM*2)
    #ax.set_zlim(top=IMG_DIM*2)

    ax.voxels(x, y, z, filled, facecolors=facecolors)
    plt.show()


parser = argparse.ArgumentParser(description="Plot pos")
parser.add_argument("folder", type=str, help="Folder")
args = parser.parse_args()

params = Params(args.folder)

nx = int(params.get("nx", 0.0))
ny = int(params.get("ny", 0.0))
nz = int(params.get("nz", 0.0))

inside_ = np.zeros((nx, ny, nz), dtype=bool)
u_ = np.zeros((nx, ny, nz), dtype=float)

data = np.loadtxt(os.path.join(args.folder, "nodal_values.dat"))
for line in data:
    ix = int(line[0])
    iy = int(line[1])
    iz = int(line[2])
    inside = bool(line[3])
    ux, uy, uz = line[4:]
    inside_[ix, iy, iz] = inside
    u_[ix, iy, iz] = np.sqrt(ux**2 + uy**2 + uz**2)

colors = np.empty
plot_cube(u_)
