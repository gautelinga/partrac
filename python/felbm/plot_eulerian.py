import argparse
import sys
sys.path.append("..")
from utils import Params, read_timestamps
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import Normalize, colorConverter, LinearSegmentedColormap

def parse_args():
    parser = argparse.ArgumentParser(description="Plot pos")
    parser.add_argument("folder", type=str, help="Folder")  
    parser.add_argument("-axis", type=int, default=0, help="Axis")
    parser.add_argument("--show", action="store_true", help="Show plot")
    parser.add_argument("--cbar", action="store_true", help="Colorbar")

    parser.add_argument("-size_x", type=float, default=5, help="Image size_y")
    parser.add_argument("-size_y", type=float, default=5, help="Image size_x")

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    #params = Params(args.folder)
    #t0 = params.get_tmin()
    #Lx = params.get("Lx", t0)

    felbm_folder = args.folder
    imgfolder = os.path.join(felbm_folder, "Images")
    if not os.path.exists(imgfolder):
        os.makedirs(imgfolder)

    timestamps = read_timestamps(os.path.join(felbm_folder, "timestamps.dat"))

    with h5py.File(os.path.join(felbm_folder, "output_is_solid.h5"), "r") as h5f:
        is_solid = np.array(h5f["is_solid"])
        nx, ny, nz = is_solid.shape
        print(nx, ny, nz)
        is_solid = is_solid.reshape((nz, ny, nx))

    is_solid = is_solid[nz // 2, :, :]
    is_fluid = np.logical_not(is_solid)
    p = np.ma.masked_where(is_solid, np.zeros_like(is_solid, dtype=float))
    x, y = np.meshgrid(np.arange(nx) + .5, np.arange(ny) + 0.5)

    figsize = (args.size_x, args.size_y)
    c_tr = colorConverter.to_rgba("blue", alpha=0.0)
    c_r = colorConverter.to_rgba("blue", alpha=0.2)
    cmap_tr = LinearSegmentedColormap.from_list("tr_cmap", [c_tr, c_r])

    for timestamp in timestamps:
        t = timestamp[0]
        h5filename = timestamp[1]
        with h5py.File(os.path.join(felbm_folder, h5filename), "r") as h5f:
            p[is_fluid] = np.array(h5f["pressure"]).reshape((nz, ny, nx))[nz//2, :, :][is_fluid]

        fig, ax = plt.subplots(figsize=figsize)

        pcm = ax.pcolormesh(x, y, p, cmap=plt.get_cmap("viridis"), shading='nearest')
        #ax.imshow(p)
        ax.set_xlim(0, nx)
        ax.set_ylim(0, ny)
        ax.set_aspect('equal')
        if args.cbar:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = fig.colorbar(pcm, ax=ax, cax=cax)
            cbar.set_label("p", rotation=270)
    
        if args.show:
            plt.show()

        plt.savefig(os.path.join(imgfolder, "p_{:06d}.png".format(int(t))))
        plt.close()
