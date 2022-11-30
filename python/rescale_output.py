import argparse
import numpy as np
import h5py


parser = argparse.ArgumentParser(description="Rescale velocity of output")
parser.add_argument("infile", type=str, help="Output file to rescale")
parser.add_argument("outfile", type=str, help="Save file as...")
parser.add_argument("-f", "--factor", type=float, default=-1.0, help="Rescale factor")
args = parser.parse_args()


h5f_in = h5py.File(args.infile, "r")
h5f_out = h5py.File(args.outfile, "w")

for field in h5f_in:
    data = np.array(h5f_in[field])
    if field in ["u_x", "u_y", "u_z"]:
        data = args.factor*data
    h5f_out.create_dataset(field, data=data)
