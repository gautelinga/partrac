import numpy as np
import os
import h5py


class Params:
    def __init__(self, folder):
        self.params_dict = read_params(folder)
        self.ts = np.array(sorted(self.params_dict.keys()))

    def get(self, key, t):
        i = np.where(self.ts <= t)[0][-1]
        tkey = self.ts[i]
        return self.params_dict[tkey][key]

    def get_tmin(self):
        return self.ts[-1]


def find_params(folder):
    paramsfiles = dict()
    for filename in os.listdir(folder):
        if "params_from_t" in filename:
            t = float(filename[13:-4])
            paramsfiles[t] = os.path.join(folder, filename)
    return paramsfiles


def parse_params(paramsfiles):
    params = dict()
    for t, paramsfile in paramsfiles.items():
        params[t] = dict()
        with open(paramsfile) as pf:
            line = pf.readline()
            cnt = 1
            while line:
                item = line.strip().split("=")
                key = item[0]
                val = item[1]
                params[t][key] = val
                line = pf.readline()
                cnt += 1
    return params


def read_params(folder):
    paramsfiles = find_params(folder)
    return parse_params(paramsfiles)


def read_timestamps(infile):
    timestamps = []
    with open(infile, "r") as tf:
        for line in tf:
            line = line.strip()
            tstr, fname = line.split("\t")
            timestamps.append((float(tstr), fname))
    return timestamps


def get_timeseries(folder, t_min=-np.inf, t_max=np.inf):
    files = os.listdir(folder)

    posf = dict()
    for file in files:
        if file[:11] == "data_from_t" and file[-3:] == ".h5":
            t = float(file[11:-3])
            posft = os.path.join(folder, file)
            try:
                with h5py.File(posft, "r") as h5f:
                    for grp in h5f:
                        posf[float(grp)] = (posft, grp)
            except:
                pass

    ts = []
    for t in list(sorted(posf.keys())):
        if t >= t_min and t <= t_max:
            ts.append(t)

    return ts, posf

def get_folders(folder):
    folders = []
    paramsfiles = find_params(folder)

    if len(paramsfiles) == 0:
        subfolders = [] 
        for a in os.listdir(folder):
            if a.isdigit():
                subfolders.append(a)
        subfolders = sorted(subfolders)
        for subfolder in subfolders:
            fullpath = os.path.join(folder, subfolder)
            paramsfiles = find_params(fullpath)
            folders.append(fullpath)
    return folders

def node_logelong(grp):
    # Mean log(dl/dl0) over each node's edges; None without edges
    if "edges" not in grp:
        return None
    edges = np.array(grp["edges"], dtype=int).flatten()
    if "logelong" in grp:
        le = np.array(grp["logelong"])[:, 0]
    else:
        le = np.log(np.array(grp["dl"])[:, 0] / np.array(grp["dl0"])[:, 0])
    n = len(grp["points"])
    le_sum = np.bincount(edges, weights=np.repeat(le, 2), minlength=n)
    count = np.bincount(edges, minlength=n)
    with np.errstate(invalid="ignore"):
        return le_sum / count


def run_folder(folder):
    # The newest numbered run when folder holds runs, not params
    if len(find_params(folder)) > 0:
        return folder
    runs = [a for a in os.listdir(folder)
            if a.isdigit() and os.path.isdir(os.path.join(folder, a))]
    if len(runs) == 0:
        return folder
    folder = os.path.join(folder, max(runs, key=int))
    print("Using run folder", folder)
    return folder
