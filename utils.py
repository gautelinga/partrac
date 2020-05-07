import numpy as np
import os


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


def read_params(folder):
    paramsfiles = dict()
    for filename in os.listdir(folder):
        if "params_from_t" in filename:
            t = float(filename[13:-4])
            paramsfiles[t] = os.path.join(folder, filename)
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
