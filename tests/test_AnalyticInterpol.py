import json
import os
import subprocess
import sys
import tempfile

import h5py
import numpy as np
import pytest


def make_temp_case(expr_params):
    tmpdir = tempfile.mkdtemp()
    prmfile = open(os.path.join(tmpdir, "expr_params.dat"), "w+")
    for prm, val in expr_params.items():
        prmfile.write("{}={}\n".format(prm, val))
    return tmpdir


def destroy_temp_case(tmpdir):
    for root, dirs, files in os.walk(tmpdir, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))


def make_batchelor_case():
    expr_params = dict(t_min=0.0,
                       t_max=10000000.0,
                       Lx=10.0,
                       Ly=10.0,
                       Lz=10.0,
                       expression="batchelor_vortex",
                       p_inf=1.0,
                       u0=1.0,
                       R1=1.0,
                       R2=1.0,
                       q=1.0,
                       rho=1.0,
                       x0=5.0,
                       y0=5.0,
                       z0=5.0,
                       x_min=0.0,
                       y_min=0.0,
                       z_min=0.0,
                       x_max=10.0,
                       y_max=10.0,
                       z_max=10.0)

    tmpdir = make_temp_case(expr_params)

    return tmpdir


@pytest.mark.parametrize("dt", [0.1, 0.2, 0.4])
def test_batchelor(dt):
    tmpdir = make_batchelor_case()
    #with open("{}/expr_params.dat".format(tmpdir), "r") as ofile:
    #    print(ofile.read())

    cmd = "partrac {}/expr_params.dat mode=analytic dt={} t0=0 T=10.0 stat_intv=0.1 dump_intv=0.1 dump_chunk_size=1000 minimal_output=true init_mode=uniform_x Dm=0"
    d = subprocess.check_output(cmd.format(tmpdir, dt),
                                shell=True).decode("utf-8").split("\n")
    #print(d)
    #assert d[0] == "AnalyticInterpol initiated."
    #assert d[1] == "BatchelorVortex selected"

    rwpath = os.path.join(tmpdir, "RandomWalkers")
    rwdir = os.listdir(rwpath)
    #print(rwdir)
    #rwkey1 = "Dm0.0000000e+00_dt{:1.7e}_Nrw100_seed0".format(dt)
    #print(rwkey1)
    #assert rwkey1 in rwdir
    rwkey1 = rwdir[0]

    #print(os.listdir(os.path.join(rwpath, rwkey1)))

    with h5py.File(os.path.join(rwpath, rwkey1,
                                "data_from_t0.000000.h5")) as h5f:
        for t in h5f:
            pass

        l = sum(h5f[t]["dl"])[0]
        l0 = sum(h5f[t]["dl0"])[0]

    #print(__file__, "test_batchelor", dt, [rwkey1, l, l0])
    #compare_reference(__file__, "test_batchelor", dt, [rwkey1, l, l0])
    data_ref = {
        0.1: {
            "l": 10.094533234085691,
            "l0": 10
        },
        0.2: {
            "l": 10.09425003346005,
            "l0": 10
        },
        0.4: {
            "l": 10.093700231880497,
            "l0": 10
        }
    }
    tol = 1e-5

    data_curr = {"l": l, "l0": l0}
    for key, val in data_curr.items():
        val_ref = data_ref[dt][key]
        assert (abs(val - val_ref) < tol)

    destroy_temp_case(tmpdir)


#def compare_reference(filename, test, key, vals, tol=1e-5):
#    folder = os.path.dirname(filename)
#    filemrk = filename.split("test_")[-1].split(".py")[0]
#    filename_ref = os.path.join(folder, "ref_{}.dat".format(filemrk))
#    print(filename_ref)
#    if not os.path.exists(filename_ref):
#        with open(filename_ref, "w"):
#            print("Creating reference file.")
#            pass
#    #struct = dict(test=dict(key=vals))
#    with open(filename_ref, "r") as infile:
#        ref = json.loads(infile.read())
#        print(ref)

if __name__ == "__main__":
    #if len(sys.argv) > 1 and sys.argv[1] == "generate":
    #    print("asdf")
    test_batchelor(0.1)
