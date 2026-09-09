"""A run stopped and resumed must be the run that was never stopped.

Deterministically, that is not a tolerance question: if the checkpoint holds
the whole state, the two agree bit for bit, and anything the checkpoint forgets
shows up as a difference. This is what catches state that was never written --
an edge's compressed time and rho_prev were absent until recently, so a
restarted diffusive-strip run silently resumed from an unmixed strip.

The final checkpoint is written one step past T, so the runs are stopped at
`n*dump_intv - dt` to put it on a round time that the uninterrupted run also
dumps at; otherwise the two never share a timestamp to compare.
"""

import os
import shutil
import subprocess

import h5py
import numpy as np
import pytest

from paths import REPO, app

PARTRAC = app("partrac")
SINE = os.path.join(REPO, "data_example", "sine_flow", "expr_params.dat")

DT, STOP, END = 0.0625, 1.0, 2.0
FIELDS = ("points", "edges", "dl", "dl0", "tau")

BASE = ("mode=analytic x0=0.5 y0=0.5 z0=0.5 Nrw=200 Nrw_max=200000 Dm=0 "
        "int_order=1 dt=%g stat_intv=1e9 random=false seed=1 dump_intv=%g "
        "integrate_tau=true tau_intv=%g tau_max=0" % (DT, STOP, DT)).split()

needs_partrac = pytest.mark.skipif(not os.path.exists(PARTRAC),
                                   reason="partrac is not built")


def call(case, extra):
    argv = [a for a in BASE if a.split("=")[0]
            not in {b.split("=")[0] for b in extra}] + extra
    r = subprocess.run([PARTRAC, str(case / "expr_params.dat")] + argv,
                       capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stdout + r.stderr


def state_at(case, t):
    """Every dumped field at time t, from whichever file holds it."""
    for f in sorted(case.rglob("data_from_t*.h5")):
        h = h5py.File(f, "r")
        for key in h:
            if abs(float(key) - t) < 1e-9:
                return {n: np.array(h[key + "/" + n]) for n in h[key]}
    raise AssertionError("no dump at t = %g under %s" % (t, case))


def continuous_and_resumed(tmp_path, extra):
    cont, split = tmp_path / "cont", tmp_path / "split"
    for d in (cont, split):
        d.mkdir(parents=True)
        shutil.copy(SINE, d / "expr_params.dat")

    call(cont, extra + ["T=%g" % (END + DT / 2), "checkpoint_intv=1e9"])
    call(split, extra + ["T=%g" % (STOP - DT), "checkpoint_intv=%g" % (STOP - DT)])
    checkpoint = list(split.rglob("edges.edge"))
    assert len(checkpoint) == 1
    call(split, extra + ["T=%g" % (END + DT / 2), "checkpoint_intv=1e9",
                         "restart_folder=" + str(checkpoint[0].parent.parent)])
    return state_at(cont, END), state_at(split, END)


@needs_partrac
def test_a_resumed_strip_is_identical_to_one_never_stopped(tmp_path):
    a, b = continuous_and_resumed(
        tmp_path, ["init_mode=strip_x", "La=0.5", "ds_max=0.01", "ds_min=0.002",
                   "refine=true", "refine_intv=%g" % DT,
                   "coarsen=true", "coarsen_intv=%g" % DT])
    assert len(a["edges"]) > 200        # it really refined past the initial 199
    assert set(FIELDS) <= set(a)
    for name in sorted(set(a) & set(b)):
        assert np.array_equal(a[name], b[name]), name


@needs_partrac
def test_a_resumed_point_cloud_is_identical_too(tmp_path):
    # no mesh, so a different path through the checkpoint
    a, b = continuous_and_resumed(
        tmp_path, ["init_mode=uniform_x", "La=0.5", "ds_max=1e9",
                   "ds_min=1e-12", "refine=false", "coarsen=false",
                   "integrate_tau=false"])
    for name in sorted(set(a) & set(b)):
        assert np.array_equal(a[name], b[name]), name
