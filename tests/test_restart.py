"""A run stopped and resumed from a checkpoint is the run that was never stopped.

partrac is deterministic here (random=false, Dm=0, sine flow), so this is not a
tolerance question: if the checkpoint holds the whole state, the continuous and
resumed runs agree bit for bit at t = END, and any state the checkpoint omits
(node positions, edges, dl0, an edge's compressed time tau and rho_prev, the
step count that sets the phase of the remeshing intervals) shows up as a
difference.

The final checkpoint is written one step past T, so the split runs stop at
`n*dump_intv - dt`; that puts the checkpoint on a round time the uninterrupted
run also dumps at, and the two runs share timestamps to compare.
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

# DT divides STOP and END exactly in binary, so step counts are exact
DT, STOP, END = 0.0625, 1.0, 2.0
FIELDS = ("points", "edges", "dl", "dl0", "tau")

BASE = ("mode=analytic x0=0.5 y0=0.5 z0=0.5 Nrw=200 Nrw_max=200000 Dm=0 "
        "int_order=1 dt=%g stat_intv=1e9 random=false seed=1 dump_intv=%g "
        "integrate_tau=true tau_intv=%g tau_max=0" % (DT, STOP, DT)).split()

needs_partrac = pytest.mark.skipif(not os.path.exists(PARTRAC),
                                   reason="partrac is not built")


def call(case, extra):
    """Run partrac in `case` with BASE, overridden by `extra`; assert success."""
    argv = [a for a in BASE if a.split("=")[0]
            not in {b.split("=")[0] for b in extra}] + extra
    r = subprocess.run([PARTRAC, str(case / "expr_params.dat")] + argv,
                       capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stdout + r.stderr


def state_at(case, t):
    """Every field dumped at time t, from whichever dump file holds it."""
    for f in sorted(case.rglob("data_from_t*.h5")):
        h = h5py.File(f, "r")
        for key in h:
            if abs(float(key) - t) < 1e-9:
                return {n: np.array(h[key + "/" + n]) for n in h[key]}
    raise AssertionError("no dump at t = %g under %s" % (t, case))


def continuous_and_resumed(tmp_path, extra, stop=STOP - DT):
    """State at END of an uninterrupted run and of one checkpointed at stop and resumed."""
    cont, split = tmp_path / "cont", tmp_path / "split"
    for d in (cont, split):
        d.mkdir(parents=True)
        shutil.copy(SINE, d / "expr_params.dat")

    # T is half a step past END so the last step lands on END despite rounding
    call(cont, extra + ["T=%g" % (END + DT / 2), "checkpoint_intv=1e9"])
    call(split, extra + ["T=%g" % stop, "checkpoint_intv=%g" % stop])
    checkpoint = list(split.rglob("edges.edge"))
    assert len(checkpoint) == 1
    call(split, extra + ["T=%g" % (END + DT / 2), "checkpoint_intv=1e9",
                         "restart_folder=" + str(checkpoint[0].parent.parent)])
    return state_at(cont, END), state_at(split, END)


@needs_partrac
def test_a_resumed_strip_is_identical_to_one_never_stopped(tmp_path):
    """A refining and coarsening strip with tau integration resumes bit for bit.
    If any edge state (dl0, tau, rho_prev) were missing from the checkpoint, a
    restarted diffusive-strip run would silently continue from a wrong, e.g.
    unmixed, strip."""
    a, b = continuous_and_resumed(
        tmp_path, ["init_mode=strip_x", "La=0.5", "ds_max=0.01", "ds_min=0.002",
                   "refine=true", "refine_intv=%g" % DT,
                   "coarsen=true", "coarsen_intv=%g" % DT])
    assert len(a["edges"]) > 200        # refined past the initial 199 edges
    assert set(FIELDS) <= set(a)
    for name in sorted(set(a) & set(b)):
        assert np.array_equal(a[name], b[name]), name


@needs_partrac
def test_a_resumed_point_cloud_is_identical_too(tmp_path):
    """A run with no mesh resumes bit for bit as well. It goes through a
    different path of the checkpoint reader than a strip does."""
    a, b = continuous_and_resumed(
        tmp_path, ["init_mode=uniform_x", "La=0.5", "ds_max=1e9",
                   "ds_min=1e-12", "refine=false", "coarsen=false",
                   "integrate_tau=false"])
    for name in sorted(set(a) & set(b)):
        assert np.array_equal(a[name], b[name]), name


@needs_partrac
def test_a_resume_that_is_not_on_a_remeshing_step_is_identical_too(tmp_path):
    """A resume between remeshing steps keeps the remeshing phase, because the
    step count is restored from the checkpoint. Otherwise the resumed run would
    remesh on different steps and diverge from the uninterrupted one."""
    # remeshing every 4 dt; the checkpoint at t = 0.875 is not a multiple of
    # 4 dt = 0.25, so it falls between remeshing steps
    a, b = continuous_and_resumed(
        tmp_path, ["init_mode=strip_x", "La=0.5", "ds_max=0.01", "ds_min=0.002",
                   "refine=true", "refine_intv=%g" % (4 * DT),
                   "coarsen=true", "coarsen_intv=%g" % (4 * DT),
                   "dump_intv=%g" % DT],
        stop=0.875 - DT)
    assert len(a["edges"]) > 200        # refined past the initial 199 edges
    for name in sorted(set(a) & set(b)):
        assert np.array_equal(a[name], b[name]), name


@needs_partrac
def test_a_resumed_run_told_not_to_remesh_does_not(tmp_path):
    """With refine_intv = coarsen_intv = 0 remeshing is off; the initial pass runs
    once before the checkpoint and is not repeated on resume. Rerunning it would
    change the mesh of a run the user asked not to remesh."""
    a, b = continuous_and_resumed(
        tmp_path, ["init_mode=strip_x", "La=0.5", "ds_max=0.01", "ds_min=0.002",
                   "refine=true", "refine_intv=0",
                   "coarsen=true", "coarsen_intv=0"])
    assert len(a["edges"]) == len(b["edges"])
    for name in sorted(set(a) & set(b)):
        assert np.array_equal(a[name], b[name]), name


# --- filaments ------------------------------------------------------------------

FILAMENTS = app("filaments")
ABC = os.path.join(REPO, "data_example", "abc_flow_unsteady", "expr_params.dat")


@pytest.mark.skipif(not os.path.exists(FILAMENTS), reason="filaments is not built")
def test_a_resumed_filament_run_keeps_the_phase_of_its_intervals(tmp_path):
    """filaments resumes its step count from the checkpoint, so intervals keep
    their phase and the resumed run is bitwise identical to a continuous one.
    Here the resize runs every 3 steps and the run stops at step 20, which 3 does
    not divide, so a step count reset to zero would resize on the wrong steps."""
    from dumps import dump_at
    pi = "3.14159265358979"
    base = ("mode=analytic init_mode=pairs_xyz Nrw=50 Nrw_max=500 int_order=1 "
            "ds_max=0.1 ds_min=0.01 ds_init=0.1 x0=%s y0=%s z0=%s Dm=0 scheme=RK4 "
            "dt=0.01 dump_intv=0.1 stat_intv=1e9 checkpoint_intv=1e9 resize_intv=0.03 "
            "random=false seed=1" % (pi, pi, pi)).split()

    def go(d, extra):
        d.mkdir(parents=True, exist_ok=True)
        shutil.copy(ABC, d / "expr_params.dat")
        keys = {a.split("=")[0] for a in extra}
        argv = [a for a in base if a.split("=")[0] not in keys] + extra
        r = subprocess.run([FILAMENTS, str(d / "expr_params.dat")] + argv,
                           capture_output=True, text=True, timeout=600)
        assert r.returncode == 0, r.stdout + r.stderr

    cont, split = tmp_path / "cont", tmp_path / "split"
    go(cont, ["T=0.4"])
    # the final checkpoint is written one step past T: t = 0.2, step 20
    go(split, ["T=0.19"])
    folder = os.path.dirname(os.path.dirname(next(split.rglob("Checkpoints/positions.pos"))))
    go(split, ["T=0.4", "restart_folder=" + folder])
    for t in (0.3, 0.4):
        a, b = dump_at(cont, t), dump_at(split, t)
        assert set(a) == set(b)
        for k in a:
            assert np.array_equal(a[k], b[k]), "%s differs at t = %g after a restart" % (k, t)
