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

import numpy as np
import pytest

from cases import ABC
from dumps import all_dumps, dump_at
from paths import REPO, app
from runs import checkpoint_folder, continuous_and_resumed, copy_example, run_app

PARTRAC = app("partrac")
TRACERS = app("tracers")
SINE = os.path.join(REPO, "data_example", "sine_flow", "expr_params.dat")

# DT divides STOP and END exactly in binary, so step counts are exact
DT, STOP, END = 0.0625, 1.0, 2.0
FIELDS = ("points", "edges", "dl", "dl0", "tau")

BASE = ("mode=analytic x0=0.5 y0=0.5 z0=0.5 Nrw=200 Nrw_max=200000 Dm=0 "
        "int_order=1 dt=%g stat_intv=1e9 random=false seed=1 dump_intv=%g "
        "integrate_tau=true tau_intv=%g tau_max=0" % (DT, STOP, DT)).split()

needs_partrac = pytest.mark.skipif(not os.path.exists(PARTRAC),
                                   reason="partrac is not built")


def continuous_and_resumed_at_end(tmp_path, extra, stop=STOP - DT, text=False):
    """Every dataset at END, as written, of an uninterrupted run and of one checkpointed at stop and resumed.

    With text the checkpoint is rewritten in the old text format before the resume.
    """
    # T is half a step past END so the last step lands on END despite rounding
    cont, split = continuous_and_resumed(
        PARTRAC, SINE, tmp_path, [BASE, extra],
        ["T=%g" % stop, "checkpoint_intv=%g" % stop],
        ["T=%g" % (END + DT / 2), "checkpoint_intv=1e9"], text=text)
    a, b = dump_at(cont, END, raw=True), dump_at(split, END, raw=True)
    assert set(a) == set(b)
    return a, b


@needs_partrac
@pytest.mark.parametrize("fmt", ["hdf5", "text"])
def test_a_resumed_strip_is_identical_to_one_never_stopped(tmp_path, fmt):
    """A refining and coarsening strip with tau integration resumes bit for bit.
    If any edge state (dl0, tau, rho_prev) were missing from the checkpoint, a
    restarted diffusive-strip run would silently continue from a wrong, e.g.
    unmixed, strip. A text checkpoint, as older runs wrote them, resumes the
    same way."""
    a, b = continuous_and_resumed_at_end(
        tmp_path, ["init_mode=strip_x", "La=0.5", "ds_max=0.01", "ds_min=0.002",
                   "refine=true", "refine_intv=%g" % DT,
                   "coarsen=true", "coarsen_intv=%g" % DT], text=fmt == "text")
    assert len(a["edges"]) > 200        # refined past the initial 199 edges
    assert set(FIELDS) <= set(a)
    for name in sorted(a):
        assert np.array_equal(a[name], b[name]), name


@needs_partrac
@pytest.mark.parametrize("scheme", ["explicit", "RK4"])
def test_a_resumed_point_cloud_is_identical_too(tmp_path, scheme):
    """A run with no mesh resumes bit for bit as well, with either scheme. It
    goes through a different path of the checkpoint reader than a strip does."""
    a, b = continuous_and_resumed_at_end(
        tmp_path, ["init_mode=uniform_x", "La=0.5", "ds_max=1e9",
                   "ds_min=1e-12", "refine=false", "coarsen=false",
                   "integrate_tau=false", "scheme=" + scheme])
    for name in sorted(a):
        assert np.array_equal(a[name], b[name]), name


@needs_partrac
def test_a_resume_that_is_not_on_a_remeshing_step_is_identical_too(tmp_path):
    """A resume between remeshing steps keeps the remeshing phase, because the
    step count is restored from the checkpoint. Otherwise the resumed run would
    remesh on different steps and diverge from the uninterrupted one."""
    # remeshing every 4 dt; the checkpoint at t = 0.875 is not a multiple of
    # 4 dt = 0.25, so it falls between remeshing steps
    a, b = continuous_and_resumed_at_end(
        tmp_path, ["init_mode=strip_x", "La=0.5", "ds_max=0.01", "ds_min=0.002",
                   "refine=true", "refine_intv=%g" % (4 * DT),
                   "coarsen=true", "coarsen_intv=%g" % (4 * DT),
                   "dump_intv=%g" % DT],
        stop=0.875 - DT)
    assert len(a["edges"]) > 200        # refined past the initial 199 edges
    for name in sorted(a):
        assert np.array_equal(a[name], b[name]), name


@needs_partrac
def test_a_resumed_run_told_not_to_remesh_does_not(tmp_path):
    """With refine_intv = coarsen_intv = 0 remeshing is off; the initial pass runs
    once before the checkpoint and is not repeated on resume. Rerunning it would
    change the mesh of a run the user asked not to remesh."""
    a, b = continuous_and_resumed_at_end(
        tmp_path, ["init_mode=strip_x", "La=0.5", "ds_max=0.01", "ds_min=0.002",
                   "refine=true", "refine_intv=0",
                   "coarsen=true", "coarsen_intv=0"])
    assert len(a["edges"]) == len(b["edges"])
    for name in sorted(a):
        assert np.array_equal(a[name], b[name]), name


# --- filaments ------------------------------------------------------------------

FILAMENTS = app("filaments")


@pytest.mark.skipif(not os.path.exists(FILAMENTS), reason="filaments is not built")
@pytest.mark.parametrize("resize,fmt", [("rescale", "hdf5"), ("doublings", "hdf5"), ("doublings", "text")])
def test_a_resumed_filament_run_keeps_the_phase_of_its_intervals(tmp_path, resize, fmt):
    """filaments resumes its step count from the checkpoint, so intervals keep
    their phase and the resumed run is bitwise identical to a continuous one.
    Here the resize runs every 3 steps and the run stops at step 20, which 3 does
    not divide, so a step count reset to zero would resize on the wrong steps.
    With resize=doublings the checkpoint also carries each edge's halvings;
    without them a resumed run's logelong would lose the stretch halved away
    before the stop. An old text checkpoint carries them too."""
    pi = "3.14159265358979"
    base = ("mode=analytic init_mode=pairs_xyz Nrw=50 Nrw_max=500 int_order=1 "
            "ds_max=0.1 ds_min=0.01 ds_init=0.1 x0=%s y0=%s z0=%s Dm=0 scheme=RK4 "
            "dt=0.01 dump_intv=0.1 stat_intv=1e9 checkpoint_intv=1e9 resize_intv=0.03 "
            "random=false seed=1 resize=%s" % (pi, pi, pi, resize))
    # the final checkpoint is written one step past T: t = 0.2, step 20
    cont, split = continuous_and_resumed(FILAMENTS, ABC, tmp_path, base, "T=0.19", "T=0.4",
                                         text=fmt == "text")
    for t in (0.3, 0.4):
        a, b = dump_at(cont, t), dump_at(split, t)
        assert set(a) == set(b)
        for k in a:
            assert np.array_equal(a[k], b[k]), "%s differs at t = %g after a restart" % (k, t)
    if resize == "doublings":
        assert dump_at(cont, 0.2)["doublings"].max() >= 1, "nothing was halved before the stop"


@pytest.mark.skipif(not os.path.exists(TRACERS), reason="tracers is not built")
def test_a_run_resumed_at_another_dt_goes_on_from_the_checkpoints_time(tmp_path):
    """Time is counted in steps from the run's start, or from the restart if dt
    has changed there: steps of 0.1 checkpointed at t = 0.6 and resumed at
    dt = 0.05 dump at 0.6, 0.7, ... up to T. Counted from the start, the
    resumed steps would be timed from 0.3."""
    params = copy_example(ABC, tmp_path)
    base = ("mode=analytic init_mode=points_xyz x0=3 y0=3 z0=3 Nrw=20 Nrw_max=20 Dm=0 int_order=1 "
            "scheme=RK4 stat_intv=1e9 dump_intv=0.1 random=false seed=1")
    # the final checkpoint is written one step past T: t = 0.6, step 6
    run_app(TRACERS, params, base, "dt=0.1 T=0.5 checkpoint_intv=0.5")
    run_app(TRACERS, params, base, "dt=0.05 T=1.0 checkpoint_intv=1e9",
            "restart_folder=%s" % checkpoint_folder(tmp_path))
    # the resumed run's own file, named by the time it starts at
    resumed = sorted(all_dumps(tmp_path, "data_from_t0.6*.h5"))
    assert len(resumed) == 5 and np.allclose(resumed, [0.6, 0.7, 0.8, 0.9, 1.0], rtol=0, atol=1e-9), resumed


# --- dump files ------------------------------------------------------------------

POISEUILLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")
POINTS = ("mode=analytic init_mode=points_xyz x0=0 y0=0 z0=0 Nrw=20 Nrw_max=20 Dm=0 int_order=1 "
          "scheme=RK4 dt=0.1 stat_intv=1e9 random=false seed=1")
needs_tracers = pytest.mark.skipif(not os.path.exists(TRACERS), reason="tracers is not built")


@needs_tracers
def test_a_dump_file_holds_dump_chunk_size_dumps(tmp_path):
    """dump_chunk_size counts dumps: with a dump every two steps (dump_intv =
    0.25 floored to 0.2) and three dumps a file, the files start at 0, 0.6 and
    1.2 and hold 0, 0.2, 0.4; 0.6, 0.8, 1.0; and 1.2."""
    h5py = pytest.importorskip("h5py")
    params = copy_example(POISEUILLE, tmp_path)
    run_app(TRACERS, params, POINTS, "T=1.2 dump_intv=0.25 dump_chunk_size=3")
    files = sorted(tmp_path.rglob("data_from_t*.h5"))
    assert [f.name for f in files] == ["data_from_t%f.h5" % t for t in (0., 0.6, 1.2)]
    sizes = []
    for f in files:
        with h5py.File(f, "r") as h:
            sizes.append(len(h))
    assert sizes == [3, 3, 1]


@needs_tracers
def test_chunked_dumps_survive_a_resume(tmp_path):
    """A run with dump files of two dumps each, stopped and resumed: the
    resumed run starts its own file and the next chunk after it, and every
    time is dumped once."""
    args = [POINTS, "dump_intv=0.1 dump_chunk_size=2"]
    params = copy_example(POISEUILLE, tmp_path)
    run_app(TRACERS, params, args, "T=0.4")
    run_app(TRACERS, params, args, "T=0.7", "restart_folder=%s" % checkpoint_folder(tmp_path))
    files = sorted(f.name for f in tmp_path.rglob("data_from_t*.h5"))
    assert files == ["data_from_t%f.h5" % t for t in (0., 0.2, 0.4, 0.5, 0.6)], files
    assert sorted(all_dumps(tmp_path)) == pytest.approx([0.1 * k for k in range(8)])

