"""Trajectory-level checks on the integrators and the restart path.

No golden numbers: these assert invariants the code must satisfy whatever the
values are. A run must be reproducible, and stopping at T/2 and restarting must
land in the same place as running straight through.
"""

import os
import shutil
import subprocess

import numpy as np
import pytest

from paths import REPO, app

PARTRAC = app("partrac")
RK4APP = app("tracervectors_analyticRK4")
EXAMPLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")

# Dm = 0 so there is no Brownian term, and random=false to pin the seed
BASE = ("mode=analytic init_mode=uniform_x Nrw=100 Nrw_max=5000 ds_max=0.4 "
        "ds_min=0.1 Dm=0 int_order=1 dt=0.01 dump_intv=0.05 stat_intv=0.05 "
        "checkpoint_intv=0.05 random=false seed=3").split()

# checkpoints are written at max_digits10, so a restart round-trips exactly


def case(tmp_path, name):
    d = tmp_path / name
    d.mkdir()
    shutil.copy(EXAMPLE, d / "expr_params.dat")
    return d


def run(binary, case_dir, args):
    r = subprocess.run([binary, str(case_dir / "expr_params.dat")] + args,
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr
    return r


def checkpoints(case_dir):
    """Checkpoint position files, outermost run first.

    A restarted run appends the rank to the folder a second time, so its
    checkpoint lands one level deeper than the run it resumed from.
    """
    found = sorted(case_dir.glob("**/Checkpoints/positions.pos"),
                   key=lambda p: len(p.parts))
    assert found, "no checkpoint written"
    return found


def final_positions(path):
    return np.loadtxt(path)


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
@pytest.mark.parametrize("scheme", ["explicit", "RK4"])
def test_run_is_reproducible(tmp_path, scheme):
    a, b = case(tmp_path, "a"), case(tmp_path, "b")
    args = BASE + ["T=0.2", "scheme=" + scheme]
    run(PARTRAC, a, args)
    run(PARTRAC, b, args)
    xa = final_positions(checkpoints(a)[0])
    xb = final_positions(checkpoints(b)[0])
    assert xa.shape == xb.shape
    assert np.array_equal(xa, xb)


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
@pytest.mark.parametrize("scheme", ["explicit", "RK4"])
def test_restart_matches_an_uninterrupted_run(tmp_path, scheme):
    direct, split = case(tmp_path, "direct"), case(tmp_path, "split")
    args = BASE + ["scheme=" + scheme]
    run(PARTRAC, direct, args + ["T=0.4"])

    run(PARTRAC, split, args + ["T=0.2"])
    resume_from = checkpoints(split)[0].parent.parent
    run(PARTRAC, split, args + ["T=0.4", "restart_folder=" + str(resume_from)])

    xa = final_positions(checkpoints(direct)[0])
    xb = final_positions(checkpoints(split)[-1])
    assert xa.shape == xb.shape
    assert np.array_equal(xa, xb)


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_command_line_beats_the_checkpoint_on_restart(tmp_path):
    d = case(tmp_path, "c")
    run(PARTRAC, d, BASE + ["T=0.2"])
    resume_from = checkpoints(d)[0].parent.parent
    # the checkpoint holds T=0.2; the command line must win
    r = run(PARTRAC, d, BASE + ["T=0.4", "restart_folder=" + str(resume_from)])
    times = [float(l.split("=")[1]) for l in r.stdout.splitlines()
             if l.startswith("Time =")]
    assert times, r.stdout
    assert max(times) > 0.2


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
@pytest.mark.parametrize("scheme,int_order", [("explicit", 1), ("explicit", 2),
                                              ("RK4", 1)])
def test_plane_poiseuille_advects_only_along_the_axis(tmp_path, scheme, int_order):
    # u is axial and depends only on x, so a particle must keep its x and y
    # whatever the scheme. Catches an integrator that mixes up components.
    d = case(tmp_path, "axis")
    args = [a for a in BASE if not a.startswith("int_order=")]
    run(PARTRAC, d, args + ["T=0.3", "scheme=" + scheme,
                            "int_order=%d" % int_order])
    h5py = pytest.importorskip("h5py")
    dumps = sorted(d.glob("**/data_from_t*.h5"))
    assert dumps
    with h5py.File(dumps[0]) as h:
        keys = sorted(h.keys(), key=float)
        first = np.array(h[keys[0]]["points"])
        last = np.array(h[keys[-1]]["points"])
    assert np.abs(last[:, 0] - first[:, 0]).max() < 1e-12
    assert np.abs(last[:, 1] - first[:, 1]).max() < 1e-12
    assert np.abs(last[:, 2] - first[:, 2]).max() > 1e-3  # it did move


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_statistics_do_not_depend_on_dump_verbosity(tmp_path):
    # minimal_output selects what is written to the h5 dump; it must not change
    # the statistics, which read the velocity whether or not it is dumped
    means = {}
    for minimal in ("true", "false"):
        d = case(tmp_path, "mo_" + minimal)
        run(PARTRAC, d, BASE + ["T=0.05", "minimal_output=" + minimal])
        stats = list(d.glob("**/tdata_from_t*.dat"))
        assert len(stats) == 1
        row = [l for l in stats[0].read_text().splitlines()
               if l and not l.startswith("#")][0].split()
        means[minimal] = [float(v) for v in row[7:10]]  # ux_mean, uy_mean, uz_mean
        assert all(np.isfinite(v) for v in means[minimal]), means[minimal]

    assert means["true"] == means["false"]
    assert abs(means["false"][2]) > 1e-6  # the flow is axial, so uz is nonzero

    # and it still trims the dump
    sizes = {}
    h5py = pytest.importorskip("h5py")
    for minimal in ("true", "false"):
        dump = sorted((tmp_path / ("mo_" + minimal)).glob("**/data_from_t*.h5"))[0]
        with h5py.File(dump) as f:
            sizes[minimal] = len(f[sorted(f.keys())[0]].keys())
    assert sizes["true"] < sizes["false"]


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_diffusive_run_is_reproducible(tmp_path):
    # each thread draws from its own generator, so a run is only reproducible
    # while the schedule is static and particle i always lands on the same
    # thread. Switching the integrator to schedule(dynamic) breaks this, and
    # the Dm=0 tests above would not notice.
    args = [a for a in BASE if not a.startswith("Dm=")]
    args += ["Dm=1e-4", "T=0.1", "num_threads=4"]
    a, b = case(tmp_path, "d1"), case(tmp_path, "d2")
    run(PARTRAC, a, args)
    run(PARTRAC, b, args)
    xa = final_positions(checkpoints(a)[0])
    xb = final_positions(checkpoints(b)[0])
    assert xa.shape == xb.shape
    assert np.array_equal(xa, xb)


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
@pytest.mark.parametrize("threads", [4, 16])
def test_results_do_not_depend_on_the_thread_count(tmp_path, threads):
    # every parallel loop writes one entry per particle, so the result must be
    # identical whatever the schedule; the statistics reductions stay serial
    ref = case(tmp_path, "t1")
    other = case(tmp_path, "t%d" % threads)
    args = [a for a in BASE if not a.startswith(("int_order=", "stat_intv="))]
    args += ["T=0.2", "int_order=2", "stat_intv=0.01"]
    run(PARTRAC, ref, args + ["num_threads=1"])
    run(PARTRAC, other, args + ["num_threads=%d" % threads])

    assert np.array_equal(final_positions(checkpoints(ref)[0]),
                          final_positions(checkpoints(other)[0]))
    sa = list(ref.glob("**/tdata_from_t*.dat"))[0].read_text()
    sb = list(other.glob("**/tdata_from_t*.dat"))[0].read_text()
    assert sa == sb


WALKERS = app("weighted_walkers")
BRINKMAN = os.path.join(REPO, "data_example", "brinkman_cylinder", "expr_params.dat")


@pytest.mark.skipif(not os.path.exists(WALKERS), reason="weighted_walkers is not built")
def test_separation_data_covers_every_exit_plane(tmp_path):
    # for dim == 2 the selection is gated on exit_plane and checks the distance
    # along the strip direction named by init_mode; exit_plane=z had no branch
    # at all and silently wrote empty separation data
    d = tmp_path / "ww"
    d.mkdir()
    shutil.copy(BRINKMAN, d / "expr_params.dat")
    args = ("Dm=1e-5 stat_intv=0.05 dt=1e-2 Nrw_max=1e4 Nrw=200 int_order=2 "
            "init_mode=strip_y_x La=0.5 Lb=0.0 T=0.1 x0=-1 y0=0 z0=0 Ln=1e9 "
            "ds_max=2.0 Lt=0 num_threads=1 dump_intv=0.05 refine_intv=0.05 "
            "seed=5").split()
    run(WALKERS, d, args + ["exit_plane=z"])

    h5py = pytest.importorskip("h5py")
    sep = sorted(d.glob("**/sepdata_from_t*.h5"))
    assert sep
    with h5py.File(sep[0]) as f:
        selected = len(f[sorted(f.keys())[0]]["w"])
    # the strip runs along y, which exit_plane=z does not exclude
    assert selected > 0


@pytest.mark.skipif(not os.path.exists(RK4APP), reason="RK4 app is not built")
def test_rk4_app_is_reproducible(tmp_path):
    # the experimental integrator, which the apps use instead of RKIntegrator
    a, b = case(tmp_path, "a"), case(tmp_path, "b")
    args = ("init_mode=points_xy Nrw=100 Nrw_max=5000 Dm=0 dt=0.005 T=0.1 "
            "dump_intv=0.05 stat_intv=0.05 random=false seed=7").split()
    run(RK4APP, a, args)
    run(RK4APP, b, args)
    fa = sorted(a.glob("**/data_from_t*.h5"))
    fb = sorted(b.glob("**/data_from_t*.h5"))
    assert fa and len(fa) == len(fb)
    h5py = pytest.importorskip("h5py")
    with h5py.File(fa[0]) as ga, h5py.File(fb[0]) as gb:
        assert sorted(ga.keys()) == sorted(gb.keys())
        for k in ga:
            for d in ga[k]:
                assert np.array_equal(np.array(ga[k][d]), np.array(gb[k][d]))
