"""The filaments app: its two edge resizes.

resize=doublings: an edge longer than the target is cut by 2^n, the least n
that brings it to the target or below; its reference length dl0 stays and n is
counted, so logelong = log(dl/dl0) + doublings log 2 is the total stretch.
resize=rescale brings the edge to the target and scales its reference with it,
so dl/dl0 carries the stretch instead. Pairs across plane Poiseuille,
u_z = 1.5 (1 - x^2), are sheared steadily, so every resize step has edges to
shorten.
"""

import os

import numpy as np
import pytest

from dumps import all_dumps, read_stats
from paths import REPO, app, needs
from runs import copy_example, run_app

FILAMENTS = app("filaments")
POISEUILLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")

# the resize target, and the initial edge length, so edges exceed it as soon as they stretch
TARGET = 0.05
BASE = ("mode=analytic init_mode=pairs_xz_x x0=0 y0=0 z0=0 Nrw=200 Nrw_max=400 int_order=1 "
        "ds_init=%g Dm=0 scheme=RK4 dt=0.01 dump_intv=0.1 stat_intv=0.1 checkpoint_intv=1e9 "
        "resize=doublings resize_target=ds_init resize_intv=0.01 random=false seed=1" % TARGET)

needs_filaments = needs(FILAMENTS)


def run(d, extra):
    """Run filaments on plane Poiseuille in d with `extra` overriding BASE; return d."""
    run_app(FILAMENTS, copy_example(POISEUILLE, d), BASE, extra)
    return d


@needs_filaments
def test_resizing_by_halvings_keeps_the_stretch(tmp_path):
    """With resize=doublings every dumped edge is at or below the target, dl0 is
    unchanged, the halving count never decreases, and logelong equals
    log(dl/dl0) + doublings log 2. The statistics report the mean of that
    logelong, so the stretch is not lost when edges are shortened."""
    d = run(tmp_path / "run", "T=0.5")
    g = all_dumps(d, raw=True)
    times = sorted(g)
    dl0 = g[times[0]]["dl0"][:, 0]
    prev = np.zeros(len(dl0))
    for t in times:
        e = g[t]
        dl, n, logelong = e["dl"][:, 0], e["doublings"][:, 0], e["logelong"][:, 0]
        assert np.array_equal(e["dl0"][:, 0], dl0)                 # the reference stays
        assert np.all(dl <= TARGET * (1 + 1e-12))                  # resized before the dump
        assert np.all(n >= prev)                                   # halvings only accumulate
        prev = n
        assert np.allclose(logelong, np.log(dl / dl0) + n * np.log(2), rtol=0, atol=1e-12)
    assert prev.max() >= 1                                         # something was halved
    st = read_stats(d)
    assert list(st)[-5:] == ["n_edges", "elong_mean", "elong2_mean", "logelong_mean", "logelong_var"]
    assert st["n_edges"][-1] == len(dl0)
    # the statistics are printed to six significant digits
    assert np.isclose(st["logelong_mean"][-1], logelong.mean(), rtol=5e-6, atol=1e-9)


@needs_filaments
def test_resizing_by_rescaling_scales_the_reference(tmp_path):
    """With resize=rescale edges stay at or below the target, no halving count or
    logelong is written, and the reference length shrinks with the edge so that
    dl/dl0 carries the stretch."""
    d = run(tmp_path / "run", "T=0.5 resize=rescale")
    g = all_dumps(d, raw=True)
    times = sorted(g)
    for t in times:
        e = g[t]
        assert "doublings" not in e and "logelong" not in e
        assert np.all(e["dl"][:, 0] <= TARGET * (1 + 1e-12))
    assert np.any(g[times[-1]]["dl0"][:, 0] < g[times[0]]["dl0"][:, 0])


@needs_filaments
@pytest.mark.parametrize("init_mode", ["points_xy", "uniform_x", "sheet_xy"])
def test_an_init_mode_other_than_pairs_is_refused(tmp_path, init_mode):
    """filaments follows pairs of points and their stretching, so it starts
    from pair_* or pairs_* only; any other init_mode is refused before the run
    starts rather than run as edges the app was not written for."""
    r = run_app(FILAMENTS, copy_example(POISEUILLE, tmp_path / "c"), BASE,
                "init_mode=" + init_mode, "T=0", check=False)
    assert r.returncode == 2, r.stdout[-1000:] + r.stderr[-1000:]
    assert "init_mode must be pair_* or pairs_*" in r.stderr, r.stderr
