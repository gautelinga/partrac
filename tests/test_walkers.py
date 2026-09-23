"""weighted_walkers against invariants that hold whatever the values are.

A walker advects with the explicit step, with noise if Dm > 0. At every
refine_intv, each walker past the exit plane is replaced by a copy of a
surviving walker j drawn with weight 2^-w_j; both copies take generation
w_j + 1, so a walker of generation w stands for 2^-w of an original, a split
conserves the weight sum of 2^-w, and only a walker that leaves loses weight.
At every stat_intv the walkers within ds_max of the source, along the strip
direction init_mode names when it lies in the exit plane, are written as
separation data.

The flow is plane Poiseuille, u = (0, 0, 1.5 (1 - x^2)): x never changes along
a path, so z(t) = z0 + u_z(x0) t exactly and every crossing of a plane z = Ln
can be predicted from the positions.
"""

import os

import numpy as np
import pytest

from dumps import all_dumps, read_stats
from paths import REPO, app
from runs import continuous_and_resumed, copy_example, run_app

WALKERS = app("weighted_walkers")
POISEUILLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")

needs_walkers = pytest.mark.skipif(not os.path.exists(WALKERS),
                                   reason="weighted_walkers is not built")

# a strip along x from wall to wall (x in [-1, 1], no width), in the plane
# Poiseuille flow u = (0, 0, 1.5 (1 - x^2)); nothing leaves unless an exit
# plane says so, and the walkers near the walls never reach one, so a
# resampling run always has survivors to copy from
BASE = ("init_mode=strip_x_y La=2.0 Lb=0.0 x0=0 y0=0 z0=0 Nrw=200 Nrw_max=2000 "
        "Dm=0 int_order=1 dt=0.01 T=0.5 dump_intv=0.05 stat_intv=0.05 "
        "checkpoint_intv=1e9 ds_max=0.1 Ln=1e9 Lt=0 exit_plane=none refine_intv=1e9 "
        "random=false seed=3").split()


def run(tmp_path, extra, name="case", env=None):
    """Run weighted_walkers on plane Poiseuille with extra overriding BASE; return the case folder."""
    run_app(WALKERS, copy_example(POISEUILLE, tmp_path / name), BASE, extra, env=env)
    return tmp_path / name


def u_z(x):
    """Return the plane Poiseuille axial velocity at x."""
    return 1.5 * (1.0 - x ** 2)      # u_inf = R = 1 in the example


# --- step -------------------------------------------------------------------

@needs_walkers
@pytest.mark.parametrize("int_order", [1, 2])
def test_advection_alone_is_exact(tmp_path, int_order):
    """x never changes and u_z depends only on x, so z(t) = z0 + u_z(x0) t
    exactly; the second-order term (a + J u) is zero here, so both orders hit
    it to round-off, and no walker gains a generation. Every other test here
    builds on this step."""
    d = run(tmp_path, "int_order=%d" % int_order)
    ds = all_dumps(d, raw=True)
    t0 = min(ds)
    x0 = ds[t0]["points"]
    assert len(ds) >= 6
    for t, g in ds.items():
        assert np.array_equal(g["points"][:, :2], x0[:, :2])
        assert np.allclose(g["points"][:, 2], x0[:, 2] + u_z(x0[:, 0]) * (t - t0),
                           rtol=0, atol=1e-12)
        assert np.all(g["w"] == 0)


@needs_walkers
def test_the_same_seed_gives_the_same_run(tmp_path):
    """With noise and resampling on, two runs with the same seed and thread
    count write bit-identical dumps, so a walker run can be reproduced."""
    extra = "Dm=1e-4 exit_plane=z Ln=0.3 refine_intv=0.01 dump_intv=0.1"
    out = []
    for name in ("a", "b"):
        d = run(tmp_path, extra, name, env={"OMP_NUM_THREADS": "4"})
        out.append(all_dumps(d, raw=True))
    assert out[0].keys() == out[1].keys()
    for t in out[0]:
        for k in out[0][t]:
            assert np.array_equal(out[0][t][k], out[1][t][k]), "%s at t = %g" % (k, t)


@needs_walkers
def test_resampling_does_not_depend_on_the_thread_count(tmp_path):
    """The resampling draws and the parents' generation updates are one serial
    pass on stream 0, so a resampled run without noise is bit-identical on 1
    and 4 threads. In a parallel pass two copies of a parent drawn twice could
    both read its old w, creating weight and tying the result to the thread
    count."""
    extra = "exit_plane=z Ln=0.3 refine_intv=0.01 dump_intv=0.05 Nrw=2000 Nrw_max=2000"
    out = {}
    for n in (1, 4):
        d = run(tmp_path, extra, str(n), env={"OMP_NUM_THREADS": str(n)})
        out[n] = all_dumps(d, raw=True)
    for t in out[1]:
        for k in out[1][t]:
            assert np.array_equal(out[1][t][k], out[4][t][k]), "%s at t = %g" % (k, t)


# --- resampling ---------------------------------------------------------------

@pytest.fixture(scope="module")
def resampled(tmp_path_factory):
    """The dumps of a run with an exit plane at z = 0.3, resampled and dumped every step."""
    if not os.path.exists(WALKERS):
        pytest.skip("weighted_walkers is not built")
    d = run(tmp_path_factory.mktemp("resampled"),
            "exit_plane=z Ln=0.3 refine_intv=0.01 dump_intv=0.01",
            env={"OMP_NUM_THREADS": "4"})
    return all_dumps(d, raw=True)


def test_the_exit_plane_is_enforced_and_weight_is_only_lost_by_leaving(resampled):
    """With the check at every step no dumped walker is past the plane, the
    count never changes and w is a non-negative integer. The weight sum of
    2^-w equals Nrw until the first walker leaves and never grows, since a
    split hands half the parent's weight to the copy and only a leaver's
    weight is lost; weight created or lost otherwise would bias every weighted
    statistic."""
    Ln = 0.3
    ds = resampled
    times = sorted(ds)
    N = len(ds[times[0]]["points"])
    total = []
    for t in times:
        g = ds[t]
        assert len(g["points"]) == N
        assert g["points"][:, 2].max() <= Ln + 1e-12
        w = g["w"][:, 0]
        assert np.array_equal(w, np.round(w)) and w.min() >= 0
        total.append(np.sum(2.0 ** -w))
    assert total[0] == N
    assert all(b <= a + 1e-12 for a, b in zip(total, total[1:]))
    # the fast walkers (x near 0, u_z = 1.5) cross z = 0.3 at t = 0.2, so by
    # the end there have been splits
    assert ds[times[-1]]["w"].max() >= 1
    assert total[-1] < N


def test_a_copy_starts_where_its_parent_is(resampled):
    """With Dm = 0 and a dump every step, a walker that would have advected
    past the plane reappears exactly on a surviving walker, inside the domain
    and at least one generation up. Only coincidence and w >= 1 are exact: a
    parent drawn twice in one pass has moved a further generation by the second
    copy, so the weight-sum test above is the conservation check."""
    Ln = 0.3
    ds = resampled
    times = sorted(ds)
    dt = times[1] - times[0]
    checked = 0
    for ta, tb in zip(times, times[1:]):
        a, b = ds[ta], ds[tb]
        # where each walker would be at tb, advected freely from ta
        z_free = a["points"][:, 2] + u_z(a["points"][:, 0]) * dt
        left = np.flatnonzero(z_free > Ln)
        stayed = np.flatnonzero(z_free <= Ln)
        for i in left:
            # replaced: it sits exactly on some surviving walker
            same = stayed[np.all(b["points"][stayed] == b["points"][i], axis=1)]
            assert len(same) >= 1, "walker %d at t = %g was not copied" % (i, tb)
            assert b["points"][i, 2] <= Ln
            assert b["w"][i] >= 1
            checked += 1
    assert checked > 0


@needs_walkers
def test_a_run_stops_when_every_walker_has_left(tmp_path):
    """A strip cut to |x| <= 0.5 has u_z >= 1.125, so it is swept out through
    z = 0.3 by t = 0.267. With no survivor there is no weight to draw a parent
    from, so the run must say so and stop there rather than copy leavers onto
    leavers, collapsing to one position with generations growing every step."""
    d = tmp_path / "swept"
    r = run_app(WALKERS, copy_example(POISEUILLE, d), BASE,
                "La=1.0 T=0.5 exit_plane=z Ln=0.3 refine_intv=0.01 dump_intv=0.01")
    assert "Every walker has crossed the exit plane" in r.stdout
    ds = all_dumps(d, raw=True)
    t_last = max(ds)
    assert 0.26 <= t_last <= 0.28, t_last                     # stopped there, not at T
    for t, g in ds.items():
        assert np.sum(2.0 ** -g["w"][:, 0]) >= 1.0            # never below one original's weight
        assert len(np.unique(g["points"], axis=0)) > 1


# --- separation data --------------------------------------------------------

@needs_walkers
@pytest.mark.parametrize("init_mode,axis,extra", [
    # a strip along x, which the exit plane z leaves as the axis
    ("strip_x_y", 0, "exit_plane=z Ln=0.3 refine_intv=0.01"),
    # a strip along y; nothing leaves, and the selection alone is checked
    ("strip_y_x", 1, "exit_plane=z Ln=1e9"),
])
def test_separation_data_is_the_dump_filtered_to_the_axis(tmp_path, init_mode, axis, extra):
    """Separation data is written at stat_intv, before the step, from the same
    state as the dump at that time: exactly the walkers within ds_max of the
    source along the strip direction init_mode names, for the exit plane set,
    in dump order, with their w. Otherwise the separation analysis would see
    different walkers or weights than the dump."""
    d = run(tmp_path, "init_mode=%s dump_intv=0.05 stat_intv=0.05 ds_max=0.2 %s" % (init_mode, extra),
               env={"OMP_NUM_THREADS": "4"})
    ds, sep = all_dumps(d, raw=True), all_dumps(d, "sepdata_from_t*.h5", raw=True)
    assert sep and set(sep) <= set(ds)
    picked = 0
    for t, s in sep.items():
        g = ds[t]
        pick = np.abs(g["points"][:, axis] - 0.0) < 0.2
        picked += pick.sum()
        assert np.array_equal(s["x"], g["points"][pick])
        assert np.array_equal(s["w"][:, 0] if s["w"].ndim == 2 else s["w"], g["w"][pick][:, 0])
    # the selection is a part of the strip, not none of it or all of it
    assert 0 < picked < sum(len(g["points"]) for t, g in ds.items() if t in sep)


# --- restart, circle init, statistics -----------------------------------------

@needs_walkers
def test_a_resumed_run_keeps_every_walkers_generation(tmp_path):
    """A checkpoint carries each walker's position and generation, so the
    resumed run's first dump equals the uninterrupted run's dump at that time,
    after the same splits. Generator state is not checkpointed, so only the
    restart time is compared; a lost generation would reset walkers' weights."""
    # the final checkpoint is written one step past T: t = 0.3, step 30
    cont, split = continuous_and_resumed(
        WALKERS, POISEUILLE, tmp_path,
        [BASE, "exit_plane=z Ln=0.3 refine_intv=0.01 dump_intv=0.05 checkpoint_intv=1e9"],
        "T=0.29", "T=0.4", env={"OMP_NUM_THREADS": "1"})
    a, b = all_dumps(cont, raw=True), all_dumps(split, raw=True)
    # the stopped run dumps up to 0.25, the resumed one from 0.3 on
    assert sorted(b) == sorted(a), sorted(b)
    assert a[0.3]["w"].max() >= 1, "no split before the restart: the test would not see a lost generation"
    for k in ("points", "w"):
        assert np.array_equal(a[0.3][k], b[0.3][k]), "%s differs at the restart" % k


@needs_walkers
def test_a_circle_spreads_only_in_the_directions_it_names(tmp_path):
    """circle_x_yz is a disc normal to x through x0, spread by Lb in y and z
    only, so every walker has x = x0 exactly; circle_x_y spreads in y only, so
    z stays on the disc. The third token of the mode names the spread
    directions, and a spread in another direction would start walkers off
    the intended source."""
    for mode, still in (("circle_x_yz", 0), ("circle_x_y", 2)):
        d = run(tmp_path, "init_mode=%s La=0.4 Lb=0.05 x0=0.1 y0=0.2 z0=0.3 T=0.01" % mode, mode)
        p = all_dumps(d, raw=True)[0.0]["points"]
        spread = [np.ptp(p[:, j]) for j in range(3)]
        if still == 0:
            assert np.all(p[:, 0] == 0.1), "x spread under %s" % mode
            assert spread[1] > 0.4 and spread[2] > 0.4
        else:
            assert np.all(np.abs(p[:, 2] - 0.3) <= 0.2 + 1e-12), "z spread under %s" % mode
            assert spread[1] > 0.4


@needs_walkers
def test_the_statistics_read_the_velocity(tmp_path):
    """The statistics evaluate the velocity whether or not it is dumped, so
    uz_mean at t = 0 is the mean of 1.5 (1 - x^2) over the walkers. Otherwise
    every velocity column of the statistics would be zero or stale."""
    d = run(tmp_path, "T=0.01")
    p = all_dumps(d, raw=True)[0.0]["points"]
    row = {k: v[0] for k, v in read_stats(d).items()}
    assert row["t"] == 0.0
    # the statistics file holds six significant digits
    assert np.isclose(row["uz_mean"], u_z(p[:, 0]).mean(), rtol=5e-6, atol=0), (row["uz_mean"], u_z(p[:, 0]).mean())
