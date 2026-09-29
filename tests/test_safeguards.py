"""What a run reports about its own trustworthiness.

- n_at_rest, the last statistics column of the tracer apps: tracers slower
  than 1e-3 of the mean speed. A field that traps tracers at walls or
  stagnation points shows nothing else: the mean velocity barely moves while
  the slow layer empties into the traps.
- logdetF_mean and logdetF_var (tracertensors): log det F, zero in an
  incompressible flow, so a drift is the field's divergence, which the
  stretching statistics would otherwise absorb.
- the step-size note: a step that crosses more than one cell cuts across
  streamlines and wall layers, which reads as trapping or false stretching.
- random = true draws one seed, prints it and names the run folder after it,
  so the run can be repeated with random = false and that seed; before, each
  thread drew its own and the folder named a seed that was never used.
"""

import os
import re

import numpy as np
import pytest

from cases import args_for
from dumps import dump_at, read_stats
from paths import REPO, app
from runs import checkpoint_folder, copy_case, copy_example, run_app

TRACERS = app("tracers")
TENSORS = app("tracertensors")
LINEAR = os.path.join(REPO, "data_example", "linear_flow", "expr_params.dat")

needs = pytest.mark.skipif(not (os.path.exists(TRACERS) and os.path.exists(TENSORS)),
                           reason="tracers or tracertensors is not built")

BASE = ("mode=analytic x0=0 y0=0 z0=0 Nrw=50 Nrw_max=50 Dm=0 int_order=1 "
        "checkpoint_intv=1e9 random=false seed=1").split()


def linear_flow(d, entries):
    """A copy of the linear-flow example in d with its gradient entries replaced by `entries`."""
    params = copy_example(LINEAR, d)
    lines = [l for l in params.read_text().splitlines() if not re.match(r"A[xyz]{2}=", l)]
    params.write_text("\n".join(lines + ["%s=%r" % kv for kv in entries.items()]) + "\n")
    return params


@needs
def test_every_tracer_on_a_stagnation_line_is_at_rest(tmp_path):
    """u = (x, 0, 0) vanishes on the y axis: 50 tracers placed there are all at rest."""
    params = linear_flow(tmp_path / "c", {"Axx": 1.0})
    run_app(TRACERS, params, BASE, "init_mode=points_y dt=0.01 T=0.1 stat_intv=0.05 dump_intv=1e9")
    st = read_stats(tmp_path / "c")
    assert list(st)[-1] == "n_at_rest"
    assert (st["n_at_rest"] == 50).all(), st["n_at_rest"]


@needs
def test_the_count_is_the_tracers_slower_than_a_thousandth_of_the_mean_speed(tmp_path):
    """At a stagnation point the count at each statistics time is the dumped
    velocities' own: |u| <= 1e-3 of their mean |u|, and not every tracer."""
    params = linear_flow(tmp_path / "c", {"Axx": 1.0, "Ayy": -1.0})
    run_app(TRACERS, params, BASE, "init_mode=points_x dt=0.01 T=0.5 stat_intv=0.25 dump_intv=0.25")
    st = read_stats(tmp_path / "c")
    for t, n in zip(st["t"], st["n_at_rest"]):
        speed = np.linalg.norm(dump_at(tmp_path / "c", t)["u"], axis=1)
        assert n == np.sum(speed <= 1e-3 * speed.mean()), (t, n)
        assert n < 50


@needs
@pytest.mark.parametrize("axx,ayy", [(1.0, -1.0), (1.0, 0.5)])
def test_log_det_F_is_the_integral_of_the_divergence(tmp_path, axx, ayy):
    """In u = A x, F = exp(A t) and log det F = tr(A) t for every particle:
    zero at the incompressible stagnation point, 1.5 t where the flow expands."""
    params = linear_flow(tmp_path / "c", {"Axx": axx, "Ayy": ayy})
    T = 0.5
    run_app(TENSORS, params, BASE, "init_mode=points_x dt=0.01 T=%g stat_intv=%g dump_intv=1e9" % (T, T))
    st = {k: v[-1] for k, v in read_stats(tmp_path / "c").items()}
    assert st["t"] == pytest.approx(T)
    assert st["logdetF_mean"] == pytest.approx((axx + ayy) * T, abs=1e-8)
    assert st["logdetF_var"] < 1e-16
    assert list(st)[-1] == "n_at_rest"


@pytest.mark.skipif(not os.path.exists(TRACERS), reason="tracers is not built")
def test_a_step_over_several_cells_is_noted(tmp_path, felbm_dir):
    """The felbm shear, |u| <= 0.05 on a unit grid: dt = 100 crosses up to 5
    cells a step and is noted; dt = 1 crosses a twentieth and is not."""
    args = args_for("tracers", "felbm")
    for dt, noted in ((100.0, True), (1.0, False)):
        d = copy_case(felbm_dir, tmp_path / ("dt%g" % dt))
        r = run_app(TRACERS, "./" + os.path.relpath(d / "felbm_params.dat"), args,
                    "dt=%g T=%g stat_intv=%g dump_intv=1e9" % (dt, 2 * dt if dt < 100 else dt, dt))
        m = re.search(r"Note: a step crosses up to ([0-9.e+-]+) cells", r.stdout)
        assert bool(m) == noted, r.stdout[-2000:]
        if m:
            assert 4.0 < float(m.group(1)) <= 5.0 + 1e-9, m.group(0)


@needs
def test_a_random_run_prints_its_seed_and_the_seed_repeats_it(tmp_path):
    """Diffusing tracers with random = true: the printed seed names the run
    folder, and random = false with that seed ends every tracer where the
    random run ended it, at the same thread count."""
    run_args = ["init_mode=points_x", "scheme=explicit", "Dm=0.01", "dt=0.01", "T=0.2",
                "stat_intv=0.1", "dump_intv=1e9", "num_threads=2"]
    first = linear_flow(tmp_path / "random", {"Axx": 1.0, "Ayy": -1.0})
    r = run_app(TRACERS, first, BASE, run_args, "random=true")
    m = re.search(r"Seed (\d+) drawn \(random=true\); random=false seed=(\d+) repeats this run", r.stdout)
    assert m and m.group(1) == m.group(2), r.stdout[-2000:]
    seed = m.group(1)
    run_dir = checkpoint_folder(tmp_path / "random")
    assert "_seed%s" % seed in str(run_dir), run_dir
    again = linear_flow(tmp_path / "again", {"Axx": 1.0, "Ayy": -1.0})
    r2 = run_app(TRACERS, again, BASE, run_args, "random=false seed=%s" % seed)
    assert "drawn" not in r2.stdout
    a = np.loadtxt(run_dir / "Checkpoints" / "positions.pos")
    b = np.loadtxt(checkpoint_folder(tmp_path / "again") / "Checkpoints" / "positions.pos")
    assert np.array_equal(a, b)
    assert not np.array_equal(a[:, 1], np.zeros(len(a)))   # the noise moved them off the axis
