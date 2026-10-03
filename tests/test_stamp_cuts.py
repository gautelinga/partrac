"""Steps cut at the field's timestamps.

A stamped field blends linearly between its stamps, so its time derivative
jumps at every stamp. The run loop cuts a step with a stamp inside into pieces
that meet at the stamp, each on its own bracket, so a step never uses the old
bracket's rate past a stamp and RK4 keeps its order.

The fields are on a P1 triangle mesh of the square [-2, 2]^2 and are exact
there: u = A x, with J = A in every cell, and uniform fields u = U. The stamps
alternate between two of them at times off the dt grid, so u(t) = A(t) x with
A(t) piecewise linear in time, and the deformation gradient is the flow map
Phi(t), the same for every particle, with x(t) = Phi(t) x0. Phi is computed
here with fine RK4 steps that meet the stamps.

partrac's n_accepted column counts every particle's accepted step, a piece
being a step, so it tells how many pieces a run took.
"""

import bisect
import os
import re

import numpy as np
import pytest

from cases import LINEAR, UNIFORM, stamp_case
from dumps import all_dumps, deformation_gradient, dump_at, read_stats
from paths import app
from runs import checkpoint_file, checkpoint_folder, run_app

PARTRAC = app("partrac")
TRACERS = app("tracers")
TENSORS = app("tracertensors")

needs_apps = pytest.mark.skipif(not all(os.path.exists(a) for a in (PARTRAC, TRACERS, TENSORS)),
                                reason="partrac, tracers or tracertensors is not built")

# Alternating stamps off the dt grid of every dt below; the last far past T
ALTERNATING = [(0.0, "a"), (0.37, "b"), (0.61, "a"), (1.13, "b"), (3.0, "a")]
T = 1.2

TENSOR_BASE = ("mode=triangle init_mode=points_x x0=0 y0=0 z0=0 Nrw=40 Nrw_max=40 "
               "Dm=0 int_order=1 stat_intv=1e9 checkpoint_intv=1e9 random=false seed=1 outside=ignore").split()

PARTRAC_BASE = ("mode=triangle init_mode=strip_x La=1 x0=0 y0=0 z0=0 Nrw=50 Nrw_max=5000 "
                "ds_max=1e9 ds_min=1e-9 refine=false coarsen=false Dm=0 int_order=1 scheme=RK4 "
                "checkpoint_intv=1e9 random=false seed=1").split()


def alternating(prefix):
    """ALTERNATING with the fields named prefix_a and prefix_b."""
    return [(t, prefix + "_" + f) for t, f in ALTERNATING]


def blended(stamps, fields):
    """f(t): the stamps' fields blended linearly in time, as the loaders do."""
    ts = [t for t, _ in stamps]
    vals = [np.array(fields[f], dtype=float) for _, f in stamps]

    def at(t):
        i = min(max(bisect.bisect_right(ts, t), 1), len(ts) - 1)
        w = (t - ts[i - 1]) / (ts[i] - ts[i - 1])
        return (1 - w) * vals[i - 1] + w * vals[i]
    return at


def flow_map(stamps, t_end, n=400):
    """Phi(t_end) for dPhi/dt = A(t) Phi, Phi(0) = I: RK4 with n steps between stamps."""
    A = blended(stamps, LINEAR)
    knots = sorted({0.0, t_end} | {t for t, _ in stamps if 0 < t < t_end})
    P = np.eye(2)
    for a, b in zip(knots[:-1], knots[1:]):
        h = (b - a) / n
        for k in range(n):
            t = a + k * h
            k1 = A(t) @ P
            k2 = A(t + h / 2) @ (P + h / 2 * k1)
            k3 = A(t + h / 2) @ (P + h / 2 * k2)
            k4 = A(t + h) @ (P + h * k3)
            P = P + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return P


def tensor_errors(stamp_mesh, root, dt, binary=TENSORS):
    """(max |F - Phi|, max |x - Phi x0|) at T of tracertensors under RK4 on the alternating linear flow."""
    params = stamp_case(stamp_mesh, root, alternating("lin"), "dt%g" % dt)
    run_app(binary, params, TENSOR_BASE, "scheme=RK4 dt=%r T=%r dump_intv=%r" % (dt, T, T + dt / 2))
    g0, gT = dump_at(params.parent, 0.0), dump_at(params.parent, T)
    # far from the walls, where no particle leaves the square
    keep = np.abs(g0["points"][:, 0]) < 0.8
    assert keep.sum() > 5
    Phi = flow_map(alternating("lin"), T)
    F = deformation_gradient(gT)[keep]
    eF = np.abs(F[:, :2, :2] - Phi).max()
    x = gT["points"][keep, :2]
    ex = np.abs(x - g0["points"][keep, :2] @ Phi.T).max()
    return eF, ex


def n_accepted(case_dir):
    """The last n_accepted of a partrac run under case_dir."""
    return int(read_stats(case_dir)["n_accepted"][-1])


@needs_apps
def test_rk4_keeps_its_order_across_stamps(stamp_mesh, tmp_path):
    """F and x converge at fourth order in dt with stamps inside steps: each
    halving of dt divides the error by about 16. A step that used the old
    bracket's rate past a stamp would leave a kink in A(t) inside it and fall
    to second order (a ratio of 4), as before steps were cut."""
    errs = [tensor_errors(stamp_mesh, tmp_path, dt) for dt in (0.1, 0.05, 0.025)]
    for (eF0, ex0), (eF1, ex1) in zip(errs[:-1], errs[1:]):
        assert eF0 / eF1 > 12, errs
        assert ex0 / ex1 > 12, errs
    assert errs[-1][0] < 1e-6


@needs_apps
def test_the_explicit_step_takes_the_new_rate_past_a_stamp(stamp_mesh, tmp_path):
    """Uniform fields alternating at stamps off the dt grid: the explicit step
    at int_order=2 adds the stamps' rate, so it integrates a velocity linear
    in time exactly, and across stamps only if the step is cut there. Every
    particle moves by the integral of U(t), to round-off."""
    stamps = alternating("uni")
    params = stamp_case(stamp_mesh, tmp_path, stamps)
    run_app(TENSORS, params, TENSOR_BASE,
            "scheme=explicit int_order=2 dt=0.1 T=%r dump_intv=%r" % (T, T + 0.05))
    g0, gT = dump_at(params.parent, 0.0), dump_at(params.parent, T)
    keep = np.abs(g0["points"][:, 0]) < 1.5
    # the integral of a piecewise linear U(t): trapezoids between the knots
    U = blended(stamps, UNIFORM)
    knots = sorted({0.0, T} | {t for t, _ in stamps if 0 < t < T})
    shift = sum((b - a) / 2 * (U(a) + U(b)) for a, b in zip(knots[:-1], knots[1:]))
    assert np.abs(gT["points"][keep, :2] - g0["points"][keep, :2] - shift).max() < 1e-12


def partrac_case(stamp_mesh, root, stamps, extra="", name="case", dt=0.1, T_end=0.5):
    """Run partrac on the strip to T_end with the stamps; return its case folder."""
    params = stamp_case(stamp_mesh, root, stamps, name)
    run_app(PARTRAC, params, PARTRAC_BASE,
            "dt=%r T=%r stat_intv=%r dump_intv=%r" % (dt, T_end, dt, T_end + dt / 2), extra)
    return params.parent


def steps_below(dt, n):
    """The first k < n with k dt just below the decimal it stands for, and that decimal."""
    for k in range(1, n):
        s = round(k * dt, 12)
        if k * dt < s:
            return k, s
    raise AssertionError("no step count below its decimal")


@needs_apps
def test_a_stamp_inside_a_step_cuts_it(stamp_mesh, tmp_path):
    """A stamp at 0.35 inside the step from 0.3 adds one piece: five steps of
    50 particles take 300 accepted steps, not 250."""
    d = partrac_case(stamp_mesh, tmp_path, [(0.0, "lin_a"), (0.35, "lin_b"), (3.0, "lin_a")])
    assert n_accepted(d) == 300


@needs_apps
def test_a_stamp_just_above_a_step_start_snaps_to_it(stamp_mesh, tmp_path):
    """t = k dt lands just below the decimal k dt, where a stamp sits: the step
    from there runs whole on the new bracket (no sliver piece) and agrees with
    a run whose stamp is a little below that t. Without the snap the whole
    step would run on the old bracket's rate, an error of order dt^2."""
    dt = 0.09
    k, s = steps_below(dt, 20)
    T_end = round((k + 2) * dt, 12)
    snapped = partrac_case(stamp_mesh, tmp_path, [(0.0, "lin_a"), (s, "lin_b"), (3.0, "lin_a")],
                           name="snapped", dt=dt, T_end=T_end)
    below = partrac_case(stamp_mesh, tmp_path, [(0.0, "lin_a"), (s - 1e-7, "lin_b"), (3.0, "lin_a")],
                         name="below", dt=dt, T_end=T_end)
    steps = k + 2
    assert n_accepted(snapped) == 50 * steps
    # the stamp below t is inside the step before: one more piece there
    assert n_accepted(below) == 50 * (steps + 1)
    a, b = dump_at(snapped, T_end)["points"], dump_at(below, T_end)["points"]
    assert np.abs(a - b).max() < 1e-6


@needs_apps
@pytest.mark.parametrize("offset", [1e-12, -1e-12])
def test_a_stamp_within_snap_of_a_step_end_is_not_a_piece(stamp_mesh, tmp_path, offset):
    """A stamp 1e-12 above the start of the step from 0.3, or 1e-12 below the
    end of the step to it, makes no sliver piece: five steps, 250 accepted
    steps, and the same positions as with the stamp at 0.3 exactly."""
    s = 3 * 0.1
    near = partrac_case(stamp_mesh, tmp_path, [(0.0, "lin_a"), (s + offset, "lin_b"), (3.0, "lin_a")],
                        name="near")
    exact = partrac_case(stamp_mesh, tmp_path, [(0.0, "lin_a"), (s, "lin_b"), (3.0, "lin_a")],
                         name="exact")
    assert n_accepted(near) == 250
    a, b = dump_at(near, 0.5)["points"], dump_at(exact, 0.5)["points"]
    assert np.abs(a - b).max() < 1e-9


@needs_apps
def test_frozen_fields_are_not_cut(stamp_mesh, tmp_path):
    """With frozen fields the stamps do not matter after the start: no piece is cut."""
    d = partrac_case(stamp_mesh, tmp_path, [(0.0, "lin_a"), (0.35, "lin_b"), (3.0, "lin_a")],
                     "frozen_fields=true t_frozen=0.2")
    assert n_accepted(d) == 250


@needs_apps
def test_tau_is_integrated_once_a_step_in_a_cut_step(stamp_mesh, tmp_path):
    """The same steady field with and without a stamp inside a step: the cut
    run takes one more piece, and its compressed time tau, integrated every
    step, is the uncut run's to RK4's step error; integrated after every piece
    it would count the cut step twice, a difference of order 0.1."""
    extra = "integrate_tau=true tau_intv=0.1 tau_max=0"
    cut = partrac_case(stamp_mesh, tmp_path, [(0.0, "lin_a"), (0.35, "lin_a"), (3.0, "lin_a")],
                       extra, name="cut")
    whole = partrac_case(stamp_mesh, tmp_path, [(0.0, "lin_a"), (3.0, "lin_a")], extra, name="whole")
    assert n_accepted(cut) == 300 and n_accepted(whole) == 250
    a, b = dump_at(cut, 0.5)["tau"], dump_at(whole, 0.5)["tau"]
    assert np.abs(b).max() > 0.1
    assert np.abs(a - b).max() < 1e-6


def outside_run(stamp_mesh, root, outside):
    """tracers on the outflow u = (x, -y), with a stamp inside the one step [0, 0.1]; return (stdout, n_declined)."""
    params = stamp_case(stamp_mesh, root, [(0.0, "out"), (0.05, "out"), (3.0, "out")], outside)
    r = run_app(TRACERS, params, TENSOR_BASE, "Nrw=200 Nrw_max=200 scheme=RK4 dt=0.1 T=0.1",
                "stat_intv=0.1 dump_intv=1e9 verbose=true outside=%s" % outside)
    return r.stdout, int(read_stats(params.parent)["n_declined"][-1])


@needs_apps
def test_a_particle_declined_in_a_piece_is_handled_before_the_next(stamp_mesh, tmp_path):
    """Particles near x = +-2 leave the square in the first piece of the step
    [0, 0.1], cut at 0.05. The outside rule runs after that piece: marked
    particles are declined again in the second piece (reported at t = 0.05),
    reinjected ones move on in it, so fewer steps are declined."""
    out_mark, declined_mark = outside_run(stamp_mesh, tmp_path, "mark")
    out_reinject, declined_reinject = outside_run(stamp_mesh, tmp_path, "reinject")
    first = re.findall(r"(\d+) nodes could not move at t = 0,", out_mark)
    second = re.findall(r"(\d+) nodes could not move at t = 0.05,", out_mark)
    assert first and second, out_mark
    assert int(second[0]) >= int(first[0]) > 0
    assert declined_mark == int(first[0]) + int(second[0])
    assert re.search(r"nodes could not move at t = 0,", out_reinject), out_reinject
    assert declined_reinject < declined_mark


@needs_apps
def test_a_resumed_run_has_the_times_of_one_never_stopped(stamp_mesh, tmp_path):
    """t comes from the step count since the run's start, so a run resumed from
    a checkpoint at 0.7 steps at the same times, to the last bit, as one never
    stopped: on the alternating field the positions depend on t's last bits
    through the blend, and they agree bit for bit at every dump, as does the
    final checkpoint's t."""
    stamps = alternating("lin")
    common = ("scheme=RK4 dt=0.1 dump_intv=0.1 stat_intv=1e9 init_mode=points_x Nrw=40 "
              "Nrw_max=40")
    end = ["T=1.1", "checkpoint_intv=1e9"]
    cont = stamp_case(stamp_mesh, tmp_path, stamps, "cont")
    split = stamp_case(stamp_mesh, tmp_path, stamps, "split")
    run_app(TENSORS, cont, TENSOR_BASE, common, end)
    run_app(TENSORS, split, TENSOR_BASE, common, ["T=0.65"])
    resume = checkpoint_folder(split.parent)
    run_app(TENSORS, split, TENSOR_BASE, common, end, ["restart_folder=%s" % resume])

    h5py = pytest.importorskip("h5py")
    a = all_dumps(cont.parent, raw=True)
    b = all_dumps(split.parent, raw=True)
    later = [t for t in a if t > 0.7]
    assert len(later) >= 4 and set(later) <= set(b)
    for t in later:
        assert np.array_equal(a[t]["points"], b[t]["points"]), t
    files = [checkpoint_file(cont.parent), checkpoint_file(split.parent)]
    ts = []
    for f in files:
        with h5py.File(f, "r") as h:
            ts.append(h.attrs["t"])
    # twelve steps: 12 * 0.1, an ulp above the 1.2 that adding 0.1 twelve times gives
    assert ts[0] == ts[1] == 12 * 0.1
