"""Transported vectors (tracervectors) and tensors (tracertensors), against
closed forms.

A material line element rhohat obeys d(rhohat)/dt = J rhohat and the
deformation gradient dF/dt = J F, with J the velocity gradient along the
path. Where J is constant on the path there is a closed form:

  plane Poiseuille  u_z(x): x never changes, so J = J(x0) is constant and
                    nilpotent (J^2 = 0). F = I + J t, exactly, for every
                    scheme -- RK4 and explicit Euler both truncate at J^2.
  linear flow       u = A x: F = exp(A t). A stagnation point stretches by
                    e^t, a solid rotation turns without stretching.

The line element is carried as its unit direction rhohat, its log-stretch w
and the rate S = rhohat . J rhohat at the current position.
"""

import os

import numpy as np
import pytest

from dumps import deformation_gradient, dump_at
from paths import REPO, app
from runs import continuous_and_resumed, copy_example, run_app

VECTORS = app("tracervectors")
TENSORS = app("tracertensors")
POISEUILLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")
LINEAR = os.path.join(REPO, "data_example", "linear_flow", "expr_params.dat")

needs_partrac = pytest.mark.skipif(not (os.path.exists(VECTORS) and os.path.exists(TENSORS)),
                                   reason="tracervectors or tracertensors is not built")

# points_x places the particles on the x axis, and a line element starts
# along the init_mode's axis, so rhohat(0) = +-e_x exactly. tracervectors
# dumps rhohat as n.
BASE = ("mode=analytic init_mode=points_x x0=0 y0=0 z0=0 Nrw=50 Nrw_max=50 "
        "Dm=0 int_order=1 stat_intv=1e9 checkpoint_intv=1e9 random=false seed=1").split()

SCHEMES = ["RK4", "explicit"]


def app_for(transport):
    """The app that carries the element: tracertensors for a tensor, tracervectors for a vector."""
    return TENSORS if transport == "tensor" else VECTORS


def run(tmp_path, example, transport, extra, name="case"):
    """Run the app for transport on a copy of example with extra overriding BASE; return the case folder."""
    d = tmp_path / name
    run_app(app_for(transport), copy_example(example, d), BASE, extra)
    return d


def poiseuille_J(x):
    """Return the one nonzero entry of J in plane Poiseuille, dU_z/dx, at x."""
    return -3.0 * 1.0 * x / 1.0 ** 2       # u_inf = R = 1 in the example


# --- plane Poiseuille: exact for every scheme (J^2 = 0) --------------------------

@needs_partrac
@pytest.mark.parametrize("scheme", SCHEMES)
def test_a_line_element_in_plane_poiseuille(tmp_path, scheme):
    """A line element starting along x becomes +-(1, 0, J_zx T), exactly for
    both schemes, so w, S and rhohat match their closed forms to round-off.
    These are the quantities the stretching statistics are built from."""
    T = 0.5
    d = run(tmp_path, POISEUILLE, "vector", "scheme=%s dt=0.01 T=%g dump_intv=%g" % (scheme, T, T))
    t0, tT = dump_at(d, 0.), dump_at(d, T)
    x = t0["points"][:, 0]
    assert np.abs(tT["points"][:, 0] - x).max() == 0          # x is invariant
    assert np.abs(np.abs(t0["n"][:, 0]) - 1).max() == 0  # starts along +-x
    # el(T) = +-(1, 0, J_zx T): |el| and rhohat . J rhohat follow
    s = poiseuille_J(x) * T
    assert np.allclose(tT["w"][:, 0], 0.5 * np.log1p(s * s), rtol=0, atol=1e-12)
    assert np.allclose(tT["S"][:, 0], poiseuille_J(x) * s / (1 + s * s), rtol=0, atol=1e-12)
    expected = np.stack([np.ones_like(s), np.zeros_like(s), s], axis=1)
    expected /= np.linalg.norm(expected, axis=1)[:, None]
    assert np.allclose(np.abs(tT["n"]), np.abs(expected), rtol=0, atol=1e-12)


@needs_partrac
@pytest.mark.parametrize("scheme", SCHEMES)
def test_the_deformation_gradient_in_plane_poiseuille(tmp_path, scheme):
    """F starts as the identity and equals I + J(x0) T to round-off for both
    schemes, with the shear in the (z, x) entry. A transposed or misindexed
    dF/dt = J F would put the shear in the wrong place."""
    T = 0.5
    d = run(tmp_path, POISEUILLE, "tensor", "scheme=%s dt=0.01 T=%g dump_intv=%g" % (scheme, T, T))
    t0, tT = dump_at(d, 0.), dump_at(d, T)
    F0 = deformation_gradient(t0)
    assert np.array_equal(F0, np.tile(np.eye(3), (len(F0), 1, 1)))
    F = deformation_gradient(tT)
    expected = np.tile(np.eye(3), (len(F), 1, 1))
    expected[:, 2, 0] = poiseuille_J(t0["points"][:, 0]) * T
    assert np.abs(F - expected).max() < 1e-12


# --- linear flows: exp(A t) -----------------------------------------------------

# RK4 is fourth order, so with dt = 0.01 over T = 1 it meets the exponential
# well inside 1e-8. Explicit Euler at int_order=1 is first order: (1 + dt)^n
# differs from e^1 by about dt/2 relative, so it is checked to its own accuracy.
TOL = {"RK4": 1e-8, "explicit": 6e-3}


def linear(tmp_path, transport, scheme, A, T=1.0):
    """Run the linear flow u = A x with the given element and scheme to T; return (case folder, T)."""
    d = tmp_path / "flow"
    d.mkdir(parents=True)
    with open(LINEAR) as f:
        lines = [l for l in f.read().splitlines() if not l.startswith(("A", "#"))]
    (d / "expr_params.dat").write_text("\n".join(lines) + "\n" + "".join(
        "A%s%s=%r\n" % ("xyz"[i], "xyz"[j], A[i][j]) for i in range(3) for j in range(3)))
    return run(tmp_path, str(d / "expr_params.dat"), transport,
               "scheme=%s dt=0.01 T=%g dump_intv=%g" % (scheme, T, T)), T


STAGNATION = [[1., 0., 0.], [0., -1., 0.], [0., 0., 0.]]
ROTATION = [[0., -1., 0.], [1., 0., 0.], [0., 0., 0.]]


@needs_partrac
@pytest.mark.parametrize("scheme", SCHEMES)
def test_a_line_element_at_a_stagnation_point(tmp_path, scheme):
    """In u = (x, -y, 0) a particle on the x axis and a line element along it
    grow as e^T, so w = T, S is the eigenvalue 1 and rhohat stays e_x. This
    checks exponential stretching to each scheme's order of accuracy."""
    d, T = linear(tmp_path, "vector", scheme, STAGNATION)
    t0, tT = dump_at(d, 0.), dump_at(d, T)
    # along the stretching axis: |el| = e^T, the rate is the eigenvalue 1
    assert np.allclose(tT["points"][:, 0], t0["points"][:, 0] * np.exp(T), rtol=TOL[scheme])
    assert np.allclose(tT["w"][:, 0], T, rtol=TOL[scheme])
    assert np.allclose(tT["S"][:, 0], 1.0, rtol=0, atol=1e-12)
    assert np.allclose(np.abs(tT["n"]), [[1., 0., 0.]], rtol=0, atol=1e-15)


@needs_partrac
@pytest.mark.parametrize("scheme", SCHEMES)
def test_the_deformation_gradient_at_a_stagnation_point(tmp_path, scheme):
    """F = exp(A T) = diag(e^T, e^-T, 1) at the stagnation point, to each
    scheme's accuracy. Unlike plane Poiseuille, F is not a polynomial in t
    here, so the schemes' truncation error shows."""
    d, T = linear(tmp_path, "tensor", scheme, STAGNATION)
    F = deformation_gradient(dump_at(d, T))
    expected = np.diag([np.exp(T), np.exp(-T), 1.0])
    assert np.abs(F - expected).max() < TOL[scheme] * np.exp(T)


@needs_partrac
def test_solid_rotation_turns_without_stretching(tmp_path):
    """In solid rotation a line element turns with the flow, rhohat =
    (cos T, sin T, 0), with no stretch (w = 0, S = 0), and F is the rotation
    matrix with det F = 1. A transport that stretched under pure rotation would
    report spurious stretching wherever the flow rotates."""
    dv, T = linear(tmp_path / "v", "vector", "RK4", ROTATION)
    tT = dump_at(dv, T)
    assert np.abs(tT["w"]).max() < 1e-8
    assert np.abs(tT["S"]).max() < 1e-12                 # rhohat . A rhohat, A antisymmetric
    assert np.allclose(np.abs(tT["n"]), [[np.cos(T), np.sin(T), 0.]], rtol=0, atol=1e-8)
    dt_, T = linear(tmp_path / "t", "tensor", "RK4", ROTATION)
    F = deformation_gradient(dump_at(dt_, T))
    R = np.array([[np.cos(T), -np.sin(T), 0.], [np.sin(T), np.cos(T), 0.], [0., 0., 1.]])
    assert np.abs(F - R).max() < 1e-8
    assert np.abs(np.linalg.det(F) - 1).max() < 1e-8


# --- restart -------------------------------------------------------------------

@needs_partrac
@pytest.mark.parametrize("fmt", ["hdf5", "text"])
@pytest.mark.parametrize("transport", ["vector", "tensor"])
def test_a_resumed_run_is_identical_to_one_never_stopped(tmp_path, transport, fmt):
    """A run stopped at t = 0.2 and resumed from its checkpoint writes a dump
    at t = 0.4 bit-identical to an uninterrupted run. The checkpoint must carry
    the element (rhohat, w, S or F's factors) as well as the position, or a resumed run
    loses the deformation accumulated before the restart. An old text
    checkpoint (Q.ten, logstretch.vec, U.vec) resumes bit for bit as well."""
    dt, stop, end = 0.01, 0.2, 0.4
    # the final checkpoint is written one step past T, so T = stop - dt checkpoints at stop
    cont, split = continuous_and_resumed(
        app_for(transport), POISEUILLE, tmp_path,
        [BASE, "scheme=RK4 dt=%g dump_intv=%g" % (dt, stop)],
        "T=%g" % (stop - dt), "T=%g" % end, text=fmt == "text")
    a, b = dump_at(cont, end), dump_at(split, end)
    assert set(a) == set(b)
    for k in a:
        assert np.array_equal(a[k], b[k]), k + " differs after a restart"


# --- F as its factors ------------------------------------------------------------

# u = A x with A^2 = I, eigenvalues +-1 along directions off the axes: every
# entry of F = cosh(t) I + sinh(t) A grows as e^t/2, so det F = 1 is a
# difference of numbers near e^(2t)/4, and at T = 30 no digit of it survives
# in F itself. Its factors keep the compressed direction.
MIXING = [[0., 2., 0.], [0.5, 0., 0.], [0., 0., 0.]]


@needs_partrac
def test_the_factors_keep_the_compressed_direction(tmp_path):
    """At T = 30 the first stretch is log |F e_x| and the second its negative,
    so log det F = 0 to 1e-8, where F itself has lost it (its determinant,
    computed from the entries, is off by more than one)."""
    d, T = linear(tmp_path, "tensor", "RK4", MIXING, T=30.0)
    g = dump_at(d, T)
    s = g["logstretch"].reshape(-1, 3)
    col0 = np.array([np.cosh(T), 0.5 * np.sinh(T), 0.])
    assert np.allclose(s[:, 0], np.log(np.linalg.norm(col0)), rtol=1e-8)
    assert np.abs(s[:, 0] + s[:, 1]).max() < 1e-8
    assert np.abs(s[:, 2]).max() < 1e-12
    F = deformation_gradient(g)
    assert np.abs(np.linalg.det(F) - 1).max() > 1        # the whole F is useless here


@needs_partrac
def test_a_checkpoint_holding_F_whole_resumes(tmp_path):
    """A text checkpoint from before the factors holds F whole in F.ten: a run
    resumed from one, made here from a new checkpoint's factors, ends where a
    run never stopped does, to round-off (F is factored on load)."""
    from runs import checkpoint_folder, write_text_checkpoint
    dt, stop, end = 0.01, 0.2, 0.4
    args = [BASE, "scheme=RK4 dt=%g dump_intv=%g" % (dt, stop)]
    cont = copy_example(POISEUILLE, tmp_path / "cont")
    old = copy_example(POISEUILLE, tmp_path / "old")
    run_app(TENSORS, cont, args, "T=%g" % end)
    run_app(TENSORS, old, args, "T=%g" % (stop - dt))
    ck = checkpoint_folder(tmp_path / "old") / "Checkpoints"
    write_text_checkpoint(ck, F_whole=True)
    assert (ck / "F.ten").exists() and not (ck / "Q.ten").exists()
    run_app(TENSORS, old, args, "T=%g" % end, "restart_folder=%s" % ck.parent)
    a, b = dump_at(tmp_path / "cont", end), dump_at(tmp_path / "old", end)
    assert np.abs(deformation_gradient(a) - deformation_gradient(b)).max() < 1e-12
    assert np.array_equal(a["points"], b["points"])


@needs_partrac
def test_a_plane_flow_stretches_past_the_range_of_a_double(tmp_path):
    """On the z axis of u = (x, -y, 0), at rest, the stretches grow as +-t without
    bound. At T = 800, where e^800 overflows a double, the factors stay finite
    and exact in their logs: the frame's third column, e_z in a plane flow,
    adds nothing to U however large the ratio of stretches it would scale."""
    d = tmp_path / "flow"
    d.mkdir(parents=True)
    with open(LINEAR) as f:
        lines = [l for l in f.read().splitlines() if not l.startswith(("A", "#"))]
    (d / "expr_params.dat").write_text("\n".join(lines) + "\nAxx=1.0\nAyy=-1.0\n")
    T = 800.0
    case = run(tmp_path, str(d / "expr_params.dat"), "tensor",
               "init_mode=points_z Nrw=10 Nrw_max=10 scheme=RK4 dt=0.1 T=%g dump_intv=%g" % (T, T))
    g = dump_at(case, T)
    s, u = g["logstretch"].reshape(-1, 3), g["U"].reshape(-1, 3)
    assert np.isfinite(s).all() and np.isfinite(u).all() and np.isfinite(g["Q"]).all()
    assert np.allclose(s[:, 0], T, rtol=1e-5) and np.allclose(s[:, 1], -T, rtol=1e-5)
    assert np.abs(s[:, 0] + s[:, 1]).max() < 1e-3       # RK4's map has det 1 - O(dt^6) a step
    assert np.array_equal(s[:, 2], np.zeros(len(s))) and np.array_equal(u, np.zeros_like(u))
