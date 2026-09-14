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
import shutil
import subprocess

import numpy as np
import pytest

from dumps import dump_at
from paths import REPO, app

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


def run(tmp_path, example, extra, name="case", env=None):
    """Run tracertensors if extra has transport=tensor, else tracervectors, on a copy of example; return the case folder."""
    d = tmp_path / name
    d.mkdir(parents=True)
    shutil.copy(example, d / "expr_params.dat")
    if isinstance(extra, str):
        extra = extra.split()
    # transport= names the app: the element is what each one carries
    binary = TENSORS if "transport=tensor" in extra else VECTORS
    extra = [a for a in extra if not a.startswith("transport=")]
    keys = {a.split("=")[0] for a in extra}
    argv = [a for a in BASE if a.split("=")[0] not in keys] + list(extra)
    r = subprocess.run([binary, str(d / "expr_params.dat")] + argv,
                       capture_output=True, text=True, timeout=900,
                       env=dict(os.environ, **(env or {})))
    assert r.returncode == 0, r.stdout + r.stderr
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
    d = run(tmp_path, POISEUILLE, "transport=vector scheme=%s dt=0.01 T=%g dump_intv=%g" % (scheme, T, T))
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
    d = run(tmp_path, POISEUILLE, "transport=tensor scheme=%s dt=0.01 T=%g dump_intv=%g" % (scheme, T, T))
    t0, tT = dump_at(d, 0.), dump_at(d, T)
    assert np.array_equal(t0["F"], np.tile(np.eye(3).ravel(), (len(t0["F"]), 1)))
    F = tT["F"].reshape(-1, 3, 3)
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
    lines = [l for l in open(LINEAR).read().splitlines()
             if not l.startswith(("A", "#"))]
    (d / "expr_params.dat").write_text("\n".join(lines) + "\n" + "".join(
        "A%s%s=%r\n" % ("xyz"[i], "xyz"[j], A[i][j]) for i in range(3) for j in range(3)))
    return run(tmp_path, str(d / "expr_params.dat"),
               "transport=%s scheme=%s dt=0.01 T=%g dump_intv=%g" % (transport, scheme, T, T)), T


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
    F = dump_at(d, T)["F"].reshape(-1, 3, 3)
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
    F = dump_at(dt_, T)["F"].reshape(-1, 3, 3)
    R = np.array([[np.cos(T), -np.sin(T), 0.], [np.sin(T), np.cos(T), 0.], [0., 0., 1.]])
    assert np.abs(F - R).max() < 1e-8
    assert np.abs(np.linalg.det(F) - 1).max() < 1e-8


# --- restart and threads -----------------------------------------------------

@needs_partrac
@pytest.mark.parametrize("transport", ["vector", "tensor"])
def test_a_resumed_run_is_identical_to_one_never_stopped(tmp_path, transport):
    """A run stopped at t = 0.2 and resumed from its checkpoint writes a dump
    at t = 0.4 bit-identical to an uninterrupted run. The checkpoint must carry
    the element (rhohat, w, S or F) as well as the position, or a resumed run
    loses the deformation accumulated before the restart."""
    dt, stop, end = 0.01, 0.2, 0.4
    common = "transport=%s scheme=RK4 dt=%g dump_intv=%g" % (transport, dt, stop)
    cont = run(tmp_path, POISEUILLE, common + " T=%g" % end, "cont")
    # the final checkpoint is written one step past T, so T = stop - dt checkpoints at stop
    split = run(tmp_path, POISEUILLE, common + " T=%g" % (stop - dt), "split")
    folder = os.path.dirname(next(split.rglob("Checkpoints/positions.pos")))
    run(tmp_path, POISEUILLE, common + " T=%g restart_folder=%s" % (end, os.path.dirname(folder)), "resume")
    a, b = dump_at(cont, end), dump_at(split, end)
    assert set(a) == set(b)
    for k in a:
        assert np.array_equal(a[k], b[k]), k + " differs after a restart"


@needs_partrac
@pytest.mark.parametrize("transport", ["vector", "tensor"])
def test_the_result_does_not_depend_on_the_thread_count(tmp_path, transport):
    """Runs on 1 and 4 threads write bit-identical dumps for both elements.
    Each particle's element evolves independently, so the result must not
    depend on how the particles are split across threads."""
    out = {}
    for n in (1, 4):
        d = run(tmp_path, POISEUILLE, "transport=%s scheme=RK4 dt=0.01 T=0.2 dump_intv=0.2 Nrw=400 Nrw_max=400" % transport,
                str(n), env={"OMP_NUM_THREADS": str(n)})
        out[n] = dump_at(d, 0.2)
    assert set(out[1]) == set(out[4])
    for k in out[1]:
        assert np.array_equal(out[1][k], out[4][k]), k + " differs between 1 and 4 threads"
