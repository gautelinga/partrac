"""StructuredInterpol across a change of timestamp.

A synthetic felbm case on a 16^3 grid, with solid walls at z = 0 and z = 15,
has three stamps at t = 0, 1, 2. Stamp k holds the uniform field u_z = k, so
with linear blending in time u_z(t) = t and z(t) = z0 + t^2/2 exactly. The
mean u_z over [1, 1.5] is 1.25 only if, at the stamp change, the old next
stamp becomes the previous one for every component and the two are blended;
z at the end checks the blend itself. With shear, stamp k also holds
u_y = k x, so the whole velocity gradient is linear in time: at t = 0.5 it is
exactly half of its value at t = 1, also in the cells next to a wall that
only int_order=2 reaches. With interpolation=constant each point takes the
nearest node's values, so u_y is k times the nearest node's x.
"""

import os

import numpy as np
import pytest

from cases import write_felbm
from dumps import all_dumps
from paths import app
from runs import run_app

FILAMENTS = app("filaments")
INTERPOL = app("interpol")
TRACERS = app("tracers")


def three_stamp_felbm(d, shear=False, extra=""):
    """Write a felbm case in d where stamp k holds u_z = k, and with shear also u_y = k x.

    extra is appended to felbm_params.dat."""
    pytest.importorskip("h5py")
    n = 16
    zero = np.zeros((n, n, n))
    x = np.arange(n, dtype=float)[:, None, None] * np.ones((n, n, n))
    solid = np.zeros((n, n, n), dtype=np.int32)
    solid[0, :, :] = 1
    solid[-1, :, :] = 1
    fields = [{"u_x": zero, "u_y": k * x if shear else zero,
               "u_z": np.full((n, n, n), float(k)),
               "density": np.ones((n, n, n)), "pressure": zero} for k in range(3)]
    write_felbm(d, fields, solid, times=(0, 1, 2), extra=extra)


@pytest.mark.skipif(not os.path.exists(FILAMENTS), reason="filaments is not built")
def test_uz_advances_with_the_timestamp(tmp_path):
    """Past stamp 1 the field is blended between stamps 1 and 2, so particles
    move at u_z > 0.9 over [1, 1.5] and z(T) follows t^2/2. If a component
    kept an older stamp, or the field were held constant between stamps, every
    structured run would advect with a field that lags the data in time."""
    pytest.importorskip("h5py")
    d = tmp_path / "felbm"
    d.mkdir()
    three_stamp_felbm(d)
    run_app(FILAMENTS, d / "felbm_params.dat",
            "mode=felbm scheme=RK4 resize=doublings resize_target=ds_init outside=reinject "
            "Dm=0 dt=0.01 T=2.0 Nrw=2 Nrw_max=100 dump_intv=0.5 stat_intv=1e9 "
            "checkpoint_intv=1e9 init_mode=pairs_xy int_order=1 ds_init=0.5 "
            "x0=8 y0=8 z0=8 random=false seed=1", timeout=600)
    assert len(list(d.rglob("data_from_t*.h5"))) == 1
    z = {t: g["points"][:, 2].mean() for t, g in all_dumps(d, raw=True).items()}
    assert 1.0 in z and 1.5 in z, sorted(z)
    u_z_after = (z[1.5] - z[1.0]) / 0.5
    assert u_z_after > 0.9, "u_z over [1, 1.5] is %.3f: stamp 1 did not become prev" % u_z_after
    # z(t) = t^2/2; first-order Euler at dt = 0.01 falls short by t*dt/2, well inside 0.05
    t = max(z)
    assert abs((z[t] - 8.0) - t * t / 2) < 0.05, (t, z[t] - 8.0, t * t / 2)


@pytest.mark.skipif(not os.path.exists(TRACERS), reason="tracers is not built")
def test_frozen_fields_hold_the_blend_at_t_frozen(tmp_path):
    """Frozen at t = 0.5, the field is u_z = 0.5 at every time, also past
    the last stamp, and the tracers move at that speed."""
    d = tmp_path / "felbm"
    d.mkdir()
    three_stamp_felbm(d)
    run_app(TRACERS, d / "felbm_params.dat",
            "mode=felbm init_mode=points_xyz x0=8 y0=8 z0=8 Nrw=4 Nrw_max=4 Dm=0 int_order=1 dt=0.1 T=3 "
            "dump_intv=1 stat_intv=1e9 checkpoint_intv=1e9 frozen_fields=true t_frozen=0.5 random=false seed=1")
    dumps = all_dumps(d, raw=True)
    assert sorted(dumps) == [0.0, 1.0, 2.0, 3.0]
    for g in dumps.values():
        assert np.allclose(g["u"][:, 2], 0.5, rtol=0, atol=1e-12)
    assert np.allclose(dumps[3.0]["points"][:, 2] - dumps[0.0]["points"][:, 2], 1.5, rtol=0, atol=1e-9)


def run_probe(d, t0, int_order, check=True):
    """Run interpol on the case in d at time t0; return the process."""
    return run_app(INTERPOL, d / "felbm_params.dat",
                   "mode=felbm Nrw=5000 int_order=%d t0=%g random=false seed=1" % (int_order, t0),
                   check=check, timeout=600)


def probe(d, t0, int_order=2):
    """Run interpol on the case in d at time t0; return the probed datasets by name."""
    h5py = pytest.importorskip("h5py")
    run_probe(d, t0, int_order)
    f = list(d.rglob("interpolation.h5part"))
    assert len(f) == 1
    with h5py.File(f[0], "r") as h:
        return {k: np.array(h["Step#0"][k]) for k in h["Step#0"]}


@pytest.mark.skipif(not os.path.exists(INTERPOL), reason="interpol is not built")
def test_gradient_blends_between_stamps(tmp_path):
    """Stamp 0 is at rest, so halfway to stamp 1 every component of the
    velocity gradient is exactly half its value at stamp 1, in every cell
    including those next to a wall. Otherwise gradU and gradA in second-order
    runs near solids would mix derivatives from different stamps."""
    a = tmp_path / "half"
    b = tmp_path / "whole"
    for d in (a, b):
        d.mkdir()
        three_stamp_felbm(d, shear=True)
    half = probe(a, 0.5)
    whole = probe(b, 1.0)
    assert np.array_equal(half["x"], whole["x"])
    for c in ("uxx", "uxy", "uxz", "uyx", "uyy", "uyz", "uzx", "uzy", "uzz"):
        assert np.allclose(half[c], 0.5 * whole[c], rtol=0, atol=1e-12), c
    assert np.abs(whole["uyx"]).max() > 0.5   # the shear is there, so the check above is not 0 == 0


@pytest.mark.skipif(not os.path.exists(INTERPOL), reason="interpol is not built")
def test_constant_takes_the_nearest_node(tmp_path):
    """interpolation=constant: halfway between stamps 0 and 1 every point has
    u_y = x_n/2 and u_z = 1/2, x_n the x of its nearest node (unit spacing,
    periodic, so x_n = 16 is node 0), and no gradient. Trilinear
    interpolation would give u_y = x/2 instead, so the check fails for
    either a trilinear evaluate or a wrong node."""
    d = tmp_path / "constant"
    d.mkdir()
    three_stamp_felbm(d, shear=True, extra="interpolation=constant\n")
    v = probe(d, 0.5, int_order=1)
    assert len(v["x"]) > 1000
    nearest = np.floor(v["x"] + 0.5) % 16
    assert np.array_equal(v["uy"], 0.5 * nearest)
    assert np.array_equal(v["uz"], np.full_like(v["uz"], 0.5))
    assert np.array_equal(v["ux"], np.zeros_like(v["ux"]))
    assert np.array_equal(v["rho"], np.ones_like(v["rho"]))
    for c in ("uxx", "uxy", "uxz", "uyx", "uyy", "uyz", "uzx", "uzy", "uzz"):
        assert not v[c].any(), c
    assert np.abs(v["uy"] - 0.5 * v["x"]).max() > 0.2   # not the trilinear field


@pytest.mark.skipif(not os.path.exists(INTERPOL), reason="interpol is not built")
def test_constant_refuses_a_gradient(tmp_path):
    """A piecewise-constant field has no gradient, so a run that needs one
    (int_order=2 here; vectors and tensors too) stops at setup with the
    parameter-error code instead of running on a zero gradient."""
    d = tmp_path / "constant"
    d.mkdir()
    three_stamp_felbm(d, shear=True, extra="interpolation=constant\n")
    r = run_probe(d, 0.5, int_order=2, check=False)
    assert r.returncode == 2, r.stdout + r.stderr
    assert "interpolation=constant" in r.stderr


@pytest.mark.skipif(not os.path.exists(INTERPOL), reason="interpol is not built")
def test_linear_is_the_default(tmp_path):
    """interpolation=linear, written out, gives the same probe as a file
    without the key, which is how every existing felbm file runs."""
    a = tmp_path / "default"
    b = tmp_path / "linear"
    for d, extra in ((a, ""), (b, "interpolation=linear\n")):
        d.mkdir()
        three_stamp_felbm(d, shear=True, extra=extra)
    va, vb = probe(a, 0.5), probe(b, 0.5)
    assert va.keys() == vb.keys()
    for k in va:
        assert np.array_equal(va[k], vb[k]), k


TRACERS = app("tracers")


def wall_felbm(d):
    """A 16^3 felbm case with walls at z = 0, 15 and a pillar 6 <= x, y <= 9,
    holding a steady field whose three components differ and are not linear."""
    n = 16
    i = np.arange(n, dtype=float)
    x, y, z = np.meshgrid(i, i, i, indexing="ij")
    k = 2 * np.pi / n
    solid = np.zeros((n, n, n), dtype=np.int32)   # [z, y, x]
    solid[0, :, :] = 1
    solid[-1, :, :] = 1
    solid[:, 6:10, 6:10] = 1
    stamp = {"u_x": np.sin(k * (x + 2 * y)) + 0.3 * np.cos(k * z),
             "u_y": 0.5 * np.cos(k * (x - y)) * np.sin(k * z),
             "u_z": np.cos(k * (2 * x + y)) + 0.4 * np.sin(2 * k * y),
             "density": np.ones((n, n, n)), "pressure": np.zeros((n, n, n))}
    write_felbm(d, [stamp, stamp], solid, times=(0, 100))


# Centres in near-solid cells, fluid side of the wall, away from the sub-cubes' mid-planes:
# the z wall, the pillar's faces at x = 6 and 9 and at y = 6, a pillar edge; the last two bulk
NEAR = [(2.3, 3.7, 0.7), (5.3, 7.3, 5.6), (9.7, 8.2, 12.3), (7.6, 5.2, 9.3), (5.3, 5.4, 4.7),
        (2.3, 12.6, 14.3)]
BULK = [(2.3, 12.6, 7.4), (12.2, 3.7, 10.6)]


@pytest.mark.skipif(not os.path.exists(TRACERS), reason="tracers is not built")
def test_lattice_gradient_is_the_velocity_gradient(tmp_path):
    """At points in cells next to solids, and in bulk cells, the dumped J is
    the derivative of the dumped u: differences of u over 1e-5 cells match
    J(i, j) = du_i/dx_j to 1e-6 of the field's gradient. Every row of J taken
    with one component's weights fails the uy and uz rows. The points are
    placed through a checkpoint and dumped one short step on, so the
    differences are taken over the dumped points."""
    from runs import checkpoint_folder, put_points, read_checkpoint
    from dumps import dump_at
    d = tmp_path / "felbm"
    d.mkdir()
    wall_felbm(d)
    delta = 1e-5
    centres = np.array(NEAR + BULK)
    m = len(centres)
    pts = [centres] + [centres + s * delta * e for e in np.eye(3) for s in (1, -1)]
    pts = np.concatenate(pts)
    base = ("mode=felbm init_mode=points_x scheme=RK4 Dm=0 dt=0.001 int_order=2 Nrw=%d Nrw_max=%d "
            "x0=2 y0=2 z0=8 random=false seed=1 stat_intv=1e9" % (len(pts), len(pts)))
    run_app(TRACERS, d / "felbm_params.dat", base, "T=0.001 dump_intv=1e9 checkpoint_intv=0.001")
    put_points(d, pts)
    ids = read_checkpoint(d)["id"][:, 0]
    ck = checkpoint_folder(d)
    run_app(TRACERS, d / "felbm_params.dat", base,
            "T=0.002 dump_intv=0.001 checkpoint_intv=1e9 output_J=true restart_folder=%s" % ck)
    g = dump_at(d, 0.002, raw=True)
    # back in the order put
    slot = np.argsort(g["id"][:, 0])[np.argsort(np.argsort(ids))]
    x = g["points"][slot].reshape(7, m, 3)
    u = g["u"][slot].reshape(7, m, 3)
    J = g["J"][slot].reshape(7, m, 3, 3)[0]
    assert np.abs(x - pts.reshape(7, m, 3)).max() < 0.01
    # du = J dx over the three pairs
    dX = np.stack([x[1 + 2 * k] - x[2 + 2 * k] for k in range(3)], axis=-1)
    dU = np.stack([u[1 + 2 * k] - u[2 + 2 * k] for k in range(3)], axis=-1)
    fd = dU @ np.linalg.inv(dX)
    assert np.abs(fd).max() > 0.1
    err = np.abs(J - fd) / np.abs(fd).max()
    rows = err.max(axis=-1)   # per centre and row
    assert rows.max() < 1e-6, rows


TENSORS = app("tracertensors")


@pytest.mark.skipif(not os.path.exists(TENSORS), reason="tracertensors is not built")
def test_steady_identity_next_to_solids(tmp_path):
    """In a steady field F carries the velocity along: F(T) u(x0) = u(x(T)).
    With tracers all over the lattice of wall_felbm, many next to its walls,
    90% hold it to 2% of |u|; a J next to solids that is not the derivative of
    u leaves a tenth of them off by more than a third."""
    from dumps import all_dumps, deformation_gradient
    d = tmp_path / "felbm"
    d.mkdir()
    wall_felbm(d)
    run_app(TENSORS, d / "felbm_params.dat",
            "mode=felbm init_mode=points_xyz scheme=RK4 Dm=0 dt=0.02 int_order=1 Nrw=1000 "
            "Nrw_max=1000 x0=0 y0=0 z0=0 random=true seed=3 T=2 dump_intv=2 stat_intv=1e9 "
            "checkpoint_intv=1e9")
    D = all_dumps(d)
    F = deformation_gradient(D[2.0])
    u0, u = D[0.0]["u"], D[2.0]["u"]
    moving = np.linalg.norm(u, axis=1) > 1e-3
    e = np.linalg.norm(np.einsum("nij,nj->ni", F, u0) - u, axis=1)[moving] / np.linalg.norm(u, axis=1)[moving]
    assert moving.sum() > 900
    assert np.percentile(e, 90) < 0.02, np.percentile(e, [50, 90, 100])
