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
import subprocess

import numpy as np
import pytest

from paths import app

FELBM = app("filaments_felbmRK4")
INTERPOL = app("interpol")


def three_stamp_felbm(d, shear=False, extra=""):
    """Write a felbm case in d where stamp k holds u_z = k, and with shear also u_y = k x.

    extra is appended to felbm_params.dat."""
    h5py = pytest.importorskip("h5py")
    n = 16
    zero = np.zeros((n, n, n))
    x = np.arange(n, dtype=float)[:, None, None] * np.ones((n, n, n))
    solid = np.zeros((n, n, n), dtype=np.int32)
    solid[0, :, :] = 1
    solid[-1, :, :] = 1
    with h5py.File(d / "output_is_solid.h5", "w") as f:
        f.create_dataset("is_solid", data=solid)
    for k in range(3):
        fields = {"u_x": zero, "u_y": k * x if shear else zero,
                  "u_z": np.full((n, n, n), float(k)),
                  "density": np.ones((n, n, n)), "pressure": zero}
        with h5py.File(d / ("output_%d.h5" % k), "w") as f:
            for name, a in fields.items():
                f.create_dataset(name, data=np.transpose(a, (2, 1, 0)).astype(float))
    (d / "timestamps.dat").write_text("".join("%d\toutput_%d.h5\n" % (k, k) for k in range(3)))
    (d / "felbm_params.dat").write_text(
        "timestamps=timestamps.dat\nis_solid_file=output_is_solid.h5\n" + extra)


@pytest.mark.skipif(not os.path.exists(FELBM), reason="filaments_felbmRK4 is not built")
def test_uz_advances_with_the_timestamp(tmp_path):
    """Past stamp 1 the field is blended between stamps 1 and 2, so particles
    move at u_z > 0.9 over [1, 1.5] and z(T) follows t^2/2. If a component
    kept an older stamp, or the field were held constant between stamps, every
    structured run would advect with a field that lags the data in time."""
    h5py = pytest.importorskip("h5py")
    d = tmp_path / "felbm"
    d.mkdir()
    three_stamp_felbm(d)
    r = subprocess.run([FELBM, str(d / "felbm_params.dat")] +
                       ("Dm=0 dt=0.01 T=2.0 Nrw=2 Nrw_max=100 dump_intv=0.5 stat_intv=1e9 "
                        "checkpoint_intv=1e9 init_mode=pairs_xy int_order=1 ds_init=0.5 "
                        "x0=8 y0=8 z0=8 random=false seed=1").split(),
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr
    f = list(d.rglob("data_from_t*.h5"))
    assert len(f) == 1
    with h5py.File(f[0], "r") as h:
        z = {float(k): np.array(h[k + "/points"])[:, 2].mean() for k in h.keys()}
    assert 1.0 in z and 1.5 in z, sorted(z)
    u_z_after = (z[1.5] - z[1.0]) / 0.5
    assert u_z_after > 0.9, "u_z over [1, 1.5] is %.3f: stamp 1 did not become prev" % u_z_after
    # z(t) = t^2/2; first-order Euler at dt = 0.01 falls short by t*dt/2, well inside 0.05
    t = max(z)
    assert abs((z[t] - 8.0) - t * t / 2) < 0.05, (t, z[t] - 8.0, t * t / 2)


def run_probe(d, t0, int_order):
    return subprocess.run([INTERPOL, str(d / "felbm_params.dat")] +
                          ("mode=felbm Nrw=5000 int_order=%d t0=%g random=false seed=1"
                           % (int_order, t0)).split(),
                          capture_output=True, text=True, timeout=600)


def probe(d, t0, int_order=2):
    """Run interpol on the case in d at time t0; return the probed datasets by name."""
    h5py = pytest.importorskip("h5py")
    r = run_probe(d, t0, int_order)
    assert r.returncode == 0, r.stdout + r.stderr
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
    r = run_probe(d, 0.5, int_order=2)
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
