"""A point outside the fluid has zero velocity, in every step and in the dumps.

The flow is the Stokes flow past a sphere (data_example/stokes_sphere), whose
interior is outside the fluid. On the symmetry axis x = y = 0 the velocity is
exactly axial, u_z = u_inf f_r(R/|z - z0|) with
f_r(eta) = 1 + eta^3/2 - 3 eta/2, so a particle on the axis stays on it and
its step is a one-dimensional closed form that can be recomputed here with the
same floating-point operations. Upstream of the sphere the flow runs into it;
with a long step, an RK4 stage point lands inside the sphere, where the formula
still returns a velocity that is not the flow's.

Two things are pinned:

- an RK4 stage whose point is outside contributes zero velocity, and a step
  that ends outside is declined; positions after one step match the reference
  for every transport element, since the position update does not depend on
  what else a particle carries;
- a particle placed outside the fluid (by editing a checkpoint) reports zero
  velocity, pressure and density in the dumps and does not move, under both
  schemes.
"""

import math
import os
import shutil
import subprocess

import h5py
import numpy as np
import pytest

from paths import REPO, app

SPHERE = os.path.join(REPO, "data_example", "stokes_sphere", "expr_params.dat")
APPS = ["tracers", "tracervectors", "tracertensors"]
DT = 8.0

# 100 particles on the z axis through the sphere's centre, noise off
BASE = ("mode=analytic init_mode=points_z x0=0 y0=0 z0=0 Nrw=100 Nrw_max=100 "
        "Dm=0 int_order=1 stat_intv=1e9 random=false seed=1 outside=ignore").split()

needs_apps = pytest.mark.skipif(not all(os.path.exists(app(a)) for a in APPS),
                                reason="the tracer apps are not built")


def params():
    """The sphere example's parameters, as floats where they parse."""
    out = {}
    with open(SPHERE) as f:
        for line in f:
            if "=" in line:
                k, v = (s.strip() for s in line.split("=", 1))
                try:
                    out[k] = float(v)
                except ValueError:
                    out[k] = v
    return out


def run(name, case_dir, extra):
    """Run app `name` on the case's parameter file with BASE and extra; assert it succeeded."""
    keys = {a.split("=")[0] for a in extra}
    argv = [a for a in BASE if a.split("=")[0] not in keys] + list(extra)
    r = subprocess.run([app(name), str(case_dir / "expr_params.dat")] + argv,
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr


def case(tmp_path, name):
    """A new case folder holding the sphere example."""
    d = tmp_path / name
    d.mkdir()
    shutil.copy(SPHERE, d / "expr_params.dat")
    return d


def dumps(case_dir):
    """{t: {dataset: array in id order}} for every dump under case_dir."""
    out = {}
    for f in sorted(case_dir.rglob("data_from_t*.h5")):
        with h5py.File(f, "r") as h:
            for key in h:
                g = h[key]
                ids = np.array(g["id"])[:, 0]
                order = np.argsort(ids, kind="stable")
                assert np.array_equal(ids[order], np.arange(len(ids)))
                out[float(key)] = {k: np.array(g[k])[order] for k in g if k != "id"}
    return out


class Axis:
    """The Stokes sphere on its axis, with the operations of Expr_StokesSphere
    in the same order, so that the reference is exact."""

    def __init__(self, p):
        self.R, self.u_inf, self.zc = p["R"], p["u_inf"], p["z0"]

    def inside(self, z):
        rz = z - self.zc
        return rz * rz >= self.R * self.R

    def u(self, z):
        rz = z - self.zc
        r2 = rz * rz
        eta = self.R / math.sqrt(r2)
        eta3 = math.pow(eta, 3)
        f_r = 1.0 + 1. / 2. * eta3 - 3. / 2. * eta
        f_theta = -1.0 + 1. / 4. * eta3 + 3. / 4. * eta
        return rz * rz / r2 * (f_r + f_theta) * self.u_inf - f_theta * self.u_inf

    def rk4(self, z, dt, gated=True):
        """One RK4 step from z; a stage outside contributes zero unless gated is
        False; a step ending outside is declined. Returns (z, an outside stage)."""
        def k(zz):
            return self.u(zz) if (self.inside(zz) or not gated) else 0.0
        k1 = k(z)
        k2 = k(z + k1 * dt / 2)
        k3 = k(z + k2 * dt / 2)
        k4 = k(z + k3 * dt)
        crossed = not all(self.inside(zz) for zz in
                          (z, z + k1 * dt / 2, z + k2 * dt / 2, z + k3 * dt))
        dz = (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6
        return (z + dz if self.inside(z + dz) else z), crossed


@needs_apps
@pytest.mark.parametrize("name", APPS)
def test_an_rk4_stage_outside_the_fluid_contributes_zero_velocity(tmp_path, name):
    """After one long RK4 step every particle is where a step that zeroes the
    velocity at its outside stage points puts it. Some particles must cross
    into the sphere at a stage and still end in the fluid, and their positions
    must differ from those of a step that uses the formula inside the sphere;
    otherwise the flow would not test the rule."""
    d = case(tmp_path, name)
    run(name, d, ["scheme=RK4", "dt=%g" % DT, "T=%g" % DT, "dump_intv=%g" % DT,
                  "checkpoint_intv=1e9"])
    D = dumps(d)
    p0, p1 = D[0.0]["points"], D[DT]["points"]
    assert np.abs(p0[:, :2]).max() == 0 and np.abs(p1[:, :2]).max() == 0

    axis = Axis(params())
    expected, discriminating = [], 0
    for z in p0[:, 2]:
        zg, crossed = axis.rk4(z, DT)
        zu, _ = axis.rk4(z, DT, gated=False)
        expected.append(zg)
        if crossed and axis.inside(zg) and abs(zg - zu) > 1e-6:
            discriminating += 1
    assert discriminating > 0, "no particle crossed the sphere at a stage"
    assert np.allclose(p1[:, 2], expected, rtol=0, atol=1e-12)


@needs_apps
@pytest.mark.parametrize("scheme", ["RK4", "explicit"])
def test_a_particle_outside_the_fluid_has_zero_fields_and_stays(tmp_path, scheme):
    """A particle moved into the sphere through a checkpoint dumps u = 0,
    p = 0 and rho = 0 and keeps its position, at every dump of the resumed run:
    nothing is evaluated outside the fluid, so it has no velocity to move with,
    and its steps are declined."""
    d = case(tmp_path, scheme)
    run("tracers", d, ["scheme=%s" % scheme, "dt=%g" % DT, "T=%g" % DT,
                       "dump_intv=%g" % DT, "checkpoint_intv=%g" % DT])
    [pos] = list(d.rglob("Checkpoints/positions.pos"))
    ids = np.loadtxt(pos.parent / "id.list", dtype=int)
    lines = pos.read_text().splitlines()
    zc = params()["z0"]
    lines[0] = "0 0 %.17g" % (zc - 0.5)                  # inside the sphere
    pos.write_text("\n".join(lines) + "\n")
    moved = ids[0]

    # the checkpoint is at t = 2 DT; the resumed run dumps from there on
    run("tracers", d, ["scheme=%s" % scheme, "dt=%g" % DT, "T=%g" % (4 * DT),
                       "dump_intv=%g" % DT, "checkpoint_intv=1e9",
                       "restart_folder=" + str(pos.parent.parent)])
    later = {t: g for t, g in dumps(d).items() if t > 1.5 * DT}
    assert later, "the resumed run dumped nothing"
    for t, g in sorted(later.items()):
        assert np.array_equal(g["points"][moved], [0., 0., zc - 0.5]), t
        assert np.all(g["u"][moved] == 0), t
        assert np.all(g["p"][moved] == 0) and np.all(g["rho"][moved] == 0), t
