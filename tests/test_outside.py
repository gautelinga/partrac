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

- an RK4 step with a stage point or its end outside is taken again in 2, 4,
  then 8 substeps whose every stage is inside, and declined if none is; a
  failed stage used to contribute zero velocity, which froze a tracer whose
  second stage left the fluid on the start point for good. Positions after
  one step match the reference for every transport element, since the
  position update does not depend on what else a particle carries;
- a particle placed outside the fluid (by editing a checkpoint) reports zero
  velocity, pressure and density in the dumps and does not move, under both
  schemes;
- with outside=reinject, a particle that no offset along init_mode's
  directions can bring back into the fluid is refused after a bounded number
  of draws, rather than drawing for ever.
"""

import math
import os

import numpy as np
import pytest

from dumps import all_dumps
from paths import REPO, app
from runs import checkpoint_folder, copy_example, put_points, read_checkpoint, run_app

SPHERE = os.path.join(REPO, "data_example", "stokes_sphere", "expr_params.dat")
POISEUILLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")
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
    run_app(app(name), case_dir / "expr_params.dat", BASE, extra, timeout=600)


def case(tmp_path, name):
    """A new case folder holding the sphere example."""
    return copy_example(SPHERE, tmp_path / name).parent


def dumps(case_dir):
    """{t: {dataset: array in id order}} for every dump under case_dir, where every id is there once."""
    out = {}
    for t, g in all_dumps(case_dir, raw=True).items():
        ids = g.pop("id")[:, 0]
        order = np.argsort(ids, kind="stable")
        assert np.array_equal(ids[order], np.arange(len(ids)))
        out[t] = {k: a[order] for k, a in g.items()}
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

    def rk4_once(self, z, dt):
        """One RK4 step from z with the formula everywhere: (end, every stage and the end inside)."""
        k1 = self.u(z)
        k2 = self.u(z + k1 * dt / 2)
        k3 = self.u(z + k2 * dt / 2)
        k4 = self.u(z + k3 * dt)
        dz = (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6
        inside = all(self.inside(zz) for zz in (z, z + k1 * dt / 2, z + k2 * dt / 2, z + k3 * dt, z + dz))
        return z + dz, inside

    def rk4(self, z, dt):
        """One step from z as the apps take it: whole, else in 2, 4, then 8
        substeps that stay inside, else declined. Returns (z, substeps used)."""
        end, ok = self.rk4_once(z, dt)
        if ok:
            return end, 1
        for m in (2, 4, 8):
            zz, ok = z, True
            for _ in range(m):
                zz, ok = self.rk4_once(zz, dt / m)
                if not ok:
                    break
            if ok:
                return zz, m
        return z, 0


@needs_apps
@pytest.mark.parametrize("name", APPS)
def test_an_rk4_step_with_a_stage_outside_the_fluid_is_taken_in_substeps(tmp_path, name):
    """After one long RK4 step every particle is where the reference puts it:
    the whole step where every stage is in the fluid, else the first of 2, 4
    or 8 substeps that stays in it, else its start. Some particles must need
    substeps and move, and some must need more than two; otherwise the flow
    would not test the rule."""
    d = case(tmp_path, name)
    run(name, d, ["scheme=RK4", "dt=%g" % DT, "T=%g" % DT, "dump_intv=%g" % DT,
                  "checkpoint_intv=1e9"])
    D = dumps(d)
    p0, p1 = D[0.0]["points"], D[DT]["points"]
    assert np.abs(p0[:, :2]).max() == 0 and np.abs(p1[:, :2]).max() == 0

    axis = Axis(params())
    expected, used = [], []
    for z in p0[:, 2]:
        zz, m = axis.rk4(z, DT)
        expected.append(zz)
        used.append(m)
    used = np.array(used)
    assert np.sum(used > 1) > 0 and np.sum(used > 2) > 0, np.bincount(used)
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
    ck = read_checkpoint(d)
    zc = params()["z0"]
    ck["points"][0] = [0., 0., zc - 0.5]                # inside the sphere
    put_points(d, ck["points"])
    moved = int(ck["id"][0, 0])

    # the checkpoint is at t = 2 DT; the resumed run dumps from there on
    run("tracers", d, ["scheme=%s" % scheme, "dt=%g" % DT, "T=%g" % (4 * DT),
                       "dump_intv=%g" % DT, "checkpoint_intv=1e9",
                       "restart_folder=" + str(checkpoint_folder(d))])
    later = {t: g for t, g in dumps(d).items() if t > 1.5 * DT}
    assert later, "the resumed run dumped nothing"
    for t, g in sorted(later.items()):
        assert np.array_equal(g["points"][moved], [0., 0., zc - 0.5]), t
        assert np.all(g["u"][moved] == 0), t
        assert np.all(g["p"][moved] == 0) and np.all(g["rho"][moved] == 0), t


@needs_apps
def test_a_particle_reinjection_cannot_reach_is_refused(tmp_path):
    """Plane Poiseuille flow, fluid where |x| <= 1: a particle moved to
    x = 1.5, in the wall, with init_mode=points_z. Reinjection offsets it
    along z only, so no draw is inside; the resumed run stops with a message
    naming the directions instead of drawing for ever. The timeout catches
    the hang."""
    d = copy_example(POISEUILLE, tmp_path / "reinject").parent
    run("tracers", d, ["scheme=RK4", "dt=0.1", "T=0.1", "dump_intv=0.1", "checkpoint_intv=0.1"])
    x = read_checkpoint(d)["points"]
    x[0] = [1.5, 0., 0.]
    put_points(d, x)
    r = run_app(app("tracers"), d / "expr_params.dat", BASE,
                ["scheme=RK4", "dt=0.1", "T=0.3", "dump_intv=0.1", "checkpoint_intv=1e9",
                 "outside=reinject", "restart_folder=" + str(checkpoint_folder(d))],
                check=False, timeout=120)
    assert r.returncode == 2, r.stdout[-1000:] + r.stderr[-1000:]
    assert "outside=reinject: no position inside the domain along z in 1000000 draws" in r.stderr, r.stderr
