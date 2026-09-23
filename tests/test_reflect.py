"""Diffusive tracers reflect off walls instead of stopping at them.

With scheme=explicit and Dm > 0, every step is walked from its start and
mirrored at every wall it meets, drift and noise together (the symmetrized
Euler scheme for the reflected SDE, whose density has zero total flux through
the wall). A step that is declined instead keeps
the particle where it was for that step, so a layer of width sqrt(2 Dm dt)
above every wall loses part of its drift. Three walls, one per mechanism:

- a triangle mesh of the unit square, periodic in y, walls at x = 0 and 1
  (the cell walk through mesh facets);
- the analytic Taylor-Couette annulus, whose inner cylinder turns at unit
  speed (a wall with tangential velocity, found through a signed distance);
- a felbm grid with solid node planes at z = 0 and 15 (the walk on the node
  lattice, whose walls are the mid-planes at z = 0.5 and 14.5).

In each, no step is declined, no particle leaves the fluid, and a uniform
density stays uniform. On the annulus the swept angle of particles started
next to the moving cylinder equals the sum of the angular velocities at the
dumped positions, which it would not if steps were lost.
"""

import os
import shutil

import numpy as np
import pytest

from dumps import dump_at, read_stats
from paths import REPO, app
from runs import run_app

TRACERS = app("tracers")

pytestmark = pytest.mark.skipif(not os.path.exists(TRACERS), reason="tracers is not built")


def run(params, args):
    """Run tracers on params with args; assert it succeeded."""
    run_app(TRACERS, params, args, timeout=600)


def restart_from(params, mode, points, args, dt):
    """Run tracers from a checkpoint holding `points` (n x 3); returns the run folder.

    A one-step run writes the checkpoint, at t = 2 dt (a step past the loop),
    and its positions are overwritten; the restart continues from there.
    """
    folder = params.parent
    for old in folder.glob("Tracers*"):
        shutil.rmtree(old)
    base = ["mode=" + mode, "init_mode=points_xyz", "int_order=1", "dt=%g" % dt,
            "Nrw=%d" % len(points), "Nrw_max=%d" % len(points), "x0=0", "y0=0", "z0=0",
            "random=false", "seed=1", "scheme=explicit"]
    run(params, base + ["Dm=0", "T=%g" % dt, "dump_intv=1000", "stat_intv=1000",
                        "checkpoint_intv=%g" % dt])
    [pos] = list(folder.rglob("Checkpoints/positions.pos"))
    np.savetxt(pos, points, fmt="%.17g")
    run(params, base + args + ["checkpoint_intv=1e9", "restart_folder=" + str(pos.parent.parent)])
    return pos.parent.parent


def declined(folder):
    """The n_declined column of the run's statistics, one row per stat time."""
    [f] = [p for p in folder.glob("tdata_from_t*.dat") if p.name != "tdata_from_t0.000000.dat"]
    return read_stats(f)["n_declined"]


def assert_uniform(x, lo, hi, bins):
    """x is spread uniformly over [lo, hi]: every bin within five binomial
    standard deviations of its expected count."""
    counts = np.histogram(x, bins=bins, range=(lo, hi))[0]
    p = 1.0 / bins
    expected = len(x) * p
    sd = np.sqrt(len(x) * p * (1 - p))
    assert np.all(np.abs(counts - expected) < 5 * sd), counts


def lattice(n, lo, hi):
    """n cell-centred points on [lo, hi]."""
    return lo + (np.arange(n) + 0.5) * (hi - lo) / n


def test_triangle_mesh_walls_reflect(mesh_dir, tmp_path):
    """Plane Poiseuille flow u_y = 6x(1-x) on a 10 x 10 mesh, sigma = 0.03,
    a third of a cell: over 400 steps every particle meets a wall many times
    (the channel is crossed by diffusion in about T = 11), none is declined or
    leaves, and the x density stays flat."""
    case = tmp_path / "channel"
    shutil.copytree(mesh_dir("triangle"), case)
    x, y = np.meshgrid(lattice(60, 0, 1), lattice(60, 0, 1), indexing="ij")
    points = np.c_[x.ravel(), y.ravel(), np.zeros(x.size)]
    dt, T = 0.01, 4.0
    Dm = 0.03**2 / (2 * dt)
    folder = restart_from(case / "dolfin_params.dat", "triangle", points,
                          ["Dm=%g" % Dm, "T=%g" % T, "dump_intv=%g" % T, "stat_intv=0.5"], dt)
    assert np.all(declined(folder) == 0), declined(folder)
    xT = dump_at(folder, T)["points"][:, 0]
    assert xT.min() >= 0 and xT.max() <= 1, (xT.min(), xT.max())
    assert_uniform(xT, 0, 1, 20)


TAYLOR_COUETTE = os.path.join(REPO, "data_example", "taylor_couette", "expr_params.dat")


def angular_velocity(g):
    """u_theta / s at each dumped position (the annulus axis is the z axis)."""
    x, u = g["points"], g["u"]
    s2 = x[:, 0]**2 + x[:, 1]**2
    return (x[:, 0] * u[:, 1] - x[:, 1] * u[:, 0]) / s2


def test_moving_cylinder_keeps_its_drag(tmp_path):
    """Particles start within 2 sigma of the inner cylinder (radius 1, turning
    at unit speed), sigma = 0.02. With int_order=1 a step moves by u dt at its
    start, so a particle's swept angle is the sum over steps of u_theta/s dt at
    the dumped positions, up to the noise (zero mean, the angle is harmonic)
    and the angle a wall bounce turns the step through, O(sigma^2). A declined
    step adds its u_theta/s to the sum without moving, so declining near the
    wall makes the swept angle fall short of the sum."""
    case = tmp_path / "tc"
    case.mkdir()
    shutil.copy(TAYLOR_COUETTE, case / "expr_params.dat")
    dt, steps, sigma = 0.01, 20, 0.02
    n = 2000
    rng = np.random.default_rng(1)
    s = 1 + 2 * sigma * (np.arange(n) + 0.5) / n
    phi = rng.uniform(0, 2 * np.pi, n)
    z = rng.uniform(-0.5, 0.5, n)
    points = np.c_[s * np.cos(phi), s * np.sin(phi), z]
    T = (steps + 2) * dt
    folder = restart_from(case / "expr_params.dat", "analytic", points,
                          ["Dm=%g" % (sigma**2 / (2 * dt)), "T=%g" % T,
                           "dump_intv=%g" % dt, "stat_intv=%g" % dt], dt)
    assert np.all(declined(folder) == 0), declined(folder)
    times = dt * np.arange(2, steps + 3)
    dumps = [dump_at(folder, t) for t in times]
    angles = np.array([np.arctan2(g["points"][:, 1], g["points"][:, 0]) for g in dumps])
    swept = np.sum(np.angle(np.exp(1j * np.diff(angles, axis=0))), axis=0)
    expected = np.sum([angular_velocity(g) * dt for g in dumps[:-1]], axis=0)
    radius = np.hypot(dumps[-1]["points"][:, 0], dumps[-1]["points"][:, 1])
    assert radius.min() >= 1 and radius.max() <= 2.5, (radius.min(), radius.max())
    assert abs(swept.mean() / expected.mean() - 1) < 0.025, (swept.mean(), expected.mean())


def test_felbm_solid_planes_reflect(felbm_dir, tmp_path):
    """The felbm case's fluid nodes are z = 1..14, so the walls are the
    mid-planes z = 0.5 and 14.5, and the shear u_y(x) runs along them.
    sigma = 0.3 node spacings over 800 steps crosses the channel about twice;
    no step is declined, every particle ends in the fluid and the z density
    stays flat."""
    case = tmp_path / "felbm"
    shutil.copytree(felbm_dir, case)
    x, z = np.meshgrid(lattice(40, 0, 16), lattice(40, 0.5, 14.5), indexing="ij")
    points = np.c_[x.ravel(), np.full(x.size, 8.0), z.ravel()]
    dt, T = 0.1, 80.0
    Dm = 0.3**2 / (2 * dt)
    folder = restart_from(case / "felbm_params.dat", "felbm", points,
                          ["Dm=%g" % Dm, "T=%g" % T, "dump_intv=%g" % T, "stat_intv=10"], dt)
    assert np.all(declined(folder) == 0), declined(folder)
    zT = dump_at(folder, T)["points"][:, 2]
    assert zT.min() >= 0.5 and zT.max() <= 14.5, (zT.min(), zT.max())
    assert_uniform(zT, 0.5, 14.5, 14)
