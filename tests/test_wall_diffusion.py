"""Diffusing tracers on the wall meshes of test_wall_p2.py (scheme=explicit, Dm > 0).

- the near-wall rule and the reflecting walk together: wall_p2=edge makes
  u.n quadratic above the walls so that advected tracers do not reach them
  and stop, and the walk mirrors a diffusive step off the walls instead of
  declining it; in the four cylinders (unit mean speed, a noise step about a
  thirteenth of a cell) the fraction at rest and the density next to the
  cylinders, with the rule and without;
- the walk across periodic faces in a run, in the sphere of the periodic
  cube, where it was tested only on single moves (test_interpol_core.cpp).
"""

import shutil
import subprocess

import numpy as np
import pytest

from dumps import all_dumps, dump_at, read_stats
from runs import checkpoint_folder, put_points
from test_wall_p2 import (DT, OBSTACLE_R, OBSTACLES, SPHERE_C, SPHERE_R, TRACERS, H,  # noqa: F401
                          ball, inside, obstacles, sphere_mesh)

DM = 1e-4


def diffusive_run(case, mode, points, T, dt=DT, Dm=DM, extra=()):
    """tracers from `points` (through a checkpoint) with Dm, explicit steps and
    wall_p2=mode to T; returns the run folder."""
    folder = case(mode)
    for old in folder.glob("Tracers*"):
        shutil.rmtree(old)
    dim = points.shape[1]
    seed = list(case.seed_point) + [0.0] * (3 - len(case.seed_point))
    base = ["mode=" + case.mode, "init_mode=points_" + "xyz"[:dim], "int_order=1", "dt=%g" % dt,
            "Nrw=%d" % len(points), "Nrw_max=%d" % len(points), "x0=%g" % seed[0], "y0=%g" % seed[1],
            "z0=%g" % seed[2], "random=false", "seed=1", "scheme=explicit", "Dm=%g" % Dm]

    def go(args):
        r = subprocess.run([TRACERS, str(folder / "dolfin_params.dat")] + base + args,
                           capture_output=True, text=True, timeout=1200)
        assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]

    go(["T=%g" % dt, "dump_intv=1000", "stat_intv=1000", "checkpoint_intv=%g" % dt])
    put_points(folder, points)
    ck = checkpoint_folder(folder)
    go(["T=%g" % T, "checkpoint_intv=1e9", "restart_folder=" + str(ck)]
       + list(extra or ["dump_intv=%g" % T, "stat_intv=1000"]))
    return ck


def near_wall(points, delta):
    """Whether each point is within delta of a cylinder (periodic images included)."""
    d = np.full(len(points), np.inf)
    for c in OBSTACLES:
        for s in ((a, b) for a in (-1, 0, 1) for b in (-1, 0, 1)):
            d = np.minimum(d, np.linalg.norm(points - c - s, axis=1))
    return d - OBSTACLE_R < delta


def measure(obstacles, mode, T, dt=DT):
    """(fraction at rest, tracer density within 0.1 h of the cylinders over the fluid's there)."""
    x = (np.arange(50) + 0.5) / 50
    P = np.stack(np.meshgrid(x, x, indexing="ij"), axis=-1).reshape(-1, 2)
    P = P[inside(obstacles.geometry, P)]
    g = dump_at(diffusive_run(obstacles, mode, P, T, dt), T)
    rest = float(np.mean(np.linalg.norm(g["u"][:, :2], axis=1) < 1e-3))
    xy = np.mod(g["points"][:, :2], 1.0)
    delta = 0.1 * H
    area = len(OBSTACLES) * np.pi * ((OBSTACLE_R + delta) ** 2 - OBSTACLE_R ** 2)
    fluid = 1.0 - len(OBSTACLES) * np.pi * OBSTACLE_R ** 2
    return rest, float(np.mean(near_wall(xy, delta)) / (area / fluid))


@pytest.mark.slow
def test_the_rule_still_evens_the_density_of_diffusing_tracers_at_walls(obstacles):
    """Diffusion frees tracers from the walls either way: next to none is at
    rest at T = 200, with the rule or without. The rule still matters for
    where they are: within a tenth of a cell of the cylinders the density
    without it is more than twice that with it and more than four times the
    uniform one, while with it the excess stays below three. The linear u.n
    of the P1 field drives tracers into the wall layer faster than the noise
    spreads them out; the quadratic one does not."""
    rest, dens = {}, {}
    for mode in ("edge", "none"):
        rest[mode], dens[mode] = measure(obstacles, mode, 200.0)
    assert max(rest.values()) < 0.01, rest
    assert dens["none"] > 4 and dens["edge"] < 3 and dens["none"] > 2 * dens["edge"], dens


def test_the_reflecting_walk_crosses_periodic_faces_in_a_run(ball):
    """The sphere in the periodic unit cube (tets of size 0.1), diffusing
    tracers with a noise step of a third of a cell to T = 4: the walk, which
    shifts a point by the box length when it leaves through a periodic face,
    carries them across the faces hundreds of times (a jump of more than half
    the box between dumps 0.1 apart), declines no step, and ends every one
    in the fluid, outside the sphere up to the facets' sagitta."""
    x = (np.arange(10) + 0.5) / 10
    P = np.stack(np.meshgrid(x, x, x, indexing="ij"), axis=-1).reshape(-1, 3)
    P = P[np.linalg.norm(P - SPHERE_C, axis=1) > SPHERE_R + 0.02]
    dt, T = 0.01, 4.0
    folder = diffusive_run(ball, "edge", P, T, dt=dt, Dm=0.03 ** 2 / (2 * dt),
                           extra=["dump_intv=0.1", "stat_intv=0.5"])
    [f] = [p for p in folder.glob("tdata_from_t*.dat") if p.name != "tdata_from_t0.000000.dat"]
    assert np.all(read_stats(f)["n_declined"] == 0)
    D = all_dumps(folder)
    t = sorted(D)
    assert len(t) > 30
    X = np.stack([D[k]["points"] for k in t])
    assert len(X[-1]) == len(P)
    jumps = np.abs(np.diff(X, axis=0)) > 0.5
    assert jumps.sum() > 200, jumps.sum()
    xT = np.mod(X[-1], 1.0)
    assert np.linalg.norm(xT - SPHERE_C, axis=1).min() > SPHERE_R - 0.01
