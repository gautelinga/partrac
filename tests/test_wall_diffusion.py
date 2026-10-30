"""Diffusing tracers on the wall meshes of test_wall_edge_rule.py (scheme=explicit, Dm > 0).

- the near-wall rule and the reflecting walk together: wall_upgrade=true makes
  u.n quadratic above the walls so that advected tracers do not reach them
  and stop, and the walk mirrors a diffusive step off the walls instead of
  declining it; in the four cylinders (unit mean speed, a noise step about a
  thirteenth of a cell) the fraction at rest and the density next to the
  cylinders, with the rule and without;
- the walk across periodic faces in a run, in the sphere of the periodic
  cube, where it was tested only on single moves (test_interpol_core.cpp).
"""

import numpy as np
import pytest

from dumps import all_dumps, dump_at, read_stats
from wall_cases import (DT, OBSTACLE_R, OBSTACLES, SPHERE_C, SPHERE_R, H, lattice, restart_from,  # noqa: F401
                        ball, obstacles, sphere_mesh)

DM = 1e-4


def diffusive_run(case, mode, points, T, dt=DT, Dm=DM, extra=()):
    """tracers from `points` (through a checkpoint) with Dm, explicit steps and
    wall_upgrade=mode to T; returns the run folder."""
    return restart_from(case, mode, points, ["T=%g" % T] + list(extra or ["dump_intv=%g" % T, "stat_intv=1000"]),
                        scheme="explicit", dt=dt, Dm=Dm)


def near_wall(points, delta):
    """Whether each point is within delta of a cylinder (periodic images included)."""
    d = np.full(len(points), np.inf)
    for c in OBSTACLES:
        for s in ((a, b) for a in (-1, 0, 1) for b in (-1, 0, 1)):
            d = np.minimum(d, np.linalg.norm(points - c - s, axis=1))
    return d - OBSTACLE_R < delta


def measure(obstacles, mode, T, dt=DT):
    """(fraction at rest, tracer density within 0.1 h of the cylinders over the fluid's there)."""
    g = dump_at(diffusive_run(obstacles, mode, lattice(50, 2, obstacles.geometry), T, dt), T)
    rest = float(np.mean(np.linalg.norm(g["u"][:, :2], axis=1) < 1e-3))
    xy = np.mod(g["points"][:, :2], 1.0)
    delta = 0.1 * H
    area = len(OBSTACLES) * np.pi * ((OBSTACLE_R + delta) ** 2 - OBSTACLE_R ** 2)
    fluid = 1.0 - len(OBSTACLES) * np.pi * OBSTACLE_R ** 2
    return rest, float(np.mean(near_wall(xy, delta)) / (area / fluid))


@pytest.mark.slow
def test_the_rule_still_evens_the_density_of_diffusing_tracers_at_walls(obstacles):
    """Diffusion frees tracers from the walls either way: next to none is at
    rest at T = 50, with the rule or without. The rule still matters for
    where they are: within a tenth of a cell of the cylinders the density
    without it is more than twice that with it and more than four times the
    uniform one, while with it the excess stays below three. The linear u.n
    of the P1 field drives tracers into the wall layer faster than the noise
    spreads them out; the quadratic one does not."""
    rest, dens = {}, {}
    for mode in ("true", "false"):
        rest[mode], dens[mode] = measure(obstacles, mode, 50.0)
    assert max(rest.values()) < 0.01, rest
    assert dens["false"] > 4 and dens["true"] < 3 and dens["false"] > 2 * dens["true"], dens


def test_the_reflecting_walk_crosses_periodic_faces_in_a_run(ball):
    """The sphere in the periodic unit cube (tets of size 0.1), diffusing
    tracers with a noise step of a third of a cell to T = 4: the walk, which
    shifts a point by the box length when it leaves through a periodic face,
    carries them across the faces hundreds of times (a jump of more than half
    the box between dumps 0.1 apart), declines no step, and ends every one
    in the fluid, outside the sphere up to the facets' sagitta."""
    P = lattice(10, 3)
    P = P[np.linalg.norm(P - SPHERE_C, axis=1) > SPHERE_R + 0.02]
    dt, T = 0.01, 4.0
    folder = diffusive_run(ball, "true", P, T, dt=dt, Dm=0.03 ** 2 / (2 * dt),
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
