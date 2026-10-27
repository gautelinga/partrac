"""The quadratic velocity XDMFTriangleInterpol and XDMFTetInterpol build near
no-slip walls.

XDMF velocities are P1. In every cell with a vertex on a wall, wall_upgrade=true
makes the field P2, each edge midpoint a fixed linear combination of the
velocity at the edge's fluid end, so the field is continuous, and the velocity
normal to a wall facet grows as the square of the distance to it, as the
no-slip asymptotics has it. An edge from a wall vertex above two or more wall
facets at an angle takes their mean normal, so near those facets the normal
velocity keeps a small linear part. wall_upgrade=false leaves the field P1, whose
normal velocity grows linearly: tracers approaching a wall reach it and stop.

Four cases, the first, second and fourth Stokes flow driven by a uniform body
force in a periodic unit box:

- one cylinder (radius 0.25) with the force along x: the P2-P1 solution
  written as vertex values the way a solver's P1 output is, on a gmsh mesh
  with periodic node pairs on opposite sides and cell size 0.04;
- four cylinders (radius 0.12) with the force 20 degrees off the x axis, so
  no streamline is symmetric and what a tracer loses passing one obstacle is
  never given back at the next;
- a synthetic field in a channel periodic in x whose bottom wall has a kink at
  x = 1/2, meshed so that the cells on both facets at the kink share their
  fluid vertex, vanishing on the walls up to round-off;
- in 3D, a sphere (radius 0.25) in the periodic unit cube, cell size 0.1,
  under a force whose cross-stream part varies in time, so tracers pass the
  sphere along ever different paths; and on the same mesh a synthetic field
  at rest on the sphere, for the local properties of the tet rule.

A velocity is probed by restarting tracers from a checkpoint whose positions
are the probe points: the restarted run dumps the velocity before it steps.
Mean velocities come from the statistics of a run seeded on a lattice, since
random seeding leaves a sampling error as large as the effects at issue.
"""

import os

import numpy as np
import pytest

from dumps import all_dumps, dump_at, read_stats
from runs import cells_counts
from wall_cases import DT, OBSTACLE_FORCE, TRACERS, inside, probe, restart_from
from wall_cases import ball, cylinder, kinked, obstacles, sphere, sphere_mesh  # noqa: F401  (fixtures)


pytestmark = pytest.mark.skipif(not os.path.exists(TRACERS), reason="tracers is not built")


@pytest.fixture(params=["cylinder", "kinked"])
def case(request):
    """Each case in turn."""
    return request.getfixturevalue(request.param)


def wall_facets(case):
    """(a, b, n, apex) for each wall facet: its end points, the unit normal into
    the fluid and the index of the cell's third vertex."""
    X, cells = case.geometry
    wall = case.wall(X)
    out = []
    for c in cells:
        for i, j in ((0, 1), (0, 2), (1, 2)):
            if wall[c[i]] and wall[c[j]]:
                a, b, v = X[c[i]], X[c[j]], X[c[3 - i - j]]
                t = b - a
                n = np.array([-t[1], t[0]]) / np.linalg.norm(t)
                out.append((a, b, n if (v - a) @ n > 0 else -n, c[3 - i - j]))
    return out


def wall_edges(case):
    """The edges of cells with a wall vertex, less the wall facets, as coordinate pairs."""
    X, cells = case.geometry
    wall = case.wall(X)
    edges = set()
    for c in cells:
        if not wall[c].any():
            continue
        for i, j in ((0, 1), (0, 2), (1, 2)):
            a, b = sorted((c[i], c[j]))
            if not (wall[a] and wall[b]):
                edges.add((a, b))
    return [(X[a], X[b]) for a, b in sorted(edges)]


def test_the_field_is_continuous_across_edges_near_the_wall(case):
    """Across every edge of a cell with a wall vertex, points a hair apart on
    either side see the same velocity with wall_upgrade=true. The quadratic field
    must differ from the P1 one at those points, or the check would pass on
    an untouched field."""
    eps = 1e-9
    pts = []
    for a, b in wall_edges(case):
        t = b - a
        n = np.array([-t[1], t[0]]) / np.linalg.norm(t)
        for s in (0.25, 0.5, 0.75):
            x = a + s * t
            pts += [x + eps * n, x - eps * n]
    pts = np.array(pts)
    u = probe(case, "true", pts)
    assert np.abs(u[0::2] - u[1::2]).max() < 1e-6
    assert np.abs(u - probe(case, "false", pts)).max() > 1e-3


def test_the_kinked_channel_has_a_kink_between_cells_sharing_a_vertex(kinked):
    """The kinked case is meant to reach the mean-normal rule for an edge
    between two wall facets at an angle; it does only if two such facets'
    cells share their third vertex."""
    facets = wall_facets(kinked)
    kinks = [(f, g) for k, f in enumerate(facets) for g in facets[k + 1:]
             if f[3] == g[3] and abs(f[2] @ g[2]) < 1 - 1e-6]
    assert kinks


def shared_side_edge(facets):
    """For each (a, b, n, apex) facet, whether an edge from one of its ends to
    its apex also rises from another facet at an angle."""
    ends = {}
    for a, b, n, apex in facets:
        for x in (a, b):
            ends.setdefault((tuple(x), apex), []).append(n)
    return [any(len(ends[(tuple(x), apex)]) > 1 and abs(ends[(tuple(x), apex)][0] @ ends[(tuple(x), apex)][1]) < 1 - 1e-12
                for x in (a, b)) for a, b, n, apex in facets]


def test_facet_normal_velocity_grows_as_the_squared_distance(case):
    """Along the inward normal of a wall facet, u.n / delta^2 is the same at
    delta = 1e-5, 1e-4 and 1e-3 with wall_upgrade=true: the velocity normal to the
    facet is exactly quadratic there. A P1 field has u.n ~ delta, so the ratio
    grows tenfold per step. Facets next to an edge above another facet at an
    angle (only at the kink) keep a linear part: at 1e-5, less than a tenth of
    the velocity along the facet, and at each point less than half the P1
    value."""
    facets = wall_facets(case)
    shared = np.repeat(shared_side_edge(facets), 3)
    deltas = np.array([1e-5, 1e-4, 1e-3])
    pts, normals, dist = [], [], []
    for a, b, n, _ in facets:
        for s in (0.2, 0.5, 0.8):
            for d in deltas:
                pts.append(a + s * (b - a) + d * n)
                normals.append(n)
                dist.append(d)
    pts, normals, dist = np.array(pts), np.array(normals), np.array(dist)
    leak = {}
    for mode in ("true", "false"):
        u = probe(case, mode, pts)
        un = np.sum(u * normals, axis=1)
        ratio = (un / dist**2).reshape(-1, len(deltas))
        spread = np.abs(ratio - ratio[:, -1:]).max(axis=1)
        if mode == "true":
            exact = ~shared
            assert spread[exact].max() <= 1e-3 * np.abs(ratio[exact]).max(), (mode, spread[exact].max())
        else:
            assert spread.max() > 10 * np.abs(ratio[:, -1]).max(), (mode, spread.max())
        ut = np.linalg.norm(u - un[:, None] * normals, axis=1)
        leak[mode] = (np.abs(un) / ut)[0::len(deltas)]
    if shared.any():
        assert leak["true"][shared].max() < 0.1, leak["true"][shared].max()
        assert (leak["true"][shared] < 0.5 * leak["false"][shared]).all(), (leak["true"][shared], leak["false"][shared])


def lattice_run(case, mode, n, T):
    """Tracers on a cell-centred lattice with n points a side, those in the
    fluid, run to T with wall_upgrade=mode: (times, mean velocity over time, the
    fraction of tracers at rest at T)."""
    dim = case.geometry[0].shape[1]
    x = (np.arange(n) + 0.5) / n
    points = np.stack(np.meshgrid(*[x] * dim, indexing="ij"), axis=-1).reshape(-1, dim)
    points = points[inside(case.geometry, points)]
    folder = restart_from(case, mode, points,
                          ["T=%g" % T, "dump_intv=%g" % T, "stat_intv=%g" % (2 * DT)])
    [f] = [p for p in folder.glob("tdata_from_t*.dat") if p.name != "tdata_from_t0.000000.dat"]
    d = np.loadtxt(f)
    u = dump_at(folder, T)["u"]
    return d[:, 0], d[:, 7:10], np.mean(np.linalg.norm(u, axis=1) < 1e-3)


@pytest.mark.slow
def test_tracers_keep_the_eulerian_mean_velocity_between_obstacles(obstacles):
    """Tracers seeded uniformly in an incompressible flow with impermeable
    walls stay uniform, so their mean velocity stays at the Eulerian mean.
    In the P1 field each pass of an obstacle moves a tracer a little towards
    its upstream face, where it stops: by t = 200 most tracers are at rest
    and the mean has collapsed. With wall_upgrade=true the mean stays within a few
    percent (its sampling noise is about 2% with 2000 tracers) and next to
    nothing stops."""
    ratio, stuck = {}, {}
    for mode in ("true", "false"):
        _, u, stuck[mode] = lattice_run(obstacles, mode, 50, T=200.0)
        along = u[:, :2] @ OBSTACLE_FORCE
        ratio[mode] = along / along[0]
    assert stuck["false"] > 0.5 and ratio["false"][-1] < 0.3, (stuck, ratio["false"][-1])
    assert stuck["true"] < 0.01, stuck
    late = ratio["true"][len(ratio["true"]) // 2:]
    assert late.mean() > 0.9 and late[-1] > 0.9, (late.mean(), late[-1])


def tet_faces(case):
    """{face (sorted vertex indices): the vertices opposite it} over all tets."""
    X, cells = case.geometry
    faces = {}
    for c in cells:
        for k in range(4):
            faces.setdefault(tuple(sorted(np.delete(c, k))), []).append(c[k])
    return faces


def test_the_tet_field_is_continuous_across_faces_near_the_wall(ball):
    """Across every face of a tet with a wall vertex, points a hair apart on
    either side see the same velocity with wall_upgrade=true, and the quadratic
    field differs from the P1 one there."""
    X, _ = ball.geometry
    wall = ball.wall(X)
    pts = []
    for f, opposite in tet_faces(ball).items():
        if len(opposite) != 2 or not wall[list(f)].any() or wall[list(f)].all():
            continue
        a, b, c = X[list(f)]
        n = np.cross(b - a, c - a)
        n /= np.linalg.norm(n)
        for w in ((1 / 3, 1 / 3, 1 / 3), (0.5, 0.25, 0.25), (0.25, 0.5, 0.25)):
            x = w[0] * a + w[1] * b + w[2] * c
            pts += [x + 1e-9 * n, x - 1e-9 * n]
    pts = np.array(pts)
    u = probe(ball, "true", pts)
    assert np.abs(u[0::2] - u[1::2]).max() < 1e-6
    assert np.abs(u - probe(ball, "false", pts)).max() > 1e-2


def test_tet_facet_normal_velocity_is_nearly_quadratic(ball):
    """Just above the sphere's facets (delta = 1e-5) the velocity through a
    facet is a small fraction of the velocity along it with wall_upgrade=true:
    below 1e-3, of order delta as the quadratic law has it, where the facet's
    cell has no side edge above another facet, and elsewhere, in median and at
    the 95th percentile, an eighth or less of what it is in the P1 field. Tets
    with all four vertices on the sphere carry no flow and are left out."""
    X, _ = ball.geometry
    wall = ball.wall(X)
    facets = [(f, o[0]) for f, o in tet_faces(ball).items()
              if len(o) == 1 and wall[list(f)].all() and not wall[o[0]]]
    apexes = {}
    for f, v in facets:
        for w in f:
            apexes.setdefault((w, v), []).append(f)
    pts, normals, single = [], [], []
    for f, v in facets:
        a, b, c = X[list(f)]
        n = np.cross(b - a, c - a)
        n /= np.linalg.norm(n)
        if (X[v] - a) @ n < 0:
            n = -n
        alone = all(len(apexes[(w, v)]) == 1 for w in f)
        for w in ((1 / 3, 1 / 3, 1 / 3), (0.6, 0.2, 0.2), (0.2, 0.6, 0.2), (0.2, 0.2, 0.6)):
            pts.append(w[0] * a + w[1] * b + w[2] * c + 1e-5 * n)
            normals.append(n)
            single.append(alone)
    pts, normals, single = np.array(pts), np.array(normals), np.array(single)
    assert single.any() and not single.all()
    leak = {}
    for mode in ("true", "false"):
        u = probe(ball, mode, pts)
        un = np.abs(np.sum(u * normals, axis=1))
        leak[mode] = un / np.linalg.norm(u - np.sum(u * normals, axis=1)[:, None] * normals, axis=1)
    assert leak["true"][single].max() < 1e-3, leak["true"][single].max()
    for q in (50, 95):
        edge, none = np.percentile(leak["true"], q), np.percentile(leak["false"], q)
        assert edge < 0.125 * none, (q, edge, none)


@pytest.mark.slow
def test_tracers_keep_the_eulerian_mean_velocity_past_a_sphere(sphere):
    """As between the cylinders, tracers in the P1 field collect on the
    sphere's upstream side and stop: by t = 100 most are at rest. With
    wall_upgrade=true few stop and the mean x velocity stays near its Eulerian
    value (u_x's mean does not depend on the cross-stream force)."""
    stuck, ux = {}, {}
    for mode in ("true", "false"):
        t, u, stuck[mode] = lattice_run(sphere, mode, 10, T=100.0)
        ux[mode] = u[:, 0] / u[0, 0]
    late = t > 50
    assert stuck["false"] > 0.5 and ux["false"][late].mean() < 0.5, (stuck, ux["false"][late].mean())
    assert stuck["true"] < 0.1, stuck
    assert ux["true"][late].mean() > 0.85, ux["true"][late].mean()


def cells_against_rk4(case, mode, points, T, dt):
    """tracers from points to T with wall_upgrade=mode at dt under RK4 and
    RK4cells: {scheme: (steps declined, fraction at rest at T, fallbacks a
    particle-step or None)}, and the positions RK4cells stored at its step
    ends that no cell holds. Such a position dumps zero u and p, which no
    position in the fluid does: p is nowhere zero."""
    res, outside = {}, None
    for scheme in ("RK4", "RK4cells"):
        out = []
        every = scheme == "RK4cells"
        args = ["T=%g" % T, "dump_intv=%g" % (dt if every else T), "stat_intv=%g" % T]
        if every:
            args.append("output_all_props=true")
        folder = restart_from(case, mode, points, args, scheme=scheme, out=out, dt=dt)
        u = dump_at(folder, T)["u"]
        [f] = [p for p in folder.glob("tdata_from_t*.dat") if p.name != "tdata_from_t0.000000.dat"]
        res[scheme] = (read_stats(f)["n_declined"][-1], np.mean(np.linalg.norm(u, axis=1) < 1e-3),
                       cells_counts(out[0]).get("fallbacks"))
        if every:
            dumps = all_dumps(folder)
            assert len(dumps) >= T / dt - 2
            outside = sum(int(np.sum((g["p"].ravel() == 0) & (g["u"] == 0).all(axis=1))) for g in dumps.values())
    return res, outside


@pytest.mark.slow
def test_rk4cells_past_the_sphere_declines_as_rk4_does(sphere):
    """RK4cells on the sphere case (a periodic box, the no-slip sphere, stamps
    every 0.25), tracers on the 10^3 lattice to T = 50, with wall_upgrade=true and
    none, dt from 0.05 to 0.5. A step meeting the sphere falls back to RK4's
    pieces and RK4's outside rule, so RK4cells declines no more steps than
    RK4. With wall_upgrade=false tracers come to rest on the sphere under both,
    within the band of its facets and moving into them: a wall never makes
    them cross before they step, so they step held in their cell and the
    fallbacks stay rare. An end below a wall is moved onto it and nudged
    inside, so every position stored at a step's end is in a cell."""
    n, T = 10, 50.0
    x = (np.arange(n) + 0.5) / n
    points = np.stack(np.meshgrid(*[x] * 3, indexing="ij"), axis=-1).reshape(-1, 3)
    points = points[inside(sphere.geometry, points)]
    runs = {(mode, dt): cells_against_rk4(sphere, mode, points, T, dt)
            for mode in ("true", "false") for dt in (DT, 0.2, 0.5)}
    table = "\n".join("%s dt=%g: %s, %d outside" % (k + v) for k, v in runs.items())
    for (mode, dt), (res, outside) in runs.items():
        (d4, rest4, _), (dc, restc, fbc) = res["RK4"], res["RK4cells"]
        assert dc <= d4 and fbc < 1e-3 and outside == 0, table
        if mode == "true":
            assert rest4 < 0.02 and restc < rest4 + 0.01, table
        else:
            assert rest4 > 0.5 and restc > 0.5, table


def test_rk4cells_between_the_cylinders_declines_none(obstacles):
    """The four cylinders in the P1 field (wall_upgrade=false), tracers on a 30^2
    lattice to T = 100 at dt = 0.5: RK4 declines steps of tracers resting on
    the walls when the flow turns them; RK4cells, which keeps every stored
    position in a cell, declines none."""
    n = 30
    x = (np.arange(n) + 0.5) / n
    points = np.stack(np.meshgrid(x, x, indexing="ij"), axis=-1).reshape(-1, 2)
    points = points[inside(obstacles.geometry, points)]
    res, outside = cells_against_rk4(obstacles, "false", points, 100.0, 0.5)
    assert res["RK4cells"][0] == 0 and outside == 0, (res, outside)
