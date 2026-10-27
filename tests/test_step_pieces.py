"""A step that leaves the fluid is taken again in pieces, and what is still
declined is counted and reported.

The fluid is a regular 12-gon of circumradius 1 (2D, mode=triangle) or the
prism over it of height 1 (3D, mode=tet), meshed as a fan about the axis, so
every wall edge is convex as seen from the fluid. The fields are linear, which
P1 holds exactly: the rigid rotation u = (-y, x[, 0]) and the outflow
u = (x, y[, 0]).

Rotating tracers start on a circle inside the 12-gon's inner circle. A step
of dt = 0.5 turns them by about 29 degrees. From radius 0.95 an RK4 step lands
inside, but its stage x + k1 dt/2 lies a factor 1.03 further out, beyond the
plane of a wall facet for every tracer whose stage falls within some 9 degrees
of a facet's middle. An explicit Euler step ends a factor 1.12 further out, so
from radius 0.92 beyond every wall. Each such step must be taken in pieces:
halves, a failed piece halved again, down to dt/64. The positions after one
step are recomputed here with the same rule from the exact field, so they are
pinned to round-off; and every tracer must move, end inside, and keep its line
element of unit length or its frame finite.

A march, tracervectors_spatial, steps by path length ds in the rotation
frozen at t, at both orders, from where a step ends beyond every wall. Each is
taken in halves of ds, every piece's time from the speed at its start, and the
ends and integration times are recomputed here the same way.

Outflowing tracers leave the fluid within the step, so a piece of dt/64 fails
too: their steps are declined and go to outside=, and the run reports the total
declined and what outside= did with them, at the end and, verbose, on every
progress line.
"""

import math
import os

import numpy as np
import pytest

from dumps import all_dumps, read_stats
from paths import app
from pieces import in_pieces
from runs import checkpoint_folder, put_points, run_app
from wall_cases import dolfin_mesh

APPS = ["tracers", "tracervectors", "tracertensors"]
N_SIDES = 12
RHO = math.cos(math.pi / N_SIDES)    # the 12-gon's inner radius
R0 = 0.95                            # the tracers' start radius, RK4
R0_EULER = 0.92                      # explicit Euler
DT = 0.5
N_TRACERS = 96

pytestmark = pytest.mark.skipif(not all(os.path.exists(app(a)) for a in APPS),
                                reason="the tracer apps are not built")


def facet_normals():
    """Outward unit normals of the 12-gon's sides, whose vertices lie at multiples of 30 degrees."""
    a = (np.arange(N_SIDES) + 0.5) * 2 * np.pi / N_SIDES
    return np.c_[np.cos(a), np.sin(a)]


NORMALS = facet_normals()


def inside(x):
    """x in the fluid: inside every side's half-plane, and 0 <= z <= 1 (z = 0 in 2D)."""
    x = np.asarray(x, float)
    return np.all(NORMALS @ x[:2] <= RHO * (1 + 1e-14)) and 0 <= x[2] <= 1


def rotation(x):
    u = np.zeros_like(x)
    u[0], u[1] = -x[1], x[0]
    return u


def outflow(x):
    u = np.zeros_like(x)
    u[0], u[1] = x[0], x[1]
    return u


def fan_mesh(dim):
    """(vertices, cells): the 12-gon as triangles about its centre, or the prism
    over it as three tets a wedge, split so neighbouring wedges share their
    radial faces' diagonals (centre to the upper rim vertex)."""
    a = np.arange(N_SIDES) * 2 * np.pi / N_SIDES
    rim = np.c_[np.cos(a), np.sin(a)]
    nxt = lambda i: (i + 1) % N_SIDES
    if dim == 2:
        X = np.vstack([[0., 0.], rim])
        cells = [(0, 1 + i, 1 + nxt(i)) for i in range(N_SIDES)]
        return X, np.array(cells)
    # 0, 1: the axis at z = 0, 1; 2 + i, 2 + N + i: rim vertex i at z = 0, 1
    X = np.vstack([[0., 0., 0.], [0., 0., 1.], np.c_[rim, np.zeros(N_SIDES)], np.c_[rim, np.ones(N_SIDES)]])
    lo = lambda i: 2 + i
    hi = lambda i: 2 + N_SIDES + i
    cells = []
    for i in range(N_SIDES):
        j = nxt(i)
        cells += [(0, lo(i), lo(j), hi(j)), (0, lo(i), hi(j), hi(i)), (0, hi(i), hi(j), 1)]
    return X, np.array(cells)


def write_fan_case(folder, dim, field):
    """A dolfin HDF5 case of the fan mesh holding the P1 field u = field(x), at
    two stamps: its parameter file."""
    df = pytest.importorskip("dolfin", reason="writing the mesh needs dolfin")
    folder.mkdir(parents=True)
    mesh = dolfin_mesh(df, *fan_mesh(dim))
    V = df.VectorFunctionSpace(mesh, "CG", 1)
    u = df.Function(V)
    vals = np.zeros(V.dim())
    vals[df.vertex_to_dof_map(V)] = np.array([field(x) for x in mesh.coordinates()]).ravel()
    u.vector().set_local(vals)
    p = df.interpolate(df.Constant(0.), df.FunctionSpace(mesh, "CG", 1))
    with df.HDF5File(mesh.mpi_comm(), str(folder / "mesh.h5"), "w") as h:
        h.write(mesh, "mesh")
    with df.HDF5File(mesh.mpi_comm(), str(folder / "up.h5"), "w") as h:
        h.write(u, "u")
        h.write(p, "p")
    (folder / "timestamps.dat").write_text("0 up.h5\n1e4 up.h5\n")
    (folder / "dolfin_params.dat").write_text(
        "velocity_space=P1\npressure_space=P1\ntimestamps=timestamps.dat\nmesh=mesh.h5\n"
        "periodic_x=false\nperiodic_y=false\nperiodic_z=false\nrho=1.0\n")
    return folder / "dolfin_params.dat"


def ring(dim, r):
    """N_TRACERS points on the circle of radius r, at z = 0.5 in 3D."""
    a = (np.arange(N_TRACERS) + 0.25) * 2 * np.pi / N_TRACERS
    return np.c_[r * np.cos(a), r * np.sin(a), np.full(N_TRACERS, 0.5 if dim == 3 else 0.)]


def rk4_once(x, h):
    """One RK4 step of the rotation: (end, every stage and the end inside)."""
    k1 = rotation(x)
    k2 = rotation(x + k1 * h / 2)
    k3 = rotation(x + k2 * h / 2)
    k4 = rotation(x + k3 * h)
    end = x + (k1 + 2 * k2 + 2 * k3 + k4) * h / 6
    return end, all(inside(p) for p in (x, x + k1 * h / 2, x + k2 * h / 2, x + k3 * h, end))


def euler_once(x, h):
    """One explicit Euler step of the rotation: (end, start and end inside)."""
    end = x + rotation(x) * h
    return end, inside(x) and inside(end)


def run_ring(tmp_path, name, dim, field, extra, r0=R0, march=None):
    """Run app name from the ring on the fan mesh holding field: a run places
    and checkpoints particles, the ring replaces their positions, and the
    resumed run dumps at its start and after one step, of DT or, for a march
    app, of path length march. Returns (the resumed process, {0 and the step:
    {dataset: array in id order}}, its statistics)."""
    params = write_fan_case(tmp_path / "case", dim, field)
    case = params.parent
    step, end = ("dxn", "Ln") if march else ("dt", "T")
    h = march or DT
    common = ["mode=" + ("triangle" if dim == 2 else "tet"), "Nrw=%d" % N_TRACERS,
              "Nrw_max=%d" % N_TRACERS, "Dm=0", "int_order=1", "%s=%g" % (step, h), "random=false", "seed=1"]
    run_app(app(name), params, common,
            ["init_mode=points_xy" + ("z" if dim == 3 else ""), "x0=0", "y0=0", "z0=0.5",
             "%s=0" % end, "dump_intv=1e9", "stat_intv=1e9", "checkpoint_intv=1e9", "outside=ignore"])
    put_points(case, ring(dim, r0))
    r = run_app(app(name), params, common,
                ["init_mode=points_xy", "%s=%g" % (end, 2 * h), "dump_intv=%g" % h, "stat_intv=%g" % h,
                 "checkpoint_intv=1e9", "restart_folder=" + str(checkpoint_folder(case))], extra)
    out = {}
    for t, g in all_dumps(case, raw=True).items():
        if t < h / 2:
            continue                                  # the placing run's
        ids = g.pop("id")[:, 0]
        order = np.argsort(ids, kind="stable")
        out[round(t - h, 12)] = {k: a[order] for k, a in g.items()}
    stats = [f for f in case.rglob("tdata_from_t*.dat") if not f.name.startswith("tdata_from_t0.000000")]
    assert len(stats) == 1, stats
    return r, out, read_stats(stats[0])


def check_element(name, g):
    """The carried element after the step: a unit line element, or a finite frame."""
    if name == "tracervectors":
        assert np.allclose(np.linalg.norm(g["n"], axis=1), 1, rtol=0, atol=1e-12)
    if name == "tracertensors":
        for k, a in g.items():
            assert np.all(np.isfinite(a)), k


# Per scheme: its arguments, the ring's radius, one step of it here, and
# whether enough tracers needed pieces for the ring to test the rule
IN_PIECES = {
    "RK4": (["scheme=RK4", "outside=ignore"], R0, rk4_once,
            lambda pieces: sum(m > 1 for m in pieces) > N_TRACERS // 4),
    "explicit": (["scheme=explicit", "int_order=1", "outside=ignore"], R0_EULER, euler_once,
                 lambda pieces: min(pieces) > 1),
}


@pytest.mark.parametrize("scheme", list(IN_PIECES))
@pytest.mark.parametrize("dim", [2, 3], ids=["2d", "3d"])
@pytest.mark.parametrize("name", APPS)
def test_a_step_leaving_the_fluid_is_taken_in_pieces(tmp_path, name, dim, scheme):
    """One step of the rotation that leaves the fluid -- RK4 from radius 0.95,
    a stage crossing a convex wall edge, or explicit Euler at Dm = 0 from
    radius 0.92, ending a factor 1.12 further out, beyond every wall -- is
    taken again in halves, a failed half halved again. Every tracer moves to
    where the rule of halves puts it, inside the fluid, and enough of them
    needed pieces (under Euler every one) that the ring tests the rule."""
    extra, r0, once, enough = IN_PIECES[scheme]
    r, D, stats = run_ring(tmp_path, name, dim, rotation, extra, r0)
    p0, p1 = D[0.0]["points"], D[DT]["points"]
    expected, pieces = [], []
    for x in p0:
        end, m = in_pieces(once, x, DT)
        assert end is not None, x
        expected.append(end)
        pieces.append(m)
    assert enough(pieces), np.bincount(pieces)
    assert np.allclose(p1, expected, rtol=0, atol=1e-12)
    assert all(inside(x) for x in p1)
    assert np.all(np.linalg.norm(p1 - p0, axis=1) > 0.4)
    check_element(name, D[DT])
    assert list(stats["n_declined"]) == [0, 0]


def spatial_once(int_order):
    """One march step of the rotation on (x, the integration time) over path
    length ds: h = ds/|u|, dx = u h, + J u h^2/2 at second order (a = 0, the
    field is steady). Returns (the end and the time plus h, the end inside)."""
    def once(xt, ds):
        x = xt[:3]
        u = rotation(x)
        h = ds / np.linalg.norm(u)
        dx = u * h
        if int_order >= 2:
            dx = dx + 0.5 * rotation(u) * h * h      # J u: the rotation of u
        return np.r_[x + dx, xt[3] + h], inside(x + dx)
    return once


# Per order: a march step's path length from radius 0.92, whose end lies a
# factor 1.12 further out, beyond every wall: h = 0.5 at first order, h = 1 at
# second, where dx = u h - x h^2/2
MARCH_DS = {1: 0.46, 2: 0.92}


@pytest.mark.skipif(not os.path.exists(app("tracervectors_spatial")),
                    reason="tracervectors_spatial is not built")
@pytest.mark.parametrize("int_order", [1, 2])
@pytest.mark.parametrize("dim", [2, 3], ids=["2d", "3d"])
def test_a_march_step_leaving_the_fluid_is_taken_in_pieces(tmp_path, dim, int_order):
    """A march step of the rotation that ends beyond every wall is taken
    again in halves of its path length, the fields frozen at t,
    each piece's time h = ds_piece/|u| from the speed at its start. Every
    tracer moves to where that rule puts it, inside the fluid, its integration
    time advanced by the pieces' h, its line element of unit length, and no
    step is declined; else a march near a convex wall would stall or lose
    tracers a smaller step keeps."""
    ds = MARCH_DS[int_order]
    r, D, stats = run_ring(tmp_path, "tracervectors_spatial", dim, rotation,
                           ["int_order=%d" % int_order, "outside=ignore"], R0_EULER, march=ds)
    a, b = D[0.0], D[ds]
    expected, pieces = [], []
    for x, t in zip(a["points"], a["tau"][:, 0]):
        end, m = in_pieces(spatial_once(int_order), np.r_[x, t], ds)
        assert end is not None, x
        expected.append(end)
        pieces.append(m)
    expected = np.array(expected)
    assert min(pieces) > 1, np.bincount(pieces)
    assert np.allclose(b["points"], expected[:, :3], rtol=0, atol=1e-12)
    assert np.allclose(b["tau"][:, 0], expected[:, 3], rtol=0, atol=1e-12)
    assert all(inside(x) for x in b["points"])
    check_element("tracervectors", b)   # the spatial march carries tracervectors' element
    assert list(stats["n_declined"]) == [0, 0]


@pytest.mark.parametrize("scheme", ["RK4", "RK4cells", "explicit"])
@pytest.mark.parametrize("dim", [2, 3], ids=["2d", "3d"])
def test_a_step_that_leaves_in_every_piece_is_declined_and_reported(tmp_path, dim, scheme):
    """In the outflow every tracer on the ring leaves the fluid within the
    step, so a piece of dt/64 fails too: each step is declined, outside=mark
    marks it (c = 2) where it started, the statistics count it, and the run
    reports the totals at the end and, verbose, on each progress line."""
    r, D, stats = run_ring(tmp_path, "tracers", dim, outflow,
                           ["scheme=" + scheme, "int_order=1", "outside=mark", "verbose=true"])
    p0, p1 = D[0.0]["points"], D[DT]["points"]
    assert np.array_equal(p0, p1)
    assert np.all(D[DT]["c"] == 2)
    assert list(stats["n_declined"]) == [0, N_TRACERS]
    # the run steps once more after its last dump
    assert "Time = %g, declined 0\n" % DT in r.stdout, r.stdout[-2000:]
    assert "Time = %g, declined %d\n" % (2 * DT, N_TRACERS) in r.stdout, r.stdout[-2000:]
    assert "Steps declined: %d (outside=mark: marked)" % (2 * N_TRACERS) in r.stdout, r.stdout[-2000:]
    assert "Warning" not in r.stderr
