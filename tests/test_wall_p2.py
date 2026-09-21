"""The quadratic velocity XDMFTriangleInterpol and XDMFTetInterpol build near
no-slip walls.

XDMF velocities are P1. In every cell with a vertex on a wall, wall_p2=edge
makes the field P2, each edge midpoint a fixed linear combination of the
velocity at the edge's fluid end, so the field is continuous, and the velocity
normal to a wall facet grows as the square of the distance to it, as the
no-slip asymptotics has it. An edge from a wall vertex above two or more wall
facets at an angle takes their mean normal, so near those facets the normal
velocity keeps a small linear part. wall_p2=none leaves the field P1, whose
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
  sphere along ever different paths.

A velocity is probed by restarting tracers from a checkpoint whose positions
are the probe points: the restarted run dumps the velocity before it steps.
Mean velocities come from the statistics of a run seeded on a lattice, since
random seeding leaves a sampling error as large as the effects at issue.
"""

import os
import shutil
import subprocess

import numpy as np
import pytest

from dumps import dump_at
from paths import app


def need_gmsh():
    """gmsh, or a skip. Its wheel loads shared libraries (libGL, libGLU) a
    container need not have; importing it then raises OSError, which
    importorskip lets through, and a machine without them should skip these
    tests rather than fail them."""
    try:
        import gmsh
        gmsh.initialize()
        gmsh.finalize()
    except (ImportError, OSError) as e:
        pytest.skip("gmsh cannot run here: %s" % e)
    return gmsh


TRACERS = app("tracers")

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(not os.path.exists(TRACERS), reason="tracers is not built"),
]

H = 0.04               # cell size of the gmsh meshes
DT = 0.05

# the four-cylinder case: centres, radius and force direction
OBSTACLES = np.array([(0.583, 0.762), (0.682, 0.319), (0.173, 0.712), (0.312, 0.276)])
OBSTACLE_R = 0.12
OBSTACLE_FORCE = np.array([np.cos(np.radians(20)), np.sin(np.radians(20))])


def periodic_mesh(gmsh, centres, radius):
    """(vertex coordinates, cells) of the box less the cylinders, with the
    nodes on opposite sides of the box paired."""
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        occ = gmsh.model.occ
        box = occ.addRectangle(0, 0, 0, 1, 1)
        disks = [(2, occ.addDisk(x, y, 0, radius, radius)) for x, y in centres]
        occ.cut([(2, box)], disks)
        occ.synchronize()

        def side(x=None, y=None):
            for _, c in gmsh.model.getEntities(1):
                b = gmsh.model.getBoundingBox(1, c)
                if ((x is None or (abs(b[0] - x) < 1e-5 and abs(b[3] - x) < 1e-5)) and
                        (y is None or (abs(b[1] - y) < 1e-5 and abs(b[4] - y) < 1e-5))):
                    return c

        shift_x = [1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
        shift_y = [1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1]
        gmsh.model.mesh.setPeriodic(1, [side(x=1)], [side(x=0)], shift_x)
        gmsh.model.mesh.setPeriodic(1, [side(y=1)], [side(y=0)], shift_y)
        gmsh.option.setNumber("Mesh.MeshSizeMax", H)
        gmsh.option.setNumber("Mesh.MeshSizeMin", H)
        gmsh.model.mesh.generate(2)
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        _, _, nodes = gmsh.model.mesh.getElements(2)
    finally:
        gmsh.finalize()
    index = {t: i for i, t in enumerate(tags)}
    X = xyz.reshape(-1, 3)[:, :2]
    cells = np.array([index[t] for t in nodes[0]]).reshape(-1, 3)
    # drop nodes no triangle uses
    used = np.unique(cells)
    renum = np.full(len(X), -1)
    renum[used] = np.arange(len(used))
    return X[used], renum[cells]


def dolfin_mesh(df, X, cells):
    """A dolfin mesh with these vertices and cells (triangles or tets), in this order."""
    dim = X.shape[1]
    mesh = df.Mesh()
    ed = df.MeshEditor()
    ed.open(mesh, "triangle" if dim == 2 else "tetrahedron", dim, dim)
    ed.init_vertices(len(X))
    ed.init_cells(len(cells))
    for i, x in enumerate(X):
        ed.add_vertex(i, x)
    for i, c in enumerate(cells):
        ed.add_cell(i, c)
    ed.close()
    return mesh


def on_cylinders(centres, radius):
    """A predicate for points on any of the cylinders."""
    def wall(X):
        d = np.min([np.hypot(X[:, 0] - x, X[:, 1] - y) for x, y in centres], axis=0)
        return np.abs(d - radius) < 1e-6
    return wall


def stokes_past_cylinders(df, X, cells, centres, radius, force):
    """The P1 vertex velocity and pressure of periodic Stokes flow past the
    cylinders under a uniform body force, scaled to unit mean velocity along
    the force, with exact zeros on the walls."""
    mesh = dolfin_mesh(df, X, cells)
    wall = on_cylinders(centres, radius)

    class Periodic(df.SubDomain):
        def inside(self, x, on):
            return bool(on and (df.near(x[0], 0) or df.near(x[1], 0))
                        and not (df.near(x[0], 1) or df.near(x[1], 1)))

        def map(self, x, y):
            y[0] = x[0] - 1 if df.near(x[0], 1) else x[0]
            y[1] = x[1] - 1 if df.near(x[1], 1) else x[1]

    class Wall(df.SubDomain):
        def inside(self, x, on):
            return bool(on and min(np.hypot(x[0] - cx, x[1] - cy) for cx, cy in centres) < radius + 1e-3)

    cell = mesh.ufl_cell()
    W = df.FunctionSpace(mesh, df.MixedElement([df.VectorElement("CG", cell, 2),
                                               df.FiniteElement("CG", cell, 1)]),
                         constrained_domain=Periodic())
    u, p = df.TrialFunctions(W)
    v, q = df.TestFunctions(W)
    a = (df.inner(df.grad(u), df.grad(v)) - p * df.div(v) - q * df.div(u)) * df.dx
    L = df.inner(df.Constant(tuple(force)), v) * df.dx
    bcs = [df.DirichletBC(W.sub(0), df.Constant((0, 0)), Wall()),
           df.DirichletBC(W.sub(1), df.Constant(0), "near(x[0], 0) && near(x[1], 0)", "pointwise")]
    w = df.Function(W)
    df.solve(a == L, w, bcs)
    u2, p1 = w.split(deepcopy=True)
    scale = (df.assemble(df.Constant(1) * df.dx(domain=mesh))
             / df.assemble((u2[0] * force[0] + u2[1] * force[1]) * df.dx))

    # the P1 output is written from a constrained space too, so a vertex and
    # its image hold one value, as the case's parameter file claims
    V1 = df.VectorFunctionSpace(mesh, "CG", 1, constrained_domain=Periodic())
    u1 = df.interpolate(u2, V1)
    vals = scale * u1.vector().get_local()
    vals[wall(V1.tabulate_dof_coordinates())] = 0.
    u1.vector().set_local(vals)
    p1 = df.interpolate(p1, df.FunctionSpace(mesh, "CG", 1, constrained_domain=Periodic()))
    p1.vector()[:] *= scale
    return {"u": u1, "p": p1}


def write_case(df, base, fields, periodic, wall, seed_point, stamps=None):
    """Write the fields as an XDMF case under base; returns case(mode), a
    folder holding it read with wall_p2=mode. case.geometry is the mesh's
    (vertex coordinates, cells), case.wall tells which points are on a wall,
    and case.seed_point is a point in the fluid.

    fields maps a name to a function, written unchanged at t = 0 and 1e4, or,
    with stamps, to a callable t -> function written at each of the stamps.
    """
    cwd = os.getcwd()
    os.chdir(base)
    try:
        for name, f in fields.items():
            xf = df.XDMFFile(name + ".xdmf")
            xf.parameters["functions_share_mesh"] = True
            xf.parameters["rewrite_function_mesh"] = False
            for t in (stamps if stamps is not None else (0.0, 1e4)):
                xf.write(f(t) if stamps is not None else f, t)
            xf.close()
    finally:
        os.chdir(cwd)

    def get(mode):
        d = base / mode
        if not d.exists():
            d.mkdir()
            for f in ("u.xdmf", "u.h5", "p.xdmf", "p.h5"):
                shutil.copy(base / f, d / f)
            (d / "dolfin_params.dat").write_text(
                "u=u.xdmf\np=p.xdmf\n" + "".join("periodic_%s=%s\n" % (a, v) for a, v in zip("xyz", periodic))
                + "wall_p2=%s\n" % mode)
        return d

    u = fields["u"](0.0) if stamps is not None else fields["u"]
    mesh = u.function_space().mesh()
    get.geometry = (mesh.coordinates().copy(), mesh.cells().copy())
    get.mode = "xdmftriangle" if mesh.geometry().dim() == 2 else "xdmftet"
    get.wall = wall
    get.seed_point = seed_point
    return get


@pytest.fixture(scope="module")
def cylinder(tmp_path_factory):
    """The single-cylinder case, as write_case returns it."""
    df = pytest.importorskip("dolfin", reason="writing XDMF needs dolfin")
    gmsh = need_gmsh()
    centres, radius = [(0.5, 0.5)], 0.25
    X, cells = periodic_mesh(gmsh, centres, radius)
    fields = stokes_past_cylinders(df, X, cells, centres, radius, np.array([1.0, 0.0]))
    return write_case(df, tmp_path_factory.mktemp("cylinder"), fields,
                      ("true", "true"), on_cylinders(centres, radius), (0.5, 0.1))


@pytest.fixture(scope="module")
def obstacles(tmp_path_factory):
    """The four-cylinder case, as write_case returns it."""
    df = pytest.importorskip("dolfin", reason="writing XDMF needs dolfin")
    gmsh = need_gmsh()
    X, cells = periodic_mesh(gmsh, OBSTACLES, OBSTACLE_R)
    fields = stokes_past_cylinders(df, X, cells, OBSTACLES, OBSTACLE_R, OBSTACLE_FORCE)
    return write_case(df, tmp_path_factory.mktemp("obstacles"), fields,
                      ("true", "true"), on_cylinders(OBSTACLES, OBSTACLE_R), (0.5, 0.1))


def kink_wall(x):
    """Height of the kinked channel's bottom wall."""
    return 0.1 * (0.5 - np.abs(x - 0.5))


def on_kinked_walls(X):
    """Which points lie on the kinked channel's walls."""
    return (np.abs(X[:, 1] - kink_wall(X[:, 0])) < 1e-12) | (np.abs(X[:, 1] - 1) < 1e-12)


def kinked_channel_mesh(nx=16, ny=8):
    """(vertex coordinates, cells) of the kinked channel: columns of nodes
    between the bottom wall and y = 1, squares split along the diagonal that
    rises towards the kink, so both cells on the kink's facets have the node
    above the kink as their third vertex."""
    X = []
    for i in range(nx + 1):
        x = i / nx
        yb = kink_wall(x)
        X += [(x, yb + (1 - yb) * j / ny) for j in range(ny + 1)]
    node = lambda i, j: i * (ny + 1) + j
    cells = []
    for i in range(nx):
        for j in range(ny):
            a, b, c, d = node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1)
            if i < nx // 2:
                cells += [(a, b, c), (a, c, d)]
            else:
                cells += [(a, b, d), (b, c, d)]
    return np.array(X), np.array(cells)


@pytest.fixture(scope="module")
def kinked(tmp_path_factory):
    """The kinked channel case, as write_case returns it. The wall velocities
    are 1e-17 rather than zero, as round-off in solver output leaves them."""
    df = pytest.importorskip("dolfin", reason="writing XDMF needs dolfin")
    X, cells = kinked_channel_mesh()
    mesh = dolfin_mesh(df, X, cells)
    V = df.VectorFunctionSpace(mesh, "CG", 1)
    x, y = X[:, 0], X[:, 1]
    g = (y - kink_wall(x)) * (1 - y)
    U = np.c_[4 * g * (1 + 0.5 * np.sin(2 * np.pi * x)), 1.2 * g * np.cos(2 * np.pi * x)]
    U[on_kinked_walls(X)] = 1e-17
    u = df.Function(V)
    vals = np.zeros(V.dim())
    vals[df.vertex_to_dof_map(V)] = U.ravel()
    u.vector().set_local(vals)
    p = df.interpolate(df.Constant(0.), df.FunctionSpace(mesh, "CG", 1))
    return write_case(df, tmp_path_factory.mktemp("kinked"), {"u": u, "p": p},
                      ("true", "false"), on_kinked_walls, (0.5, 0.5))


SPHERE_C = np.array([0.5, 0.5, 0.5])
SPHERE_R = 0.25
SPHERE_H = 0.1


def crossflow(t):
    """The sphere case's force at time t."""
    return np.array([1.0, 0.5 * np.sin(2 * np.pi * t / 7), 0.5 * np.sin(2 * np.pi * t / 11 + 1)])


def periodic_sphere_mesh(gmsh):
    """(vertex coordinates, cells) of the unit cube less the sphere, with the
    nodes on opposite faces paired."""
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        occ = gmsh.model.occ
        box = occ.addBox(0, 0, 0, 1, 1, 1)
        ball = occ.addSphere(*SPHERE_C, SPHERE_R)
        occ.cut([(3, box)], [(3, ball)])
        occ.synchronize()

        def face(axis, value):
            for _, f in gmsh.model.getEntities(2):
                b = gmsh.model.getBoundingBox(2, f)
                if abs(b[axis] - value) < 1e-5 and abs(b[axis + 3] - value) < 1e-5:
                    return f

        for axis in range(3):
            shift = np.eye(4)
            shift[axis, 3] = 1
            gmsh.model.mesh.setPeriodic(2, [face(axis, 1)], [face(axis, 0)], shift.ravel().tolist())
        gmsh.option.setNumber("Mesh.MeshSizeMax", SPHERE_H)
        gmsh.option.setNumber("Mesh.MeshSizeMin", SPHERE_H)
        gmsh.model.mesh.generate(3)
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        _, _, nodes = gmsh.model.mesh.getElements(3)
    finally:
        gmsh.finalize()
    index = {t: i for i, t in enumerate(tags)}
    X = xyz.reshape(-1, 3)
    cells = np.array([index[t] for t in nodes[0]]).reshape(-1, 4)
    used = np.unique(cells)
    renum = np.full(len(X), -1)
    renum[used] = np.arange(len(used))
    return X[used], renum[cells]


def on_sphere(X):
    """Which points lie on the sphere."""
    return np.abs(np.linalg.norm(X - SPHERE_C, axis=1) - SPHERE_R) < 1e-6


@pytest.fixture(scope="module")
def sphere(tmp_path_factory):
    """The sphere case, as write_case returns it: the Stokes solutions for a
    unit force along each axis, combined at each stamp (every 0.25 up to
    t = 110) with the weights of crossflow(t), and scaled to unit mean u_x
    under a unit force along x."""
    df = pytest.importorskip("dolfin", reason="writing XDMF needs dolfin")
    gmsh = need_gmsh()
    X, cells = periodic_sphere_mesh(gmsh)
    mesh = dolfin_mesh(df, X, cells)

    class Periodic(df.SubDomain):
        def inside(self, x, on):
            return bool(on and any(df.near(x[k], 0) for k in range(3))
                        and not any(df.near(x[k], 1) for k in range(3)))

        def map(self, x, y):
            for k in range(3):
                y[k] = x[k] - 1 if df.near(x[k], 1) else x[k]

    class Wall(df.SubDomain):
        def inside(self, x, on):
            return bool(on and np.linalg.norm(np.array(x) - SPHERE_C) < SPHERE_R + 1e-3)

    cell = mesh.ufl_cell()
    W = df.FunctionSpace(mesh, df.MixedElement([df.VectorElement("CG", cell, 2),
                                               df.FiniteElement("CG", cell, 1)]),
                         constrained_domain=Periodic())
    u, p = df.TrialFunctions(W)
    v, q = df.TestFunctions(W)
    a = (df.inner(df.grad(u), df.grad(v)) - p * df.div(v) - q * df.div(u)) * df.dx
    bcs = [df.DirichletBC(W.sub(0), df.Constant((0, 0, 0)), Wall()),
           df.DirichletBC(W.sub(1), df.Constant(0), "near(x[0], 0) && near(x[1], 0) && near(x[2], 0)",
                          "pointwise")]
    # the P1 output is written from a constrained space too, so a vertex and
    # its image hold one value, as the case's parameter file claims
    V1 = df.VectorFunctionSpace(mesh, "CG", 1, constrained_domain=Periodic())
    Q1 = df.FunctionSpace(mesh, "CG", 1, constrained_domain=Periodic())
    wall = on_sphere(V1.tabulate_dof_coordinates())
    solver = None
    U, P = [], []
    for k in range(3):
        L = df.inner(df.Constant(tuple(np.eye(3)[k])), v) * df.dx
        A, b = df.assemble_system(a, L, bcs)
        solver = solver or df.LUSolver(A, "mumps")
        w = df.Function(W)
        solver.solve(w.vector(), b)
        u2, p2 = w.split(deepcopy=True)
        vals = df.interpolate(u2, V1).vector().get_local()
        vals[wall] = 0.
        U.append(vals)
        P.append(df.interpolate(p2, Q1).vector().get_local())
        if k == 0:
            vol = df.assemble(df.Constant(1) * df.dx(domain=mesh))
            scale = vol / df.assemble(u2[0] * df.dx)
    uf, pf = df.Function(V1), df.Function(Q1)

    def velocity(t):
        uf.vector().set_local(scale * crossflow(t) @ np.array(U))
        return uf

    def pressure(t):
        pf.vector().set_local(scale * crossflow(t) @ np.array(P))
        return pf

    return write_case(df, tmp_path_factory.mktemp("sphere"), {"u": velocity, "p": pressure},
                      ("true", "true", "true"), on_sphere, (0.5, 0.1, 0.5),
                      stamps=np.arange(0, 110.01, 0.25))


@pytest.fixture(params=["cylinder", "kinked"])
def case(request):
    """Each case in turn."""
    return request.getfixturevalue(request.param)


def run(case, args):
    """Run tracers on the case's parameter file with args; assert it succeeded."""
    r = subprocess.run([TRACERS, str(case / "dolfin_params.dat")] + args,
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr


def restart_from(case, mode, points, args):
    """Run tracers on case(mode) from a checkpoint holding `points`; returns
    the run folder.

    A short run with positions drawn around case.seed_point (from a fixed
    seed) writes the checkpoint, at t = 2 DT and step 2, and its positions
    are overwritten; the restart continues from step 2.
    """
    folder = case(mode)
    if (folder / "Tracers").exists():
        shutil.rmtree(folder / "Tracers")
    dim = points.shape[1]
    seed = list(case.seed_point) + [0.0] * (3 - dim)
    base = ["mode=" + case.mode, "init_mode=points_" + "xyz"[:dim], "int_order=1", "Dm=0",
            "dt=%g" % DT, "Nrw=%d" % len(points), "Nrw_max=%d" % len(points),
            "x0=%g" % seed[0], "y0=%g" % seed[1], "z0=%g" % seed[2],
            "random=false", "seed=1", "scheme=RK4"]
    run(folder, base + ["T=%g" % DT, "dump_intv=1000", "stat_intv=1000",
                        "checkpoint_intv=%g" % DT])
    [pos] = list(folder.rglob("Checkpoints/positions.pos"))
    np.savetxt(pos, np.c_[points, np.zeros((len(points), 3 - dim))], fmt="%.17g")
    run(folder, base + args + ["checkpoint_intv=1e9", "restart_folder=" + str(pos.parent.parent)])
    return pos.parent.parent


def probe(case, mode, points):
    """The velocity at each point, in the order given, with wall_p2=mode."""
    folder = restart_from(case, mode, points,
                          ["T=%g" % (2 * DT), "dump_intv=%g" % (2 * DT), "stat_intv=1000"])
    ids = np.loadtxt(folder / "Checkpoints" / "id.list", dtype=int)
    g = dump_at(folder, 2 * DT)
    # dump_at orders by id; the probe points were written in id.list order
    dim = points.shape[1]
    assert np.allclose(g["points"][ids][:, :dim], points, rtol=0, atol=1e-12)
    return g["u"][ids][:, :dim]


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


def inside(geometry, points):
    """Which points lie in a mesh cell."""
    X, cells = geometry
    found = np.zeros(len(points), dtype=bool)
    for c in cells:
        T = (X[c[1:]] - X[c[0]]).T
        lam = np.linalg.solve(T, (points - X[c[0]]).T)
        found |= (lam >= 0).all(axis=0) & (lam.sum(axis=0) <= 1)
    return found


def test_the_field_is_continuous_across_edges_near_the_wall(case):
    """Across every edge of a cell with a wall vertex, points a hair apart on
    either side see the same velocity with wall_p2=edge. The quadratic field
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
    u = probe(case, "edge", pts)
    assert np.abs(u[0::2] - u[1::2]).max() < 1e-6
    assert np.abs(u - probe(case, "none", pts)).max() > 1e-3


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
    delta = 1e-5, 1e-4 and 1e-3 with wall_p2=edge: the velocity normal to the
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
    for mode in ("edge", "none"):
        u = probe(case, mode, pts)
        un = np.sum(u * normals, axis=1)
        ratio = (un / dist**2).reshape(-1, len(deltas))
        spread = np.abs(ratio - ratio[:, -1:]).max(axis=1)
        if mode == "edge":
            exact = ~shared
            assert spread[exact].max() <= 1e-3 * np.abs(ratio[exact]).max(), (mode, spread[exact].max())
        else:
            assert spread.max() > 10 * np.abs(ratio[:, -1]).max(), (mode, spread.max())
        ut = np.linalg.norm(u - un[:, None] * normals, axis=1)
        leak[mode] = (np.abs(un) / ut)[0::len(deltas)]
    if shared.any():
        assert leak["edge"][shared].max() < 0.1, leak["edge"][shared].max()
        assert (leak["edge"][shared] < 0.5 * leak["none"][shared]).all(), (leak["edge"][shared], leak["none"][shared])


def lattice_run(case, mode, n, T):
    """Tracers on a cell-centred lattice with n points a side, those in the
    fluid, run to T with wall_p2=mode: (times, mean velocity over time, the
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


def test_tracers_keep_the_eulerian_mean_velocity_between_obstacles(obstacles):
    """Tracers seeded uniformly in an incompressible flow with impermeable
    walls stay uniform, so their mean velocity stays at the Eulerian mean.
    In the P1 field each pass of an obstacle moves a tracer a little towards
    its upstream face, where it stops: by t = 200 most tracers are at rest
    and the mean has collapsed. With wall_p2=edge the mean stays within a few
    percent (its sampling noise is about 2% with 2000 tracers) and next to
    nothing stops."""
    ratio, stuck = {}, {}
    for mode in ("edge", "none"):
        _, u, stuck[mode] = lattice_run(obstacles, mode, 50, T=200.0)
        along = u[:, :2] @ OBSTACLE_FORCE
        ratio[mode] = along / along[0]
    assert stuck["none"] > 0.5 and ratio["none"][-1] < 0.3, (stuck, ratio["none"][-1])
    assert stuck["edge"] < 0.01, stuck
    late = ratio["edge"][len(ratio["edge"]) // 2:]
    assert late.mean() > 0.9 and late[-1] > 0.9, (late.mean(), late[-1])


def tet_faces(case):
    """{face (sorted vertex indices): the vertices opposite it} over all tets."""
    X, cells = case.geometry
    faces = {}
    for c in cells:
        for k in range(4):
            faces.setdefault(tuple(sorted(np.delete(c, k))), []).append(c[k])
    return faces


def test_the_tet_field_is_continuous_across_faces_near_the_wall(sphere):
    """Across every face of a tet with a wall vertex, points a hair apart on
    either side see the same velocity with wall_p2=edge, and the quadratic
    field differs from the P1 one there."""
    X, _ = sphere.geometry
    wall = sphere.wall(X)
    pts = []
    for f, opposite in tet_faces(sphere).items():
        if len(opposite) != 2 or not wall[list(f)].any() or wall[list(f)].all():
            continue
        a, b, c = X[list(f)]
        n = np.cross(b - a, c - a)
        n /= np.linalg.norm(n)
        for w in ((1 / 3, 1 / 3, 1 / 3), (0.5, 0.25, 0.25), (0.25, 0.5, 0.25)):
            x = w[0] * a + w[1] * b + w[2] * c
            pts += [x + 1e-9 * n, x - 1e-9 * n]
    pts = np.array(pts)
    u = probe(sphere, "edge", pts)
    assert np.abs(u[0::2] - u[1::2]).max() < 1e-6
    assert np.abs(u - probe(sphere, "none", pts)).max() > 1e-3


def test_tet_facet_normal_velocity_is_nearly_quadratic(sphere):
    """Just above the sphere's facets (delta = 1e-5) the velocity through a
    facet is a small fraction of the velocity along it with wall_p2=edge:
    below 1e-3, of order delta as the quadratic law has it, where the facet's
    cell has no side edge above another facet, and elsewhere, in median and at
    the 95th percentile, a quarter or less of what it is in the P1 field
    (about a seventh on this mesh). Tets with all four vertices on the sphere
    carry no flow and are left out."""
    X, _ = sphere.geometry
    wall = sphere.wall(X)
    facets = [(f, o[0]) for f, o in tet_faces(sphere).items()
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
    for mode in ("edge", "none"):
        u = probe(sphere, mode, pts)
        un = np.abs(np.sum(u * normals, axis=1))
        leak[mode] = un / np.linalg.norm(u - np.sum(u * normals, axis=1)[:, None] * normals, axis=1)
    assert leak["edge"][single].max() < 1e-3, leak["edge"][single].max()
    for q in (50, 95):
        edge, none = np.percentile(leak["edge"], q), np.percentile(leak["none"], q)
        assert edge < 0.25 * none, (q, edge, none)


def test_tracers_keep_the_eulerian_mean_velocity_past_a_sphere(sphere):
    """As between the cylinders, tracers in the P1 field collect on the
    sphere's upstream side and stop: by t = 100 most are at rest. With
    wall_p2=edge few stop and the mean x velocity stays near its Eulerian
    value (u_x's mean does not depend on the cross-stream force)."""
    stuck, ux = {}, {}
    for mode in ("edge", "none"):
        t, u, stuck[mode] = lattice_run(sphere, mode, 10, T=100.0)
        ux[mode] = u[:, 0] / u[0, 0]
    late = t > 50
    assert stuck["none"] > 0.5 and ux["none"][late].mean() < 0.5, (stuck, ux["none"][late].mean())
    assert stuck["edge"] < 0.1, stuck
    assert ux["edge"][late].mean() > 0.85, ux["edge"][late].mean()
