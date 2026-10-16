"""The near-wall cases the wall tests share: meshes, the fields on them,
written as XDMF cases read with and without wall_upgrade, and the runs that
probe them. The cases are described in test_wall_edge_rule's docstring."""

import os
import shutil
import subprocess

import numpy as np
import pytest

from dumps import dump_at
from paths import app
from runs import checkpoint_folder, put_points, read_checkpoint


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
    folder holding it read with wall_upgrade=mode. case.geometry is the mesh's
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
                + "wall_upgrade=%s\n" % mode)
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
def sphere_mesh():
    """periodic_sphere_mesh, built once."""
    return periodic_sphere_mesh(need_gmsh())


@pytest.fixture(scope="module")
def ball(tmp_path_factory, sphere_mesh):
    """A synthetic field on the sphere case's mesh, as write_case returns it:
    the distance to the sphere times a smooth periodic field, so at rest on the
    sphere up to round-off, and periodic since the sphere is centred in the
    box."""
    df = pytest.importorskip("dolfin", reason="writing XDMF needs dolfin")
    X, cells = sphere_mesh
    mesh = dolfin_mesh(df, X, cells)
    V = df.VectorFunctionSpace(mesh, "CG", 1)
    x, y, z = X.T
    d = np.linalg.norm(X - SPHERE_C, axis=1) - SPHERE_R
    U = d[:, None] * np.c_[1 + 0.5 * np.sin(2 * np.pi * y), 0.6 * np.cos(2 * np.pi * z),
                           0.4 * np.sin(2 * np.pi * x)]
    u = df.Function(V)
    vals = np.zeros(V.dim())
    vals[df.vertex_to_dof_map(V)] = U.ravel()
    u.vector().set_local(vals)
    p = df.interpolate(df.Constant(0.), df.FunctionSpace(mesh, "CG", 1))
    return write_case(df, tmp_path_factory.mktemp("ball"), {"u": u, "p": p},
                      ("true", "true", "true"), on_sphere, (0.5, 0.1, 0.5))


@pytest.fixture(scope="module")
def sphere(tmp_path_factory, sphere_mesh):
    """The sphere case, as write_case returns it: the Stokes solutions for a
    unit force along each axis, combined at each stamp (every 0.25 up to
    t = 110) with the weights of crossflow(t), and scaled to unit mean u_x
    under a unit force along x."""
    df = pytest.importorskip("dolfin", reason="writing XDMF needs dolfin")
    from petsc4py import PETSc
    X, cells = sphere_mesh
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
    bc = df.DirichletBC(W.sub(0), df.Constant((0, 0, 0)), Wall())
    # MINRES preconditioned by GAMG on the velocity Laplacian and Jacobi on the
    # pressure mass; the pressure's constant is the null space
    L = df.inner(df.Constant((1., 0., 0.)), v) * df.dx
    A, _ = df.assemble_system(a, L, bc)
    B, _ = df.assemble_system((df.inner(df.grad(u), df.grad(v)) + p * q) * df.dx, L, bc)
    ones = df.Function(W).vector()
    W.sub(1).dofmap().set(ones, 1.0)
    df.as_backend_type(A).set_nullspace(df.VectorSpaceBasis([ones / ones.norm("l2")]))
    ksp = PETSc.KSP().create()
    ksp.setOptionsPrefix("sphere_")
    ksp.setOperators(df.as_backend_type(A).mat(), df.as_backend_type(B).mat())
    ksp.setType("minres")
    ksp.setTolerances(rtol=1e-8, max_it=1000)
    pc = ksp.getPC()
    pc.setType("fieldsplit")
    pc.setFieldSplitIS(*[(f, PETSc.IS().createGeneral(W.sub(i).dofmap().dofs())) for i, f in enumerate("up")])
    opts = PETSc.Options("sphere_")
    opts["pc_fieldsplit_type"] = "additive"
    opts["fieldsplit_u_ksp_type"] = opts["fieldsplit_p_ksp_type"] = "preonly"
    opts["fieldsplit_u_pc_type"], opts["fieldsplit_p_pc_type"] = "gamg", "jacobi"
    ksp.setFromOptions()
    # the P1 output is written from a constrained space too, so a vertex and
    # its image hold one value, as the case's parameter file claims
    V1 = df.VectorFunctionSpace(mesh, "CG", 1, constrained_domain=Periodic())
    Q1 = df.FunctionSpace(mesh, "CG", 1, constrained_domain=Periodic())
    wall = on_sphere(V1.tabulate_dof_coordinates())
    U, P = [], []
    for k in range(3):
        b = df.assemble(df.inner(df.Constant(tuple(np.eye(3)[k])), v) * df.dx)
        bc.apply(b)
        w = df.Function(W)
        ksp.solve(df.as_backend_type(b).vec(), df.as_backend_type(w.vector()).vec())
        assert ksp.getConvergedReason() > 0
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


def run(case, args):
    """Run tracers on the case's parameter file with args; assert it succeeded; its stdout."""
    r = subprocess.run([TRACERS, str(case / "dolfin_params.dat")] + args,
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr
    return r.stdout


def restart_from(case, mode, points, args, scheme="RK4", out=None, dt=DT):
    """Run tracers on case(mode) from a checkpoint holding `points`, stepped
    by scheme at dt; returns the run folder, and appends the run's stdout to out.

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
            "dt=%g" % dt, "Nrw=%d" % len(points), "Nrw_max=%d" % len(points),
            "x0=%g" % seed[0], "y0=%g" % seed[1], "z0=%g" % seed[2],
            "random=false", "seed=1", "scheme=" + scheme]
    run(folder, base + ["T=%g" % dt, "dump_intv=1000", "stat_intv=1000",
                        "checkpoint_intv=%g" % dt])
    put_points(folder, points)
    ck = checkpoint_folder(folder)
    stdout = run(folder, base + args + ["checkpoint_intv=1e9", "restart_folder=" + str(ck)])
    if out is not None:
        out.append(stdout)
    return ck


def probe(case, mode, points):
    """The velocity at each point, in the order given, with wall_upgrade=mode."""
    folder = restart_from(case, mode, points,
                          ["T=%g" % (2 * DT), "dump_intv=%g" % (2 * DT), "stat_intv=1000"])
    ids = read_checkpoint(folder)["id"][:, 0].astype(int)
    g = dump_at(folder, 2 * DT)
    # dump_at orders by id; the probe points were written in the checkpoint's id order
    dim = points.shape[1]
    assert np.allclose(g["points"][ids][:, :dim], points, rtol=0, atol=1e-12)
    return g["u"][ids][:, :dim]


def inside(geometry, points):
    """Which points lie in a mesh cell."""
    X, cells = geometry
    found = np.zeros(len(points), dtype=bool)
    for c in cells:
        T = (X[c[1:]] - X[c[0]]).T
        lam = np.linalg.solve(T, (points - X[c[0]]).T)
        found |= (lam >= 0).all(axis=0) & (lam.sum(axis=0) <= 1)
    return found
