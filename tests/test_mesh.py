"""The mesh interpolators (mode=triangle, tet, fenics, xdmftriangle, xdmftet).

data_example ships generate_up.py but not the mesh.h5 it writes, so the meshes
and fields, and the XDMF cases, are generated at test time by conftest.
Generating them needs python dolfin; reading them does not, except in
mode=fenics. On both meshes |u| <= sqrt(3); the triangle case
is plane Poiseuille with u along y depending only on x, so x is conserved
exactly along every trajectory.
"""

import os
import shutil

import numpy as np
import pytest

from dumps import all_dumps, dump_at
from paths import app, built_with_dolfin
from runs import copy_case, read_checkpoint, run_app

PARTRAC = app("partrac")
INTERPOL = app("interpol")

# Every mode here but fenics is read from the file's arrays, so only the
# fenics cases need a build with dolfin; writing the fixtures needs python
# dolfin, which conftest's mesh_dir asks for.
needs_fenics = [pytest.mark.fenics,
                pytest.mark.skipif(not built_with_dolfin(),
                                   reason="mode=fenics needs a build with dolfin")]

# mode -> the mesh kind that conftest generates
CASES = [("triangle", "triangle"), ("tet", "tet")]

BASE = ("init_mode=uniform_x Nrw=200 Nrw_max=5000 ds_max=0.4 ds_min=0.1 "
        "Dm=0 int_order=1 dt=0.005 T=0.05 dump_intv=0.05 stat_intv=0.05 "
        "checkpoint_intv=0.05 random=false seed=1").split()



@pytest.fixture(params=CASES, ids=[c[0] for c in CASES])
def mesh_case(request, mesh_dir):
    """(mode, directory of the generated mesh case) for each mesh kind."""
    mode, kind = request.param
    return mode, mesh_dir(kind)


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_particles_move_and_stay_in_the_periodic_domain(mesh_case, tmp_path):
    """Particles advect on the triangle and tet meshes, stay finite, and keep
    unwrapped positions across the periodic boundaries. A periodic run has to
    keep the true displacement, or dispersion computed from the dumps is wrong.

    On the triangle mesh the field is plane Poiseuille, u along y depending
    only on x, so x must also stay fixed to round-off while y advances. Drift
    in x would mean the mesh interpolation does not reproduce the field's
    direction.
    """
    mode, src = mesh_case
    d = copy_case(src, tmp_path / "run")
    run_app(PARTRAC, d / "dolfin_params.dat", BASE, "mode=" + mode)
    pytest.importorskip("h5py")
    dumps = all_dumps(d)
    first, last = dumps[min(dumps)]["points"], dumps[max(dumps)]["points"]
    assert np.abs(last - first).max() > 1e-3
    assert np.isfinite(last).all()
    # |u| <= sqrt(3), so over T = 0.05 nothing can travel 0.2; a wrap into the
    # periodic cell would show as a jump of a whole period
    assert np.abs(last - first).max() < 0.2
    if mode == "triangle":
        assert np.abs(last[:, 0] - first[:, 0]).max() < 1e-12
        assert np.abs(last[:, 1] - first[:, 1]).max() > 1e-3


# mode -> the mesh kind it runs on
DEGENERATE = {"triangle": "triangle", "fenics": "triangle", "tet": "tet",
              "xdmftriangle": "triangle", "xdmftet": "tet"}
DEGENERATE_PARAMS = [pytest.param(m, marks=needs_fenics) if m == "fenics" else m
                     for m in DEGENERATE]


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
@pytest.mark.parametrize("mode", DEGENERATE_PARAMS)
@pytest.mark.parametrize("stamps", ["last", "single"])
def test_fields_hold_on_a_degenerate_stamp_bracket(mesh_dir, xdmf, tmp_path, mode, stamps):
    """On the last timestamp, or when there is only one, both ends of the time
    bracket are the same stamp, and every interpolator must return that stamp's
    field, finite everywhere. A 0/0 time weight there gives NaN velocities and
    accelerations, and particles can no longer take a step at the end of a run
    or in a steady single-snapshot case."""
    pytest.importorskip("h5py")
    kind = DEGENERATE[mode]
    d = tmp_path / "run"
    if mode.startswith("xdmf"):
        copy_case(xdmf(kind, 1 if stamps == "single" else 2), d)
        t_last = 1.
    else:
        copy_case(mesh_dir(kind), d)
        t_last = max(float(l.split()[0]) for l in
                     (d / "timestamps.dat").read_text().splitlines() if l.strip())
        if stamps == "single":
            (d / "timestamps.dat").write_text("0 up_0.h5\n")
    # dt = 1/8 is exact in binary, so the run lands on the stamp exactly
    if stamps == "single":
        extra, t_end = ["t0=0", "T=1"], 0.
    else:
        extra, t_end = ["t0=%r" % (t_last - 0.25), "T=%r" % t_last], t_last
    run_app(PARTRAC, d / "dolfin_params.dat", BASE,
            ["mode=" + mode, "int_order=2", "dt=0.125", "dump_intv=0.125"], extra)

    dumps = all_dumps(d, raw=True)
    t = max(dumps)
    fields = dumps[t]
    assert t == t_end
    for name, a in fields.items():
        assert np.isfinite(a).all(), name
    # the checkpoint is written after the step from the last dump, so particles
    # that could still step have moved away from the dumped positions
    x_end = read_checkpoint(d)["points"]
    assert np.abs(x_end - fields["points"]).max() > 1e-6


@pytest.mark.fenics
@pytest.mark.skipif(not built_with_dolfin(), reason="mode=fenics needs a build with dolfin")
def test_fenics_ignores_the_pressure_when_asked(mesh_dir, tmp_path):
    """ignore_pressure=true lets mode=fenics run on a stamp file without a
    pressure field and dump p = 0; without it the missing field stops the run.
    ignore_pressure is a key of every mesh loader's file, and mode=fenics must
    honour it as the others do. (The shipped test_tet_p2 case sets it.)"""
    h5py = pytest.importorskip("h5py")
    src = mesh_dir("tet")
    for name, ignore in (("ignored", "true"), ("read", "false")):
        d = copy_case(src, tmp_path / name)
        with h5py.File(d / "up_0.h5", "a") as h:
            del h["p"]
        params = [l for l in (d / "dolfin_params.dat").read_text().splitlines()
                  if not l.startswith("ignore_pressure")]
        (d / "dolfin_params.dat").write_text("\n".join(params + ["ignore_pressure=" + ignore]) + "\n")
        r = run_app(PARTRAC, d / "dolfin_params.dat", BASE, "mode=fenics", check=False)
        if ignore == "false":
            assert r.returncode != 0, "a missing pressure field went unnoticed"
            continue
        assert r.returncode == 0, r.stdout + r.stderr
        ps = [g["p"] for g in all_dumps(d, raw=True).values() if "p" in g]
        assert ps and all(not p.any() for p in ps)


@pytest.mark.fenics
@pytest.mark.skipif(not built_with_dolfin(), reason="mode=fenics needs a build with dolfin")
def test_fenics_blends_two_stamp_files_as_the_tet_loader_does(mesh_dir, tmp_path):
    """Between two stamp files mode=fenics reads both and blends them in time,
    and in the next bracket keeps the later one and reads the one after: a run
    across three stamps, the middle one reversed and halved, traces what
    mode=tet traces, to round-off."""
    h5py = pytest.importorskip("h5py")
    out = {}
    for mode in ("fenics", "tet"):
        d = copy_case(mesh_dir("tet"), tmp_path / mode)
        shutil.copy(d / "up_0.h5", d / "up_1.h5")
        with h5py.File(d / "up_1.h5", "r+") as h:
            h["u/vector_0"][...] *= -0.5
        (d / "timestamps.dat").write_text("0 up_0.h5\n0.25 up_1.h5\n0.5 up_0.h5\n")
        run_app(app("tracers"), d / "dolfin_params.dat",
                "mode=%s init_mode=points_xyz x0=0.5 y0=0.5 z0=0.5 Nrw=50 Nrw_max=50 Dm=0 int_order=1 "
                "scheme=RK4 dt=0.025 T=0.5 dump_intv=0.5 stat_intv=1e9 checkpoint_intv=1e9 random=false "
                "seed=1" % mode)
        out[mode] = (dump_at(d, 0.0)["points"], dump_at(d, 0.5)["points"])
    (x0, xf), (_, xt) = out["fenics"], out["tet"]
    assert np.abs(xf - x0).max() > 1e-2
    assert np.abs(xf - xt).max() < 1e-12


@pytest.fixture
def fenics_case(tmp_path):
    """(dolfin, the folder to write a case in) for interpol's mode=fenics; a
    skip without a build with dolfin, python dolfin or interpol."""
    if not built_with_dolfin():
        pytest.skip("mode=fenics needs a build with dolfin")
    df = pytest.importorskip("dolfin", reason="writing the case needs dolfin")
    if not os.path.exists(INTERPOL):
        pytest.skip("interpol is not built")
    return df, tmp_path / "case"


@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.fenics
def test_fenics_reads_a_p3_p2_field_exactly(fenics_case, dim):
    """mode=fenics takes P1 to P3, and nothing else ran the P3 velocity or the
    P2 pressure spaces. A cubic velocity and a quadratic pressure are held
    exactly by P3-P2, so the probed values must match them to round-off; a
    wrong dof order or a wrong basis would not."""
    df, d = fenics_case
    h5py = pytest.importorskip("h5py")
    d.mkdir()
    mesh = df.UnitSquareMesh(3, 3) if dim == 2 else df.UnitCubeMesh(2, 2, 2)
    V = df.VectorFunctionSpace(mesh, "CG", 3)
    P = df.FunctionSpace(mesh, "CG", 2)
    u_expr = ["x[1]*x[1]*x[1] + 0.5*x[0]", "x[0]*x[0] - x[1]"] + (["0.25*x[2]*x[0]"] if dim == 3 else [])
    u = df.interpolate(df.Expression(u_expr, degree=3), V)
    p = df.interpolate(df.Expression("x[0]*x[1] + 1.0", degree=2), P)
    with df.HDF5File(mesh.mpi_comm(), str(d / "mesh.h5"), "w") as f:
        f.write(mesh, "mesh")
    with df.HDF5File(mesh.mpi_comm(), str(d / "up_0.h5"), "w") as f:
        f.write(u, "u")
        f.write(p, "p")
    (d / "timestamps.dat").write_text("0.0\tup_0.h5\n")
    (d / "dolfin_params.dat").write_text(
        "velocity_space=P3\npressure_space=P2\ntimestamps=timestamps.dat\nmesh=mesh.h5\n"
        "periodic_x=false\nperiodic_y=false\nperiodic_z=false\nrho=1.0\n")
    run_app(INTERPOL, d / "dolfin_params.dat",
            "mode=fenics Nrw=400 int_order=2 t0=0 random=false seed=1", timeout=600)
    out = list(d.rglob("interpolation.h5part"))
    assert len(out) == 1
    with h5py.File(out[0], "r") as h:
        v = {k: np.array(h["Step#0"][k]) for k in h["Step#0"]}
    x, y = v["x"], v["y"]
    z = v["z"] if dim == 3 else np.zeros_like(x)
    assert len(x) > 100
    assert np.abs(v["ux"] - (y**3 + 0.5 * x)).max() < 1e-10
    assert np.abs(v["uy"] - (x**2 - y)).max() < 1e-10
    if dim == 3:
        assert np.abs(v["uz"] - 0.25 * z * x).max() < 1e-10
    assert np.abs(v["p"] - (x * y + 1.0)).max() < 1e-10
    # the gradient of a cubic is held too: du_x/dy = 3 y^2
    assert np.abs(v["uxy"] - 3 * y**2).max() < 1e-9


def periodic_unit_case(df, d, dim, degree, shuffled, constrained, offset=0.0):
    """A unit square or cube periodic in every direction, with a smooth
    periodic velocity and pressure in P<degree> written to d; the vertices as
    dolfin generates them or shuffled (reversing seam edges against their
    images), the fields interpolated with dolfin's periodic constraint or
    without (constrained: both, or a pair for the velocity and the pressure).
    A nonzero offset is added to the pressure and is the velocity's first
    component. Returns the written velocity and pressure."""
    mesh = df.UnitSquareMesh(4, 4) if dim == 2 else df.UnitCubeMesh(2, 2, 2)
    if shuffled:
        X, cells = mesh.coordinates(), mesh.cells()
        perm = np.random.default_rng(1).permutation(len(X))
        mesh = df.Mesh()
        ed = df.MeshEditor()
        ed.open(mesh, "triangle" if dim == 2 else "tetrahedron", dim, dim)
        ed.init_vertices(len(X))
        ed.init_cells(len(cells))
        for i, x in enumerate(X[perm]):
            ed.add_vertex(i, x)
        for i, c in enumerate(np.argsort(perm)[cells]):
            ed.add_cell(i, c)
        ed.close()

    class Periodic(df.SubDomain):
        def inside(self, x, on):
            return bool(on and any(df.near(x[k], 0) for k in range(dim))
                        and not any(df.near(x[k], 1) for k in range(dim)))

        def map(self, x, y):
            for k in range(dim):
                y[k] = x[k] - 1 if df.near(x[k], 1) else x[k]

    cu, cp = constrained if isinstance(constrained, tuple) else (constrained, constrained)
    V = df.VectorFunctionSpace(mesh, "CG", degree, **({"constrained_domain": Periodic()} if cu else {}))
    Q = df.FunctionSpace(mesh, "CG", degree, **({"constrained_domain": Periodic()} if cp else {}))
    u_expr = ["sin(2*pi*x[1]) + 0.3*cos(2*pi*x[0])", "cos(2*pi*x[0])*sin(2*pi*x[1])"]
    p_expr = "cos(2*pi*x[0])*sin(2*pi*x[1])"
    if dim == 3:
        u_expr = [e + " + 0.2*sin(2*pi*x[2])" for e in u_expr] + ["cos(2*pi*x[0] + 2*pi*x[2])"]
        p_expr += " + cos(2*pi*x[2])"
    if offset:
        u_expr[0] = "%r" % offset
        p_expr = "%r + %s" % (offset, p_expr)
    u = df.interpolate(df.Expression(u_expr, degree=degree + 2), V)
    p = df.interpolate(df.Expression(p_expr, degree=degree + 2), Q)
    d.mkdir()
    with df.HDF5File(mesh.mpi_comm(), str(d / "mesh.h5"), "w") as f:
        f.write(mesh, "mesh")
    with df.HDF5File(mesh.mpi_comm(), str(d / "up_0.h5"), "w") as f:
        f.write(u, "u")
        f.write(p, "p")
    (d / "timestamps.dat").write_text("0.0\tup_0.h5\n")
    (d / "dolfin_params.dat").write_text(
        "velocity_space=P%d\npressure_space=P%d\ntimestamps=timestamps.dat\nmesh=mesh.h5\n"
        "periodic_x=true\nperiodic_y=true\nperiodic_z=%s\nrho=1.0\n"
        % (degree, degree, "true" if dim == 3 else "false"))
    return u, p


def seam_jump(df, f, dim):
    """The largest difference of f between seam points and their images."""
    s = np.linspace(0.013, 0.987, 17)
    pts = np.stack(np.meshgrid(*[s] * (dim - 1), indexing="ij"), -1).reshape(-1, dim - 1)
    jump = 0.
    for k in range(dim):
        for q in pts:
            a, b = np.insert(q, k, 0.), np.insert(q, k, 1.)
            jump = max(jump, np.abs(f(df.Point(*a)) - f(df.Point(*b))).max())
    return jump


def probe_fenics(d, args=()):
    """The run of interpol in mode=fenics on the case in d, and its probes."""
    r = run_app(INTERPOL, d / "dolfin_params.dat",
                "mode=fenics Nrw=1500 int_order=2 t0=0 random=false seed=1", list(args),
                check=False, timeout=600)
    if r.returncode != 0:
        return r, None
    h5py = pytest.importorskip("h5py")
    out = list(d.rglob("interpolation.h5part"))
    assert len(out) == 1
    with h5py.File(out[0], "r") as h:
        return r, {k: np.array(h["Step#0"][k]) for k in h["Step#0"]}


def assert_probes_match(df, v, u, p, dim):
    """The probed velocity and pressure are u's and p's at the probed points."""
    xs = np.stack([v[a] for a in "xyz"[:dim]], -1)
    assert len(xs) > 1000
    pts = [df.Point(*x) for x in xs]
    uu = np.array([u(q) for q in pts])
    pp = np.array([p(q) for q in pts])
    got = np.stack([v["u" + a] for a in "xyz"[:dim]], -1)
    assert np.abs(got - uu).max() < 1e-10
    assert np.abs(v["p"] - pp).max() < 1e-10


@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.fenics
def test_fenics_reads_a_periodic_p3_field_on_reversed_seam_edges(fenics_case, dim):
    """dolfin's periodic constraint pairs the two dofs of a P3 seam edge with
    its image's by each edge's own vertex order, so on a numbering that
    reverses an edge against its image it ties the 1/3 node of one to the 2/3
    node of the other. A periodic field written without that constraint is
    good, and mode=fenics must read it as written: the same at every probe as
    dolfin's own evaluation of it."""
    df, d = fenics_case
    u, p = periodic_unit_case(df, d, dim, 3, shuffled=True, constrained=False)
    assert seam_jump(df, u, dim) < 1e-12
    r, v = probe_fenics(d)
    assert r.returncode == 0, r.stdout + r.stderr
    assert_probes_match(df, v, u, p, dim)


@pytest.mark.fenics
def test_fenics_reads_a_good_periodic_p3_field_offset_by_1e7(fenics_case):
    """A good periodic P3 field beside values of 1e7 passes the seam check:
    images apart by the values' round-off (1e-8 here) are no jump."""
    df, d = fenics_case
    h5py = pytest.importorskip("h5py")
    periodic_unit_case(df, d, 2, 3, shuffled=True, constrained=False, offset=1e7)
    rng = np.random.default_rng(1)
    with h5py.File(d / "up_0.h5", "r+") as f:
        for name in ("u/vector_0", "p/vector_0"):
            f[name][...] += rng.uniform(-1e-8, 1e-8, f[name].shape)
    r, _ = probe_fenics(d)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "periodic image" not in r.stdout + r.stderr


@pytest.mark.fenics
def test_fenics_reads_a_periodic_p3_field_named_by_its_dataset(fenics_case):
    """A field may be named by its vector dataset, as dolfin's read takes it
    (u/vector_0): read by cells all the same."""
    df, d = fenics_case
    u, p = periodic_unit_case(df, d, 2, 3, shuffled=True, constrained=False)
    with open(d / "dolfin_params.dat", "a") as f:
        f.write("velocity_field=u/vector_0\npressure_field=p/vector_0\n")
    r, v = probe_fenics(d)
    assert r.returncode == 0, r.stdout + r.stderr
    assert_probes_match(df, v, u, p, 2)


@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.fenics
def test_fenics_refuses_a_periodic_p3_field_wrong_at_the_seam(fenics_case, dim):
    """The same field interpolated into dolfin's periodic P3 space on that
    numbering is wrong at the seam, and so is any periodic P3 solve there:
    the file differs between a seam point and its image. mode=fenics must
    refuse it, naming the cause, not trace it."""
    df, d = fenics_case
    u, _ = periodic_unit_case(df, d, dim, 3, shuffled=True, constrained=True)
    assert seam_jump(df, u, dim) > 1e-2
    r, _ = probe_fenics(d)
    assert r.returncode != 0, "a P3 field wrong at the periodic seam went unnoticed"
    assert "periodic image" in r.stdout + r.stderr


@pytest.mark.parametrize("field", ["velocity", "pressure"])
@pytest.mark.fenics
def test_fenics_refuses_a_seam_jump_small_beside_the_values(fenics_case, field):
    """The seam check measures a jump against the range of each component, not
    its largest value: a pressure offset by 1e7, or a velocity whose other
    component is 1e7, hides nothing. The wrong field is the one in the
    periodic space, the other written good."""
    df, d = fenics_case
    constrained = (field == "velocity", field == "pressure")
    u, p = periodic_unit_case(df, d, 2, 3, shuffled=True, constrained=constrained, offset=1e7)
    wrong = u if field == "velocity" else p
    assert seam_jump(df, wrong, 2) > 1e-2
    r, _ = probe_fenics(d)
    assert r.returncode != 0, "a seam jump beside values of 1e7 went unnoticed"
    assert "the %s differs between periodic images" % field in r.stdout + r.stderr


@pytest.mark.parametrize("dim,degree,shuffled",
                         [(2, 1, False), (2, 1, True), (2, 2, False), (2, 2, True),
                          (2, 3, False), (3, 2, True), (3, 3, False)])
@pytest.mark.fenics
def test_fenics_reads_good_periodic_fields(fenics_case, dim, degree, shuffled):
    """A periodic field in dolfin's periodic space is good in P1 and P2 on any
    numbering, and in P3 on one that keeps seam edges' orientation: read as
    written, and the seam check stays quiet."""
    df, d = fenics_case
    u, p = periodic_unit_case(df, d, dim, degree, shuffled=shuffled, constrained=True)
    assert seam_jump(df, u, dim) < 1e-12
    r, v = probe_fenics(d)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "periodic image" not in r.stdout + r.stderr
    assert_probes_match(df, v, u, p, dim)
