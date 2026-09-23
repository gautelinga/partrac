"""The mesh interpolators (mode=triangle, tet, fenics, xdmftriangle, xdmftet).

data_example ships generate_up.py but not the mesh.h5 it writes, so the meshes
and fields, and the XDMF cases, are generated at test time by conftest.
Generating them needs python dolfin; reading them does not, except in
mode=fenics. On both meshes |u| <= sqrt(3); the triangle case
is plane Poiseuille with u along y depending only on x, so x is conserved
exactly along every trajectory.
"""

import os

import numpy as np
import pytest

from dumps import all_dumps
from paths import app, built_with_dolfin
from runs import copy_case, run_app

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
    x_end = np.loadtxt(sorted(d.glob("**/Checkpoints/positions.pos"))[0])
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


@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.fenics
@pytest.mark.skipif(not built_with_dolfin(), reason="mode=fenics needs a build with dolfin")
def test_fenics_reads_a_p3_p2_field_exactly(tmp_path, dim):
    """mode=fenics takes P1 to P3, and nothing else ran the P3 velocity or the
    P2 pressure spaces. A cubic velocity and a quadratic pressure are held
    exactly by P3-P2, so the probed values must match them to round-off; a
    wrong dof order or a wrong basis would not."""
    df = pytest.importorskip("dolfin", reason="writing the case needs dolfin")
    h5py = pytest.importorskip("h5py")
    if not os.path.exists(INTERPOL):
        pytest.skip("interpol is not built")
    d = tmp_path / "p3"
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
