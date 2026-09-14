"""The mesh interpolators (mode=triangle, tet, fenics, xdmftriangle, xdmftet).

data_example ships generate_up.py but not the mesh.h5 it writes, so the meshes
and fields are generated at test time by conftest (and, for XDMF, by the
xdmf_steps fixture below). Generating them needs dolfin, so the whole module
skips on a build without it. On both meshes |u| <= sqrt(3); the triangle case
is plane Poiseuille with u along y depending only on x, so x is conserved
exactly along every trajectory.
"""

import os
import shutil
import subprocess

import numpy as np
import pytest

from dumps import by_id
from paths import REPO, app, built_with_dolfin

PARTRAC = app("partrac")
DATA = os.path.join(REPO, "data_example")

pytestmark = pytest.mark.skipif(not built_with_dolfin(),
                                reason="partrac was built without dolfin")

# mode -> the mesh kind that conftest generates
CASES = [("triangle", "triangle"), ("tet", "tet")]

BASE = ("init_mode=uniform_x Nrw=200 Nrw_max=5000 ds_max=0.4 ds_min=0.1 "
        "Dm=0 int_order=1 dt=0.005 T=0.05 dump_intv=0.05 stat_intv=0.05 "
        "checkpoint_intv=0.05 random=false seed=1").split()

FILES = ("dolfin_params.dat", "mesh.h5", "up_0.h5", "timestamps.dat")


@pytest.fixture(params=CASES, ids=[c[0] for c in CASES])
def mesh_case(request, mesh_dir):
    """(mode, directory of the generated mesh case) for each mesh kind."""
    mode, kind = request.param
    return mode, mesh_dir(kind)


def run(case_dir, args):
    """Run partrac on case_dir/dolfin_params.dat with args; returns the process result."""
    r = subprocess.run([PARTRAC, str(case_dir / "dolfin_params.dat")] + args,
                       capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stdout + r.stderr
    return r


def make_run(mesh_case, dest, extra=()):
    """Copy the mesh case into dest and run BASE on it; returns dest."""
    mode, src = mesh_case
    dest.mkdir()
    for f in FILES:
        shutil.copy(src / f, dest / f)
    run(dest, BASE + ["mode=" + mode] + list(extra))
    return dest


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_particles_move_and_stay_in_the_periodic_domain(mesh_case, tmp_path):
    """Particles advect on the triangle and tet meshes, stay finite, and keep
    unwrapped positions across the periodic boundaries. A periodic run has to
    keep the true displacement, or dispersion computed from the dumps is wrong."""
    d = make_run(mesh_case, tmp_path / "run")
    h5py = pytest.importorskip("h5py")
    dumps = sorted(d.glob("**/data_from_t*.h5"))
    assert dumps
    with h5py.File(dumps[0]) as h:
        keys = sorted(h.keys(), key=float)
        first = by_id(h[keys[0]])["points"]
        last = by_id(h[keys[-1]])["points"]
    assert np.abs(last - first).max() > 1e-3
    assert np.isfinite(last).all()
    # |u| <= sqrt(3), so over T = 0.05 nothing can travel 0.2; a wrap into the
    # periodic cell would show as a jump of a whole period
    assert np.abs(last - first).max() < 0.2


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_mesh_run_is_reproducible(mesh_case, tmp_path):
    """Two identical runs on a mesh give bit-identical checkpoints, so mesh
    runs can be repeated and compared exactly."""
    a = make_run(mesh_case, tmp_path / "a")
    b = make_run(mesh_case, tmp_path / "b")
    xa = np.loadtxt(sorted(a.glob("**/Checkpoints/positions.pos"))[0])
    xb = np.loadtxt(sorted(b.glob("**/Checkpoints/positions.pos"))[0])
    assert np.array_equal(xa, xb)


@pytest.fixture(scope="session")
def xdmf_steps(tmp_path_factory):
    """xdmf_steps(kind, n) -> a dolfin-written XDMF case on a triangle or tet mesh
    with n timesteps (t = 0, then t = 1), u = (0, sin 2 pi x[, 0]) and p = 0."""
    df = pytest.importorskip("dolfin", reason="writing XDMF needs dolfin")
    built = {}

    def get(kind, n):
        if (kind, n) not in built:
            d = tmp_path_factory.mktemp("xdmf_%s_%d" % (kind, n))
            if kind == "triangle":
                mesh, uexpr = df.UnitSquareMesh(8, 8), ("0.0", "sin(2*M_PI*x[0])")
            else:
                mesh, uexpr = df.UnitCubeMesh(4, 4, 4), ("0.0", "sin(2*M_PI*x[0])", "0.0")
            fields = {
                "u": df.interpolate(df.Expression(uexpr, degree=1),
                                    df.VectorFunctionSpace(mesh, "CG", 1)),
                "p": df.interpolate(df.Expression("0.0", degree=1),
                                    df.FunctionSpace(mesh, "CG", 1)),
            }
            cwd = os.getcwd()
            os.chdir(d)
            try:
                for name, f in fields.items():
                    xf = df.XDMFFile(name + ".xdmf")
                    xf.parameters["functions_share_mesh"] = True
                    xf.parameters["rewrite_function_mesh"] = False
                    for t in (0.0, 1.0)[:n]:
                        xf.write(f, t)
                    xf.close()
            finally:
                os.chdir(cwd)
            (d / "dolfin_params.dat").write_text(
                "u=u.xdmf\np=p.xdmf\n" + "".join("periodic_%s=false\n" % a for a in "xyz"))
            built[(kind, n)] = d
        return built[(kind, n)]

    return get


# mode -> the mesh kind it runs on
DEGENERATE = {"triangle": "triangle", "fenics": "triangle", "tet": "tet",
              "xdmftriangle": "triangle", "xdmftet": "tet"}


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
@pytest.mark.parametrize("mode", list(DEGENERATE))
@pytest.mark.parametrize("stamps", ["last", "single"])
def test_fields_hold_on_a_degenerate_stamp_bracket(mesh_dir, xdmf_steps, tmp_path, mode, stamps):
    """On the last timestamp, or when there is only one, both ends of the time
    bracket are the same stamp, and every interpolator must return that stamp's
    field, finite everywhere. A 0/0 time weight there gives NaN velocities and
    accelerations, and particles can no longer take a step at the end of a run
    or in a steady single-snapshot case."""
    h5py = pytest.importorskip("h5py")
    kind = DEGENERATE[mode]
    d = tmp_path / "run"
    if mode.startswith("xdmf"):
        shutil.copytree(xdmf_steps(kind, 1 if stamps == "single" else 2), d)
        t_last = 1.
    else:
        src = mesh_dir(kind)
        d.mkdir()
        for f in FILES:
            shutil.copy(src / f, d / f)
        t_last = max(float(l.split()[0]) for l in
                     (d / "timestamps.dat").read_text().splitlines() if l.strip())
        if stamps == "single":
            (d / "timestamps.dat").write_text("0 up_0.h5\n")
    # dt = 1/8 is exact in binary, so the run lands on the stamp exactly
    if stamps == "single":
        extra, t_end = ["t0=0", "T=1"], 0.
    else:
        extra, t_end = ["t0=%r" % (t_last - 0.25), "T=%r" % t_last], t_last
    args = {a.split("=")[0]: a for a in
            BASE + ["mode=" + mode, "int_order=2", "dt=0.125", "dump_intv=0.125"] + extra}
    run(d, list(args.values()))

    with h5py.File(sorted(d.glob("**/data_from_t*.h5"))[-1]) as h:
        key = max(h.keys(), key=float)
        fields = {k: np.array(h[key][k]) for k in h[key]}
    assert float(key) == t_end
    for name, a in fields.items():
        assert np.isfinite(a).all(), name
    # the checkpoint is written after the step from the last dump, so particles
    # that could still step have moved away from the dumped positions
    x_end = np.loadtxt(sorted(d.glob("**/Checkpoints/positions.pos"))[0])
    assert np.abs(x_end - fields["points"]).max() > 1e-6


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_plane_poiseuille_on_a_mesh_conserves_x(mesh_case, tmp_path):
    """On the triangle mesh the field is plane Poiseuille, u along y depending
    only on x, so x must stay fixed to round-off while y advances. Drift in x
    would mean the mesh interpolation does not reproduce the field's direction."""
    if mesh_case[0] != "triangle":
        pytest.skip("only the plane Poiseuille case has an invariant axis")
    d = make_run(mesh_case, tmp_path / "axis")
    h5py = pytest.importorskip("h5py")
    with h5py.File(sorted(d.glob("**/data_from_t*.h5"))[0]) as h:
        keys = sorted(h.keys(), key=float)
        first = by_id(h[keys[0]])["points"]
        last = by_id(h[keys[-1]])["points"]
    assert np.abs(last[:, 0] - first[:, 0]).max() < 1e-12
    assert np.abs(last[:, 1] - first[:, 1]).max() > 1e-3
