"""Shared fixtures that build the non-analytic inputs the apps read.

The mesh examples ship a generator rather than the mesh, so the cases that need
a mesh build it here, once per session and shared by the xdist workers, and
skip when dolfin is missing. The FELBM and XDMF cases are small synthetic inputs
written directly in the formats StructuredInterpol and the XDMF loaders read.
"""

import os
import shutil
import subprocess

import pytest

from cases import shared_dir, write_felbm
from paths import REPO

DATA = os.path.join(REPO, "data_example")

# Under pytest-xdist the workers share the cores: an app's OpenMP threads,
# unless a test sets them, get the machine divided by the workers. Spinning
# OpenMP threads on an oversubscribed machine make the suite slower than serial.
if "PYTEST_XDIST_WORKER_COUNT" in os.environ:
    os.environ.setdefault("OMP_NUM_THREADS", str(max(1, (os.cpu_count() or 1)
                                                     // int(os.environ["PYTEST_XDIST_WORKER_COUNT"]))))


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: left out of the quick run, pytest -m 'not slow'")
    # mode=fenics is the only mode a build without dolfin does not have; a run
    # against such a build names what it leaves out, pytest -m 'not fenics'
    config.addinivalue_line("markers", "fenics: needs mode=fenics, so a build with dolfin")


# kind -> (example folder, generator arguments)
MESH_KINDS = {
    "triangle": ("ppf_triangle_p2", ["-dim", "1"]),
    "tet": ("test_tet_p2", []),
    "trianglefreq": ("sine_trianglefreq_p2", []),
    "tetfreq": ("sine_tetfreq_p2", []),
}


@pytest.fixture(scope="session")
def mesh_dir(tmp_path_factory):
    """mesh_dir(kind) -> a folder holding dolfin_params.dat and its generated mesh.

    Each mesh is generated once per session, by whichever worker asks first;
    the folder is shared, so a test copies what it runs on.
    """
    pytest.importorskip("dolfin", reason="mesh generation needs dolfin")

    def get(kind):
        folder, args = MESH_KINDS[kind]

        def build(d):
            for f in os.listdir(os.path.join(DATA, folder)):
                shutil.copy(os.path.join(DATA, folder, f), d / f)
            r = subprocess.run(["python3", "generate_up.py"] + args,
                               cwd=d, capture_output=True, text=True, timeout=900)
            assert r.returncode == 0, r.stdout + r.stderr
            assert (d / "mesh.h5").exists(), "the generator wrote no mesh"

        return shared_dir(tmp_path_factory, folder, build)

    return get


@pytest.fixture(scope="session")
def felbm_dir(tmp_path_factory):
    """A folder holding a synthetic FELBM case: felbm_params.dat and its h5 files.

    The FELBM solver is a separate project and its output is not shipped, but
    the format is small. felbm_params.dat names a timestamps file and an
    is_solid file; is_solid fixes the grid size, and each timestep holds u_x,
    u_y, u_z, density and pressure as doubles indexed nx*ny*iz + nx*iy + ix.
    The grid is cubic so that the index order is unambiguous.
    """
    np = pytest.importorskip("numpy")
    pytest.importorskip("h5py")

    def build(d):
        # 16^3 cell centres; u_y varies with x only, so the flow is a steady shear
        n = 16
        x = np.arange(n) + 0.5
        X = np.meshgrid(x, x, x, indexing="ij")[0]
        zero = np.zeros((n, n, n))
        fields = {"u_x": zero, "u_y": 0.05 * np.sin(2 * np.pi * X / n), "u_z": zero,
                  "density": np.ones((n, n, n)), "pressure": zero}
        # solid walls at both z ends, open elsewhere
        solid = np.zeros((n, n, n), dtype=np.int32)
        solid[0, :, :] = 1
        solid[-1, :, :] = 1
        # two identical timesteps, so the field is steady
        write_felbm(d, [fields, fields], solid)

    return shared_dir(tmp_path_factory, "felbm", build)


@pytest.fixture(scope="session")
def xdmf(tmp_path_factory):
    """xdmf(kind, n) -> a folder holding a dolfin-written XDMF case, the input the XDMF loaders read.

    dolfin_params.dat names one xdmf per field; each xdmf carries the topology
    and geometry paths into its h5 and one Grid per timestep. The mesh is the
    unit square (kind "triangle") or the unit cube ("tet"), the velocity
    u = (0, sin(2 pi x)[, 0]) and p = 0, written unchanged at t = 0 and, for
    n = 2, at t = 1.
    """
    df = pytest.importorskip("dolfin", reason="writing XDMF needs dolfin")

    def get(kind, n):
        def build(d):
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
                    # one mesh in the h5, referenced by every timestep's Grid
                    xf.parameters["functions_share_mesh"] = True
                    xf.parameters["rewrite_function_mesh"] = False
                    for t in (0.0, 1.0)[:n]:
                        xf.write(f, t)
                    xf.close()
            finally:
                os.chdir(cwd)
            (d / "dolfin_params.dat").write_text(
                "u=u.xdmf\np=p.xdmf\n" + "".join("periodic_%s=false\n" % a for a in "xyz"))

        return shared_dir(tmp_path_factory, "xdmf_%s_%d" % (kind, n), build)

    return get


@pytest.fixture(scope="session")
def xdmf_dir(xdmf):
    """The two-stamp XDMF case on the unit square."""
    return xdmf("triangle", 2)
