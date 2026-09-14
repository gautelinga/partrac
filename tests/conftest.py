"""Shared fixtures that build the non-analytic inputs the apps read.

The mesh examples ship a generator rather than the mesh, so the cases that need
a mesh build it here, once per session, and skip when dolfin is missing. The
FELBM and XDMF cases are small synthetic inputs written directly in the formats
StructuredInterpol and XDMFTriangleInterpol read.
"""

import os
import shutil
import subprocess

import pytest

from paths import REPO

DATA = os.path.join(REPO, "data_example")

# kind -> (example folder, generator arguments)
MESH_KINDS = {
    "triangle": ("ppf_triangle_p2", ["-dim", "1"]),
    "tet": ("test_tet_p2", []),
    "trianglefreq": ("sine_trianglefreq_p2", []),
}


@pytest.fixture(scope="session")
def mesh_dir(tmp_path_factory):
    """mesh_dir(kind) -> a folder holding dolfin_params.dat and its generated mesh."""
    pytest.importorskip("dolfin", reason="mesh generation needs dolfin")
    built = {}

    def get(kind):
        if kind not in built:
            folder, args = MESH_KINDS[kind]
            d = tmp_path_factory.mktemp(folder)
            for f in os.listdir(os.path.join(DATA, folder)):
                shutil.copy(os.path.join(DATA, folder, f), d / f)
            r = subprocess.run(["python3", "generate_up.py"] + args,
                               cwd=d, capture_output=True, text=True, timeout=900)
            assert r.returncode == 0, r.stdout + r.stderr
            assert (d / "mesh.h5").exists(), "the generator wrote no mesh"
            built[kind] = d
        return built[kind]

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
    h5py = pytest.importorskip("h5py")
    d = tmp_path_factory.mktemp("felbm")

    # 16^3 cell centres; u_y varies with x only, so the flow is a steady shear
    n = 16
    x = np.arange(n) + 0.5
    X = np.meshgrid(x, x, x, indexing="ij")[0]
    zero = np.zeros((n, n, n))
    fields = {"u_x": zero, "u_y": 0.05 * np.sin(2 * np.pi * X / n), "u_z": zero,
              "density": np.ones((n, n, n)), "pressure": zero}

    # solid walls at both x ends, open elsewhere
    solid = np.zeros((n, n, n), dtype=np.int32)
    solid[0, :, :] = 1
    solid[-1, :, :] = 1
    with h5py.File(d / "output_is_solid.h5", "w") as f:
        f.create_dataset("is_solid", data=solid)

    # two identical timesteps, so the field is steady
    for name in ("output_0.h5", "output_1.h5"):
        with h5py.File(d / name, "w") as f:
            for field, a in fields.items():
                f.create_dataset(field, data=np.transpose(a, (2, 1, 0)).astype(float))

    (d / "timestamps.dat").write_text("0\toutput_0.h5\n100\toutput_1.h5\n")
    (d / "felbm_params.dat").write_text(
        "timestamps=timestamps.dat\nis_solid_file=output_is_solid.h5\n")
    return d

@pytest.fixture(scope="session")
def xdmf_dir(tmp_path_factory):
    """A folder holding a dolfin-written XDMF case, the input XDMFTriangleInterpol reads.

    dolfin_params.dat names one xdmf per field; each xdmf carries the topology
    and geometry paths into its h5 and one Grid per timestep. The velocity is
    u = (0, sin(2 pi x)) on the unit square, written unchanged at t = 0 and 1.
    """
    df = pytest.importorskip("dolfin", reason="writing XDMF needs dolfin")
    d = tmp_path_factory.mktemp("xdmf")
    cwd = os.getcwd()
    os.chdir(d)
    try:
        mesh = df.UnitSquareMesh(8, 8)
        V = df.VectorFunctionSpace(mesh, "CG", 1)
        P = df.FunctionSpace(mesh, "CG", 1)
        fields = {
            "u": df.interpolate(df.Expression(("0.0", "sin(2*M_PI*x[0])"), degree=1), V),
            "p": df.interpolate(df.Expression("0.0", degree=1), P),
        }
        for name, f in fields.items():
            xf = df.XDMFFile(name + ".xdmf")
            # one mesh in the h5, referenced by every timestep's Grid
            xf.parameters["functions_share_mesh"] = True
            xf.parameters["rewrite_function_mesh"] = False
            for t in (0.0, 1.0):
                xf.write(f, t)
            xf.close()
        with open("dolfin_params.dat", "w") as f:
            f.write("u=u.xdmf\np=p.xdmf\nperiodic_x=false\nperiodic_y=false\n")
    finally:
        os.chdir(cwd)
    return d
