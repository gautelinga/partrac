"""mode=openfoam on the checked-in OpenFOAM cases, through the apps.

With the inverse-distance node values (nodes=inverse_distance, OpenFOAM's
volPointInterpolation) on the W12 split, which is cellPoint's own
decomposition, interpol must reproduce OpenFOAM's cellPoint at the points it
probes: cellpoint_U.txt beside each case holds OpenFOAM's values at
interpol's own 200 points (seed 1, one thread), so the whole path -- the
reader, the split, the node values, the tables, the periodic pairing, the
evaluation -- is checked to round-off against OpenFOAM itself. The 2D cavity
checks the triangle path and the lid corners (a moving and a fixed wall at
one point), the jittered channel the tet path and the cyclic pairs (no facet
left unpaired). With the default least-squares nodes the channel's field,
linear in z, must come out exact in value and gradient, which cellPoint's
does not; the cavity loads either way, zero winning at the lid's corners;
a stamp whose conditions differ from the first's is refused. The load
classes each boundary patch by its velocity condition (no-slip wall, moving
wall, cyclic, other) and says so; a condition built only from its dictionary
is read too. A 2D case may have its empty pair normal to any axis: the
cavity turned onto x and onto y is the same field and traces the same. A
steady case is traced frozen at one stamp: frozen_fields must give that
stamp's field at every step, through tracers and partrac. The refusals are
checked by their messages. A phase_field is P1 on the split, its nodes
cellPoint's, each mixed fanned cell's centre moved within [0, 1] so the
volume above 1/2 is the cells' volume fractions', with a gradient of its
own; a scalar condition OpenFOAM cannot build is not data, a velocity one is
refused. A tracers run on the channel is a smoke test, a case with cells
whose every point is on a no-slip wall is warned of, a build without the
reader refuses the mode, and, marked slow, an installed tree reads the cavity
with its own reader.

The tests after those run on solved OpenFOAM cases too large to check in:
the cyclic pairing complete on four of them; Poiseuille and Couette against
their closed forms; the pressure, and a symmetry plane, against OpenFOAM's
own cellPoint; and, marked slow, the tracers at rest near no-slip obstacles
in 2D and 3D with and without the near-wall rule, and diffusing tracers kept
out of the obstacles by the reflecting walk. PARTRAC_OPENFOAM_CASES lists
the folders that hold them (os.pathsep between), each case a solved
OpenFOAM case under its name: rsa2d, spheres3d, channel3d_jitter/n32,
channel3d_graded/n32, couette, cavity_p and halfchannel, with the
postProcessing probes the tests compare against. The make.sh scripts under
harness/openfoam/ make them. Without the variable, or without its case, a
test skips.

The reader is built with PARTRAC_ENABLE_OPENFOAM=ON; without it
(build_features.txt says openfoam=off) these tests skip, and fail instead
under PARTRAC_REQUIRE_OPENFOAM=1, as a job that has OpenFOAM runs them."""

import gzip
import os
import re
import shutil
import subprocess

import numpy as np
import pytest

from dumps import all_dumps
from paths import BIN, REPO, app
from runs import run_app

INTERPOL = app("interpol")
TRACERS = app("tracers")
CAVITY = os.path.join(REPO, "data_example", "openfoam_cavity")
CHANNEL = os.path.join(REPO, "data_example", "openfoam_channel3d")


def built_with_openfoam():
    features = os.path.join(BIN, "build_features.txt")
    return os.path.exists(features) and "openfoam=on" in open(features).read().split()


def _needs_openfoam():
    if built_with_openfoam():
        return pytest.mark.skipif(not os.path.exists(INTERPOL), reason="interpol is not built")
    if os.environ.get("PARTRAC_REQUIRE_OPENFOAM"):
        return pytest.mark.skipif(False, reason="")
    return pytest.mark.skip(reason="built without PARTRAC_ENABLE_OPENFOAM")


needs_openfoam = _needs_openfoam()
PROBE = ["mode=openfoam", "Nrw=200", "int_order=2", "random=false", "seed=1", "num_threads=1"]


def case_copy(src, d, params=None):
    """The case copied to d, its parameter file replaced by `params` lines if given."""
    shutil.copytree(src, d, ignore=shutil.ignore_patterns("generate_case.sh", "*.py", "split_expected.h5"))
    if params is not None:
        (d / "partrac_params.dat").write_text("".join(p + "\n" for p in params))
    return d / "partrac_params.dat"


def probe(params, *args, gradient=False):
    """interpol's points (x, y, z) and velocities, and its output, stdout
    then stderr; with gradient, the velocity gradients J[:, i, j] = du_i/dx_j
    after them."""
    import h5py
    r = run_app(INTERPOL, params, PROBE, list(args), env={"OMP_NUM_THREADS": "1"}, timeout=300)
    [out] = list(params.parent.rglob("interpolation.h5part"))
    with h5py.File(out, "r") as h:
        g = h[list(h.keys())[0]]
        X, U = np.c_[g["x"], g["y"], g["z"]], np.c_[g["ux"], g["uy"], g["uz"]]
        if gradient:
            J = np.stack([np.c_[g["u%sx" % a], g["u%sy" % a], g["u%sz" % a]] for a in "xyz"], 1)
            return X, U, J, r.stdout + r.stderr
        return X, U, r.stdout + r.stderr


def against_cellpoint(src, d, dim, *args):
    """The largest difference between interpol and OpenFOAM's cellPoint over
    the points both find inside, and how many they are; args after t0 go to
    the parameter file."""
    params = case_copy(src, d, ["velocity_field=U", "pressure_field=p", "wall_p2=none",
                                "nodes=inverse_distance"] + list(args[1:]))
    args = args[:1]
    X, U, log = probe(params, *args)
    ref = np.loadtxt(os.path.join(src, "cellpoint_U.txt"))
    inside = np.abs(ref[:, 3]) < 1e200
    # interpol writes its points in its own order
    order = np.lexsort(X[:, :dim].T[::-1])
    ro = np.lexsort(ref[:, :dim].T[::-1])
    assert np.abs(X[order, :dim] - ref[ro, :dim]).max() == 0.0
    diff = np.abs(U[order, :dim] - ref[ro, 3:3 + dim])[inside[ro]]
    return diff.max(), int(inside.sum()), log


@needs_openfoam
def test_the_cavity_is_cellpoint(tmp_path):
    """2D: the front plane fanned, the lid corners included (OpenFOAM's 1e-13
    at seven significant digits of a unit lid speed is round-off)."""
    worst, n, log = against_cellpoint(CAVITY, tmp_path / "c", 2, "t0=10")
    assert n == 200
    assert worst < 1e-13, worst


@needs_openfoam
def test_the_jittered_channel_is_cellpoint_and_every_cyclic_facet_is_paired(tmp_path):
    worst, n, log = against_cellpoint(CHANNEL, tmp_path / "c", 3, "t0=0")
    assert n == 200
    assert worst < 1e-13, worst
    assert "Split tets: 20736 for 1728 cells" in log
    assert "1152 cyclic (left right front back)" in log, log


@needs_openfoam
@pytest.mark.parametrize("src,dim,t0", [(CAVITY, 2, "10"), (CHANNEL, 3, "0")])
def test_w6_loads(tmp_path, src, dim, t0):
    """split=6: Dompierre's tets (two triangles a quad in 2D) on the same
    nodes, without the cell values: near cellPoint, not equal to it. Two of
    the jittered channel's hexes have face diagonals that meet at no corner,
    and are fanned from their centres."""
    worst, n, log = against_cellpoint(src, tmp_path / "c", dim, "t0=" + t0, "split=6")
    assert ("Split tets: 10380 for 1728 cells, 2199 nodes; 2 cells fanned" if dim == 3
            else "Split triangles: 800") in log, log
    assert 1e-6 < worst < (0.3 if dim == 2 else 0.05), worst


@needs_openfoam
def test_the_channel_is_exact_for_its_linear_field(tmp_path):
    """The default least-squares nodes are exact for a linear field on any
    mesh: the channel's u = z a + b (cyclic in x and y, fixed walls), jittered,
    comes out exact at every probed point, value and gradient, where cellPoint
    is off by the jitter."""
    params = case_copy(CHANNEL, tmp_path / "c")
    X, U, J, log = probe(params, "t0=0", gradient=True)
    a, b = np.array([2.0, 0.7, -1.1]), np.array([1.0, -0.5, 0.3])
    assert len(X) == 200
    assert "Nodes by least squares, U" in log, log
    assert np.abs(U - (X[:, 2:3] * a + b)).max() < 1e-13
    grad = np.zeros((3, 3))
    grad[:, 2] = a
    assert np.abs(J - grad).max() < 1e-11
    ref = np.loadtxt(os.path.join(CHANNEL, "cellpoint_U.txt"))
    assert np.abs(ref[:, 3:6] - (ref[:, 2:3] * a + b)).max() > 1e-4


@needs_openfoam
@pytest.mark.parametrize("nodes", ["least_squares", "inverse_distance"])
def test_the_cavity_loads_with_either_node_rule(tmp_path, nodes):
    """Least squares fixes the lid's two corners at rest (zero wins against
    the no-slip walls) and fits the pressure's zeroGradient walls on the
    second ring; inverse distance is cellPoint's."""
    params = case_copy(CAVITY, tmp_path / "c", ["velocity_field=U", "pressure_field=p", "nodes=" + nodes])
    X, U, log = probe(params, "t0=10")
    assert len(X) == 200 and np.isfinite(U).all()
    if nodes == "least_squares":
        assert "Nodes by least squares, U: 441 rows" in log, log
        assert "80 fixed, 0 on the second ring (0 still deficient), 2 where zero wins" in log, log
        assert "Nodes by least squares, p: 441 rows" in log and "0 fixed, 80 on the second ring" in log, log
    else:
        assert "Nodes by inverse distance" in log, log


@needs_openfoam
def test_a_condition_that_changes_between_stamps_is_refused(tmp_path):
    """The node values are built from the first stamp's conditions: a lid
    that stops fixing the value at a later stamp is refused, naming both."""
    params = case_copy(CAVITY, tmp_path / "c")
    shutil.copytree(tmp_path / "c" / "10", tmp_path / "c" / "20")
    u = tmp_path / "c" / "20" / "U"
    text = u.read_text()
    u.write_text(text.replace("type            fixedValue;\n        value           uniform (1 0 0);",
                              "type            zeroGradient;", 1))
    assert "zeroGradient" in u.read_text()
    r = run_app(INTERPOL, params, PROBE, ["t0=20"], check=False, timeout=120)
    assert r.returncode == 2, r.stdout[-1000:] + r.stderr
    assert "U at 20: patch movingWall is zeroGradient" in r.stderr, r.stderr


@needs_openfoam
def test_tracers_run_on_the_channel(tmp_path):
    """A smoke test: two identical stamps (times 1 and 2), 250 RK4 steps."""
    params = case_copy(CHANNEL, tmp_path / "c")
    for t in ("1", "2"):
        shutil.copytree(tmp_path / "c" / "0", tmp_path / "c" / t)
    r = run_app(TRACERS, params, ["mode=openfoam", "Nrw=500", "Nrw_max=500", "int_order=1", "Dm=0",
                                  "dt=0.002", "t0=1", "T=1.5", "x0=0.5", "y0=0.5", "z0=0.5",
                                  "init_mode=points_xyz", "random=true", "seed=1", "dump_intv=0.5",
                                  "stat_intv=0.5", "checkpoint_intv=1e9"],
                env={"OMP_NUM_THREADS": "2"}, timeout=300)
    assert "Time = 1.5" in r.stdout, r.stdout[-2000:]


@needs_openfoam
def test_the_load_classes_the_boundary_facets_by_the_velocitys_condition(tmp_path):
    """The log counts the boundary facets by the class of their patch: the
    cavity's lid fixes a velocity along itself (a moving wall), its other
    walls none (no-slip); the channel's cyclic sides pair, and its walls,
    where the manufactured field crosses them, are neither kind of wall."""
    _, _, log = probe(case_copy(CAVITY, tmp_path / "a"), "t0=10")
    assert "Boundary facets: 60 no-slip wall (fixedWalls), 20 moving wall (movingWall), 0 cyclic, 0 other" in log, log
    _, _, log = probe(case_copy(CHANNEL, tmp_path / "b"), "t0=0")
    assert ("Boundary facets: 0 no-slip wall, 0 moving wall, 1152 cyclic (left right front back), "
            "576 other (bottom top)") in log, log


VECTOR = re.compile(r"\(([-+.\deE]+) ([-+.\deE]+) ([-+.\deE]+)\)")


def permute_vectors(path, perm):
    """Every (a b c) of an ascii OpenFOAM file with its components taken in the order perm."""
    text = path.read_text()
    path.write_text(VECTOR.sub(lambda m: "(%s)" % " ".join(m.group(1 + i) for i in perm), text))


def turned_cavity(d, normal):
    """The cavity with its empty pair normal to x or y: the points and the
    velocities rotated by a cyclic permutation of the axes (x) or mirrored
    by swapping y and z with the faces reversed to point out of their owners
    again (y); either way the in-plane axes are the original x and y."""
    params = case_copy(CAVITY, d)
    perm = {"x": (2, 0, 1), "y": (0, 2, 1)}[normal]
    for f in ("constant/polyMesh/points", "10/U", "0/U"):
        permute_vectors(d / f, perm)
    if normal == "y":
        faces = d / "constant" / "polyMesh" / "faces"
        faces.write_text(re.sub(r"(\d+)\(([\d ]+)\)",
                                lambda m: "%s(%s)" % (m.group(1), " ".join([m.group(2).split()[0]] +
                                                                          m.group(2).split()[:0:-1])),
                                faces.read_text()))
    return params


@needs_openfoam
@pytest.mark.parametrize("normal", ["x", "y"])
def test_a_2d_case_on_any_axis_is_the_same_field_and_traces(tmp_path, normal):
    """The empty pair's axis picks the plane: the turned cavity loads as
    triangles in the other two axes, interpol draws the same points there
    and finds the same velocity (to round-off: the geometry is computed in
    another component order), and tracers run on it as on the original."""
    X0, U0, _ = probe(case_copy(CAVITY, tmp_path / "z"), "t0=10")
    params = turned_cavity(tmp_path / normal, normal)
    X, U, log = probe(params, "t0=10")
    assert "Split triangles: 1600 for 400 cells" in log, log
    assert np.array_equal(X, X0)
    assert np.abs(U - U0).max() < 1e-14, np.abs(U - U0).max()
    assert np.abs(U0).max() > 0.1
    run = ["mode=openfoam", "Nrw=200", "Nrw_max=200", "int_order=1", "Dm=0", "dt=0.001", "t0=10", "T=10.2",
           "x0=0.05", "y0=0.05", "init_mode=points_xy", "random=false", "seed=1", "dump_intv=0.2",
           "stat_intv=0.2", "checkpoint_intv=1e9", "frozen_fields=true", "t_frozen=10"]
    ends = []
    for p in (tmp_path / "z" / "partrac_params.dat", params):
        run_app(TRACERS, p, run, env={"OMP_NUM_THREADS": "1"}, timeout=300)
        ends.append(all_dumps(p.parent)[10.2]["points"])
    assert np.abs(ends[0] - ends[1]).max() < 1e-12
    assert np.abs(ends[0] - all_dumps(tmp_path / "z")[10.0]["points"]).max() > 1e-3


@needs_openfoam
def test_a_3d_case_with_empty_patches_more_than_one_cell_thick_is_refused(tmp_path):
    """Empty sides on a 3D mesh: the case is not a 2D one and is not read as one."""
    params = case_copy(CHANNEL, tmp_path / "c")
    b = tmp_path / "c" / "constant" / "polyMesh" / "boundary.gz"
    text = gzip.decompress(b.read_bytes()).decode()
    for wall in ("bottom", "top"):
        i = text.index(wall)
        text = text[:i] + text[i:].replace("type            wall;", "type            empty;", 1)
    b.write_bytes(gzip.compress(text.encode()))
    refused(params, "is not one cell thick")


@needs_openfoam
def test_an_empty_pair_off_the_axes_is_refused(tmp_path):
    """The cavity turned by 30 degrees about x: its empty faces normal to no axis."""
    params = case_copy(CAVITY, tmp_path / "c")
    pts = tmp_path / "c" / "constant" / "polyMesh" / "points"
    c, s = np.cos(np.radians(30)), np.sin(np.radians(30))
    pts.write_text(VECTOR.sub(lambda m: "(%r %r %r)" % (float(m.group(1)),
                                                        c * float(m.group(2)) - s * float(m.group(3)),
                                                        s * float(m.group(2)) + c * float(m.group(3))),
                              pts.read_text()))
    r = run_app(INTERPOL, params, PROBE, ["t0=10"], check=False, timeout=120)
    assert r.returncode == 2 and "not normal to one axis" in r.stderr, r.stdout[-1000:] + r.stderr


def two_stamp_cavity(d):
    """The cavity with a second, earlier stamp at t = 5: every velocity halved."""
    params = case_copy(CAVITY, d)
    shutil.copytree(d / "10", d / "5")
    u = d / "5" / "U"
    u.write_text(VECTOR.sub(lambda m: "(%r %r %r)" % tuple(0.5 * float(m.group(i)) for i in (1, 2, 3)),
                            u.read_text()))
    return params


FROZEN = {"tracers": ["mode=openfoam", "init_mode=points_xy", "Nrw=200", "Nrw_max=200", "int_order=1", "Dm=0",
                      "dt=0.002", "x0=0.05", "y0=0.05", "random=false", "seed=1", "dump_intv=0.1",
                      "stat_intv=0.1", "checkpoint_intv=1e9", "frozen_fields=true", "num_threads=2"],
          "partrac": ["mode=openfoam", "init_mode=uniform_x", "Nrw=100", "Nrw_max=2000", "ds_max=0.01",
                      "ds_min=0.001", "La=0.08", "int_order=1", "Dm=0", "dt=0.002", "x0=0.05", "y0=0.07",
                      "random=false", "seed=1", "dump_intv=0.1", "stat_intv=0.1", "checkpoint_intv=1e9",
                      "frozen_fields=true", "num_threads=2"]}


@needs_openfoam
@pytest.mark.parametrize("name", ["tracers", "partrac"])
def test_frozen_fields_hold_the_stamp_at_t_frozen_whatever_the_time(tmp_path, name):
    """A steady case is traced as a series frozen at one stamp: run from t = 5
    on stamps at t = 5 and 10 with the fields frozen at either, the tracers
    see that stamp's field at every step, bit for bit what they see on a case
    holding that stamp alone, and not the time weight's blend of the two. A
    time outside the stamps is refused."""
    app_ = app(name)
    ends = {}
    for tf in (5, 10):
        two = two_stamp_cavity(tmp_path / ("two%d" % tf))
        run_app(app_, two, FROZEN[name], ["t0=5", "T=5.2", "t_frozen=%d" % tf], timeout=300)
        one = two_stamp_cavity(tmp_path / ("one%d" % tf))
        shutil.rmtree(tmp_path / ("one%d" % tf) / str(15 - tf))
        run_app(app_, one, FROZEN[name], ["t0=%d" % tf, "T=%g" % (tf + 0.2), "t_frozen=%d" % tf], timeout=300)
        a, b = all_dumps(tmp_path / ("two%d" % tf))[5.2], all_dumps(tmp_path / ("one%d" % tf))[tf + 0.2]
        assert len(a["points"]) >= 100
        for key in ("points", "u"):
            assert np.array_equal(a[key], b[key]), (tf, key)
        ends[tf] = a["points"]
    assert np.abs(ends[5] - ends[10]).max() > 1e-4
    # between the stamps: their blend, here three quarters of the t = 10 field
    run_app(app_, two_stamp_cavity(tmp_path / "mid"), FROZEN[name], ["t0=5", "T=5.2", "t_frozen=7.5"], timeout=300)
    one = two_stamp_cavity(tmp_path / "quarter")
    shutil.rmtree(tmp_path / "quarter" / "5")
    u = tmp_path / "quarter" / "10" / "U"
    u.write_text(VECTOR.sub(lambda m: "(%r %r %r)" % tuple(0.75 * float(m.group(i)) for i in (1, 2, 3)),
                            u.read_text()))
    run_app(app_, one, FROZEN[name], ["t0=10", "T=10.2", "t_frozen=10"], timeout=300)
    a, b = all_dumps(tmp_path / "mid")[5.2], all_dumps(tmp_path / "quarter")[10.2]
    assert np.abs(a["points"] - b["points"]).max() < 1e-12
    r = run_app(app_, two, FROZEN[name], ["t0=5", "T=5.2", "t_frozen=20"], check=False, timeout=300)
    assert r.returncode == 2 and "cannot be frozen at t = 20, outside the stamps' 5 to 10" in r.stderr, r.stderr


def refused(params, message, *args):
    r = run_app(INTERPOL, params, PROBE, ["t0=0"] + list(args), check=False, timeout=120)
    assert r.returncode == 2, r.stdout[-1000:] + r.stderr
    assert message in r.stderr, r.stderr


@needs_openfoam
def test_a_periodic_key_is_refused(tmp_path):
    refused(case_copy(CAVITY, tmp_path / "c", ["periodic_x=true"]), "cyclic patches say which axes")


@needs_openfoam
def test_a_case_without_system_is_refused(tmp_path):
    params = case_copy(CAVITY, tmp_path / "c")
    shutil.rmtree(tmp_path / "c" / "system")
    refused(params, "no system/ directory; not an OpenFOAM case")


@needs_openfoam
def test_a_decomposed_case_without_a_mesh_is_refused(tmp_path):
    params = case_copy(CAVITY, tmp_path / "c")
    shutil.move(str(tmp_path / "c" / "constant" / "polyMesh"), str(tmp_path / "c" / "processor0"))
    refused(params, "is decomposed and has no reconstructed mesh; run reconstructPar")


@needs_openfoam
def test_a_moving_mesh_is_refused(tmp_path):
    params = case_copy(CAVITY, tmp_path / "c")
    shutil.copytree(tmp_path / "c" / "constant" / "polyMesh", tmp_path / "c" / "10" / "polyMesh")
    refused(params, "10/polyMesh/points: a moving or changed mesh is not handled")


@needs_openfoam
@pytest.mark.parametrize("kind", ["wedge", "cyclicAMI", "processor"])
def test_a_patch_type_not_handled_is_refused(tmp_path, kind):
    params = case_copy(CAVITY, tmp_path / "c")
    b = tmp_path / "c" / "constant" / "polyMesh" / "boundary"
    text = b.read_text()
    b.write_text(text.replace("type            wall;\n        inGroups        List<word> 1(wall);",
                              "type            %s;" % kind, 1))
    assert kind in b.read_text()
    r = run_app(INTERPOL, params, PROBE, ["t0=10"], check=False, timeout=120)
    assert r.returncode == 2 and kind in r.stderr, r.stdout[-1000:] + r.stderr


@needs_openfoam
def test_binary_with_64_bit_labels_is_refused(tmp_path):
    """OpenFOAM.org reads such a file as 32-bit; the message names the conversion."""
    params = case_copy(CHANNEL, tmp_path / "c")
    owner = tmp_path / "c" / "constant" / "polyMesh" / "owner.gz"
    data = gzip.decompress(owner.read_bytes()).replace(b"format      binary;",
                                                        b"format      binary;\n    arch        \"LSB;label=64;scalar=64\";", 1)
    owner.write_bytes(gzip.compress(data))
    refused(params, "foamFormatConvert")


@needs_openfoam
def test_a_baffle_is_refused(tmp_path):
    """Two boundary faces on one set of points: the facet table would pair them."""
    params = case_copy(CAVITY, tmp_path / "c")
    poly = tmp_path / "c" / "constant" / "polyMesh"
    faces = poly / "faces"
    text = faces.read_text()
    # the last front face made a copy of the first boundary face
    lines = text.splitlines()
    i0 = next(i for i, l in enumerate(lines) if l.startswith("4(")) + 760
    last = max(i for i, l in enumerate(lines) if l.startswith("4("))
    lines[last] = lines[i0][:2] + " ".join(reversed(lines[i0][2:-1].split())) + ")"
    faces.write_text("\n".join(lines) + "\n")
    r = run_app(INTERPOL, params, PROBE, ["t0=10"], check=False, timeout=120)
    assert r.returncode == 2 and "baffle" in r.stderr, r.stdout[-1000:] + r.stderr


@needs_openfoam
def test_a_missing_pressure_names_the_way_out(tmp_path):
    """A case without p (a p_rgh solver's) is refused with the two settings
    that read it anyway."""
    params = case_copy(CAVITY, tmp_path / "c")
    os.remove(tmp_path / "c" / "10" / "p")
    r = run_app(INTERPOL, params, PROBE, ["t0=10"], check=False, timeout=120)
    assert r.returncode == 2, r.stdout[-1000:] + r.stderr
    assert "10/p: no such file; set pressure_field, or ignore_pressure=true" in r.stderr, r.stderr


ALPHA = (1.0, 3.0, 5.0)   # alpha.water = 1 + 3 x + 5 y on the cavity's [0, 0.1]^2, water everywhere
WALLS = ("    movingWall\n    {\n        type            zeroGradient;\n    }\n"
         "    fixedWalls\n    {\n        type            zeroGradient;\n    }\n")


def write_alpha(case, time, values=None, walls=WALLS):
    """alpha.water in the cell centres of time, the linear ALPHA or values(C)
    of the centres C, the walls' conditions as given; the mesh."""
    import openfoam_expected as oe
    m = oe.read_mesh(str(case))
    oe.geometry(m)
    C = np.array(m["C"])
    a = values(C) if values else ALPHA[0] + ALPHA[1] * C[:, 0] + ALPHA[2] * C[:, 1]
    (case / time / "alpha.water").write_text(
        "FoamFile\n{\n    format      ascii;\n    class       volScalarField;\n    location    \"%s\";\n"
        "    object      alpha.water;\n}\n\ndimensions      [0 0 0 0 0 0 0];\n\n"
        "internalField   nonuniform List<scalar> %d\n(\n%s\n);\n\nboundaryField\n{\n%s"
        "    frontAndBack\n    {\n        type            empty;\n    }\n}\n"
        % (time, len(a), "\n".join(repr(float(v)) for v in a), walls))
    m["alpha"] = a
    return m


def probe_phi(params, *args, nrw=200):
    """interpol's points, its phase field and the phase field's gradient where
    written, and its stdout."""
    import h5py
    r = run_app(INTERPOL, params, PROBE, ["Nrw=%d" % nrw] + list(args), env={"OMP_NUM_THREADS": "1"}, timeout=300)
    [out] = list(params.parent.rglob("interpolation.h5part"))
    with h5py.File(out, "r") as h:
        g = h[list(h.keys())[0]]
        G = np.c_[g["dphi_dx"], g["dphi_dy"], g["dphi_dz"]] if "dphi_dx" in g else None
        return np.c_[g["x"], g["y"], g["z"]], np.array(g["phi"]), G, r.stdout


def openfoam_bashrc():
    """OpenFOAM's etc/bashrc, the build's or the environment's; None without."""
    dirs = [os.environ.get("WM_PROJECT_DIR", "")]
    cache = os.path.join(os.path.dirname(BIN), "CMakeCache.txt")
    if os.path.exists(cache):
        m = re.search(r"^OPENFOAM_DIR:PATH=(.*)$", open(cache).read(), re.M)
        if m:
            dirs.insert(0, m.group(1))
    return next((os.path.join(d, "etc", "bashrc") for d in dirs
                 if d and os.path.exists(os.path.join(d, "etc", "bashrc"))), None)


def cellpoint(case, field, time, X):
    """OpenFOAM's own cellPoint of a scalar field at the points X, by its probes."""
    pts = "\n".join("        (%.17g %.17g %.17g)" % tuple(p) for p in X)
    (case / "system" / "probesCP").write_text(
        "FoamFile\n{\n    format      ascii;\n    class       dictionary;\n    object      probesCP;\n}\n"
        "type            probes;\nlibs            (\"libsampling.so\");\nwriteControl    timeStep;\n"
        "writeInterval   1;\nfields          (%s);\nfixedLocations  true;\ninterpolationScheme cellPoint;\n"
        "probeLocations\n    (\n%s\n    );\n" % (field, pts))
    r = subprocess.run(["bash", "-c", "source %s > /dev/null 2>&1; foamPostProcess -func probesCP -time %s"
                        % (openfoam_bashrc(), time)], cwd=case, capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    [out] = list((case / "postProcessing" / "probesCP").rglob(field))
    rows = [l for l in open(out).read().splitlines() if l.strip() and not l.startswith("#")]
    return np.array(rows[-1].split()[1:], dtype=float)


@needs_openfoam
@pytest.mark.parametrize("nodes", ["least_squares", "inverse_distance"])
def test_the_phase_fields_nodes_are_cellpoints_whatever_the_node_rule(tmp_path, nodes):
    """nodes= is for U and p: the phase field's node values are inverse
    distance's, cellPoint's (OpenFOAM's probes, where OpenFOAM runs), and with
    water in every cell no centre needs moving."""
    params = case_copy(CAVITY, tmp_path / "c", ["velocity_field=U", "pressure_field=p", "phase_field=alpha.water",
                                                "nodes=" + nodes])
    m = write_alpha(tmp_path / "c", "10")
    X, phi, G, log = probe_phi(params, "t0=10")
    assert len(X) == 200
    assert ("Nodes by inverse distance, alpha.water" in log) == (nodes == "least_squares"), log
    assert "Gradient of alpha.water by least squares: 441 rows" in log, log
    assert phase_line(log)[4:] == [400, 0, 0, 0, 0], log   # every cell pure, none moved
    assert np.abs(phi - (ALPHA[0] + ALPHA[1] * X[:, 0] + ALPHA[2] * X[:, 1])).max() > 1e-4
    assert np.abs(G[:, :2] - ALPHA[1:]).max() < 1e-10 and np.all(G[:, 2] == 0)
    if openfoam_bashrc() is None:
        return
    z = [x[2] for x in m["X"]]
    P = X.copy()
    P[:, 2] = 0.5 * (min(z) + max(z))
    assert np.abs(phi - cellpoint(tmp_path / "c", "alpha.water", "10", P)).max() < 1e-12


# A water layer below a tilted line, with a round bubble in it, as a tanh step
# 1.2 cells wide across the signed distance to the nearer interface
LINE_N, LINE_C = np.array([0.3, 0.954]) / np.hypot(0.3, 0.954), 0.07
BUBBLE_X, BUBBLE_R, EPS = np.array([0.05, 0.03]), 0.018, 0.006


def layer_distance(X):
    """The signed distance to the interface, water positive, and its gradient."""
    d_line = LINE_C - X[:, :2] @ LINE_N
    r = X[:, :2] - BUBBLE_X
    rn = np.hypot(r[:, 0], r[:, 1])
    d_bub = rn - BUBBLE_R
    line = d_line < d_bub
    grad = np.where(line[:, None], -LINE_N[None, :], r / rn[:, None])
    return np.minimum(d_line, d_bub), grad, np.abs(d_line - d_bub)


def layer(X):
    return 0.5 * (1 + np.tanh(layer_distance(X)[0] / EPS))


def phase_line(log):
    """The numbers of the log's phase line: the volume, the target, the
    relative errors at the nodes and corrected, and the cells pure,
    corrected, clipped, unreachable and without a centre, 0 where the line
    leaves a count out."""
    m = re.search(r"Phase alpha\.water at 10: volume above 1/2 (\S+) for (\S+), relative error (\S+) at the nodes, "
                  r"(\S+) corrected; .*; cells: (.*)", log)
    assert m, log
    counts = [re.search(r"(?:^| )(\d+) " + re.escape(what) + r"(?:,|$)", m.group(5))
              for what in ["pure", "corrected", "clipped to [0, 1]", "unreachable", "off without a centre"]]
    return [float(v) for v in m.groups()[:4]] + [int(c.group(1)) if c else 0 for c in counts]


@needs_openfoam
@pytest.mark.parametrize("split", ["12", "6"])
def test_the_volume_above_one_half_is_the_cells_volume_fractions(tmp_path, split):
    """W12 moves the centre of each mixed cell, within [0, 1], until its
    volume above 1/2 is V_c alpha_c; a cell within 1e-6 of 0 or 1 is pure and
    keeps its centre. Here every mixed cell reaches its target (none clipped
    or unreachable), so what is left of the sum's error is the pure cells':
    under 1e-6 of the sum V_c alpha_c, against 1e-2 at the nodes; the pure
    count is the mesh's, the target the mesh's sum. W6, with no centre, logs
    its error and moves nothing. The phase gradient, P1 on the least-squares
    fit's, points along the exact normal near the interface, within 5
    degrees."""
    params = case_copy(CAVITY, tmp_path / "c", ["velocity_field=U", "pressure_field=p", "phase_field=alpha.water",
                                                "split=" + split])
    m = write_alpha(tmp_path / "c", "10", layer)
    X, phi, G, log = probe_phi(params, "t0=10", nrw=4000)
    z = np.array(m["X"])[:, 2]
    target = np.dot(m["V"], m["alpha"]) / (z.max() - z.min())
    above, logged, before, after, pure, corrected, clipped, unreachable, unfanned = phase_line(log)
    # the log's six digits
    assert abs(logged - target) < 1e-5 * target
    assert abs(before) > 1e-2
    a = np.array(m["alpha"])
    if split == "6":
        assert unfanned == len(a) and pure + corrected + clipped + unreachable == 0 and after == before
        return
    assert pure == np.sum((a < 1e-6) | (a > 1 - 1e-6)) > 0
    assert (corrected, clipped, unreachable, unfanned) == (len(a) - pure, 0, 0, 0)
    assert abs(after) < 1e-6
    assert np.all((phi >= 0) & (phi <= 1))
    d, normal, kink = layer_distance(X)
    near = (np.abs(d) < EPS) & (kink > 3 * EPS)
    assert near.sum() > 100
    unit = G[:, :2] / np.linalg.norm(G[:, :2], axis=1)[:, None]
    angle = np.degrees(np.arccos(np.clip(np.sum(unit * normal, axis=1), -1, 1)))
    assert angle[near].max() < 5, np.sort(angle[near])[-5:]
    assert np.all(G[:, 2] == 0)


@needs_openfoam
def test_a_scalar_condition_that_cannot_be_built_is_not_data(tmp_path):
    """contactAngle's library is not loaded in the cavity: alpha.water's lid
    is taken as not fixing the value, said once, and the field read as with a
    zeroGradient lid, the condition's written value where it gives one."""
    ref = case_copy(CAVITY, tmp_path / "a", ["velocity_field=U", "pressure_field=p", "phase_field=alpha.water"])
    write_alpha(tmp_path / "a", "10", layer)
    X0, phi0, G0, _ = probe_phi(ref, "t0=10")
    params = case_copy(CAVITY, tmp_path / "b", ["velocity_field=U", "pressure_field=p", "phase_field=alpha.water"])
    lid = WALLS.replace("zeroGradient;\n    }\n    fixedWalls", "contactAngle;\n        theta0          90;\n"
                        "        limit           gradient;\n    }\n    fixedWalls", 1)
    write_alpha(tmp_path / "b", "10", layer, lid)
    X, phi, G, log = probe_phi(params, "t0=10")
    assert log.count("alpha.water: patch movingWall is contactAngle, which OpenFOAM cannot build here "
                     "(its library is not loaded); its faces are not data") == 1, log
    assert np.array_equal(X, X0) and np.array_equal(phi, phi0) and np.array_equal(G, G0)


@needs_openfoam
def test_a_velocity_condition_that_cannot_be_built_is_refused_naming_its_library(tmp_path):
    params = case_copy(CAVITY, tmp_path / "c")
    u = tmp_path / "c" / "10" / "U"
    u.write_text(u.read_text().replace("type            fixedValue;", "type            waveVelocity;", 1))
    refused(params, "10/U: patch movingWall is waveVelocity, which OpenFOAM cannot build here: its library is not "
                    "loaded, probably libwaves.so; add libs (\"libwaves.so\"); to system/controlDict", "t0=10")


@needs_openfoam
def test_without_phase_field_no_phase_field_is_read(tmp_path):
    """No key: no phi written, alpha.water in the case or not."""
    import h5py
    params = case_copy(CAVITY, tmp_path / "c")
    write_alpha(tmp_path / "c", "10")
    _, _, log = probe(params, "t0=10")
    assert "alpha.water" not in log
    [out] = list(params.parent.rglob("interpolation.h5part"))
    with h5py.File(out, "r") as h:
        assert "phi" not in h[list(h.keys())[0]]


@needs_openfoam
def test_a_phase_field_that_is_not_a_cell_scalar_is_refused(tmp_path):
    """OpenFOAM's phi, the face flux, and a vector."""
    params = case_copy(CAVITY, tmp_path / "c", ["phase_field=phi"])
    shutil.copy(tmp_path / "c" / "10" / "p", tmp_path / "c" / "10" / "phi")
    f = tmp_path / "c" / "10" / "phi"
    f.write_text(f.read_text().replace("class       volScalarField;", "class       surfaceScalarField;"))
    refused(params, "10/phi is a surfaceScalarField, values on the faces, as OpenFOAM's phi, the volume flux", "t0=10")
    params = case_copy(CAVITY, tmp_path / "v", ["phase_field=U"])
    refused(params, "10/U has 3 components, not a scalar; phase_field names a volScalarField such as alpha.water "
                    "(OpenFOAM's phi is the volume flux through the faces)", "t0=10")


@needs_openfoam
def test_a_time_without_the_phase_field_is_refused_at_load(tmp_path):
    """Refused at load, naming the time."""
    params = case_copy(CAVITY, tmp_path / "c", ["phase_field=alpha.water"])
    write_alpha(tmp_path / "c", "10")
    for f in ("U", "p"):
        os.makedirs(tmp_path / "c" / "20", exist_ok=True)
        shutil.copy(tmp_path / "c" / "10" / f, tmp_path / "c" / "20" / f)
    refused(params, "/20/alpha.water: no such file; the phase field must be in every time directory", "t0=10")


@needs_openfoam
def test_include_phi_is_refused_for_phase_field(tmp_path):
    refused(case_copy(CAVITY, tmp_path / "c", ["include_phi=true"]),
            "mode openfoam reads the phase field phase_field names, as phase_field=alpha.water")


@needs_openfoam
def test_a_reader_of_another_build_is_refused(tmp_path):
    """The apps and the reader pass C++ structs: a reader whose layout stamp
    is not the app's is refused before any call, as a plugin compiled against
    another openfoam_load.hpp would crash. Here a stand-in library exporting
    a wrong stamp, loaded through PARTRAC_OPENFOAM_PLUGIN."""
    cc = shutil.which("cc") or shutil.which("gcc")
    if not cc:
        pytest.skip("no C compiler for the stand-in reader")
    src = tmp_path / "stamp.c"
    src.write_text("unsigned long long partrac_openfoam_layout(void) { return 7; }\n")
    lib = tmp_path / "libstamp.so"
    subprocess.run([cc, "-shared", "-fPIC", "-o", str(lib), str(src)], check=True)
    params = case_copy(CAVITY, tmp_path / "c")
    r = run_app(INTERPOL, params, PROBE, ["t0=10"], env={"PARTRAC_OPENFOAM_PLUGIN": str(lib)},
                check=False, timeout=120)
    assert r.returncode == 2, r.stdout[-1000:] + r.stderr
    assert "libstamp.so: an OpenFOAM reader of another build (layout 7, this build's " in r.stderr, r.stderr


@needs_openfoam
@pytest.mark.slow
def test_an_installed_tree_reads_the_cavity_with_its_own_reader(tmp_path):
    """cmake --install into a fresh prefix: the installed interpol, run with
    an empty environment (no OpenFOAM sourced, no library path), loads the
    installed reader, not the build tree's beside it, finds OpenFOAM through
    the reader's RPATH, and probes the cavity exactly as the build does."""
    build = os.path.dirname(BIN)
    if not shutil.which("cmake") or not os.path.exists(os.path.join(build, "CMakeCache.txt")):
        pytest.skip("no cmake, or %s is not a build tree" % build)
    prefix = tmp_path / "inst"
    r = subprocess.run(["cmake", "--install", build, "--prefix", str(prefix)], capture_output=True, text=True,
                       timeout=300)
    assert r.returncode == 0, r.stdout + r.stderr
    X0, U0, _ = probe(case_copy(CAVITY, tmp_path / "b"), "t0=10")
    params = case_copy(CAVITY, tmp_path / "c")
    r = subprocess.run([str(prefix / "bin" / "interpol"), str(params)] + PROBE + ["t0=10"], capture_output=True,
                       text=True, timeout=300, cwd=params.parent,
                       env={"HOME": os.environ.get("HOME", str(tmp_path)), "LD_DEBUG": "files"})
    assert r.returncode == 0, r.stdout[-1000:] + r.stderr[-3000:]
    [installed] = list(prefix.rglob("libpartrac_openfoam.so"))
    loaded = set(re.findall(r"file=(\S*libpartrac_openfoam\S*)", r.stderr))
    assert loaded == {str(installed)}, loaded
    import h5py
    [out] = list(params.parent.rglob("interpolation.h5part"))
    with h5py.File(out, "r") as h:
        g = h[list(h.keys())[0]]
        assert np.array_equal(np.c_[g["x"], g["y"], g["z"]], X0)
        assert np.array_equal(np.c_[g["ux"], g["uy"], g["uz"]], U0)


@pytest.mark.skipif(built_with_openfoam() or not os.path.exists(INTERPOL),
                    reason="the reader is built" if built_with_openfoam() else "interpol is not built")
def test_the_mode_needs_the_reader(tmp_path):
    params = case_copy(CAVITY, tmp_path / "c")
    refused(params, "mode openfoam needs a build with PARTRAC_ENABLE_OPENFOAM=ON.")


@needs_openfoam
def test_a_cyclic_pair_that_is_not_a_translation_is_refused(tmp_path):
    """One point of a cyclic side moved along its plane: its image misses it."""
    params = case_copy(CHANNEL, tmp_path / "c")
    path = tmp_path / "c" / "constant" / "polyMesh" / "points.gz"
    data = bytearray(gzip.decompress(path.read_bytes()))
    head = data.index(b"2197\n(") + len(b"2197\n(")
    X = np.frombuffer(bytes(data[head:head + 2197 * 24]), dtype="<f8").reshape(-1, 3).copy()
    i = np.flatnonzero((X[:, 0] == 1.0) & (X[:, 1] > 0.1) & (X[:, 1] < 0.9) & (X[:, 2] > 0.1) & (X[:, 2] < 0.9))[0]
    X[i, 1] += 1e-3
    data[head:head + 2197 * 24] = X.astype("<f8").tobytes()
    path.write_bytes(gzip.compress(bytes(data)))
    refused(params, "is not a translation")


@needs_openfoam
def test_a_condition_without_a_null_constructor_is_read(tmp_path):
    """Whether a condition fixes the value is asked of the condition itself;
    one that must be built from its dictionary (movingWallVelocity) is, and
    the lid stays a moving wall with the same field."""
    X0, U0, _ = probe(case_copy(CAVITY, tmp_path / "a"), "t0=10")
    params = case_copy(CAVITY, tmp_path / "b")
    u = tmp_path / "b" / "10" / "U"
    u.write_text(u.read_text().replace("type            fixedValue;", "type            movingWallVelocity;", 1))
    X, U, log = probe(params, "t0=10")
    assert "20 moving wall (movingWall)" in log, log
    assert np.array_equal(X, X0) and np.array_equal(U, U0)


FOAM_HEADER = "FoamFile\n{\n    format      ascii;\n    class       %s;\n    object      %s;\n}\n\n"


def strip_case(d, split, nx=5):
    """A 2D ascii case of nx x 1 unit cells, 0.1 thick with empty front and
    back: no-slip walls left, at the bottom and over the first cell's top
    (the lid), open (zeroGradient) at the rest of the top and the right, a
    uniform (1 0 0) inside. The first cell has every point on a wall, the
    others a free point at the top. The cavity's system and constant are
    reused."""
    params = case_copy(CAVITY, d, ["velocity_field=U", "ignore_pressure=true", "split=" + split])
    for t in ("0", "10"):
        shutil.rmtree(d / t)
    poly = d / "constant" / "polyMesh"
    shutil.rmtree(poly)
    poly.mkdir()
    pid = lambda i, j, k: i + (nx + 1) * (j + 2 * k)
    X = [(i, j, 0.1 * k) for k in range(2) for j in range(2) for i in range(nx + 1)]
    xface = lambda i: [pid(i, 0, 0), pid(i, 1, 0), pid(i, 1, 1), pid(i, 0, 1)]   # normal +x
    yface = lambda i, j: [pid(i, j, 0), pid(i, j, 1), pid(i + 1, j, 1), pid(i + 1, j, 0)]   # normal +y
    zface = lambda i, k: [pid(i, 0, k), pid(i + 1, 0, k), pid(i + 1, 1, k), pid(i, 1, k)]   # normal +z
    rev = lambda f: [f[0]] + f[:0:-1]
    faces = [xface(i + 1) for i in range(nx - 1)]
    owner = list(range(nx - 1))
    neighbour = list(range(1, nx))
    patches = [("left", "wall", [rev(xface(0))], [0]),
               ("bottom", "wall", [rev(yface(i, 0)) for i in range(nx)], list(range(nx))),
               ("lid", "wall", [yface(0, 1)], [0]),
               ("top", "patch", [yface(i, 1) for i in range(1, nx)], list(range(1, nx))),
               ("right", "patch", [xface(nx)], [nx - 1]),
               ("frontAndBack", "empty", [rev(zface(i, 0)) for i in range(nx)] + [zface(i, 1) for i in range(nx)],
                list(range(nx)) * 2)]
    bounds = []
    for name, kind, fs, cells in patches:
        bounds.append("%s\n{\n    type %s;\n    nFaces %d;\n    startFace %d;\n}\n" % (name, kind, len(fs), len(faces)))
        faces += fs
        owner += cells
    lists = {"points": ("vectorField", ["(%r %r %r)" % x for x in X]),
             "faces": ("faceList", ["4(%s)" % " ".join(map(str, f)) for f in faces]),
             "owner": ("labelList", list(map(str, owner))),
             "neighbour": ("labelList", list(map(str, neighbour))),
             "boundary": ("polyBoundaryMesh", bounds)}
    for name, (cls, rows) in lists.items():
        (poly / name).write_text(FOAM_HEADER % (cls, name) + "%d\n(\n%s)\n" % (len(rows), "\n".join(rows)))
    conds = {"left": "noSlip", "bottom": "noSlip", "lid": "noSlip", "top": "zeroGradient", "right": "zeroGradient",
             "frontAndBack": "empty"}
    (d / "0").mkdir()
    (d / "0" / "U").write_text(FOAM_HEADER % ("volVectorField", "U") +
                               "dimensions [0 1 -1 0 0 0 0];\ninternalField uniform (1 0 0);\nboundaryField\n{\n" +
                               "".join("%s { type %s; }\n" % kv for kv in conds.items()) + "}\n")
    return params


@needs_openfoam
@pytest.mark.parametrize("split", ["12", "6"])
def test_cells_with_every_point_on_a_wall_are_warned_of(tmp_path, split):
    """A cell whose every point is on a no-slip wall (a gap one cell wide)
    has a zero velocity throughout under split=6 and only its centre free
    under split=12: the load says how many there are, and goes on."""
    X, U, log = probe(strip_case(tmp_path / "c", split), "t0=0")
    assert len(X) == 200 and np.isfinite(U).all()
    assert "Warning: 1 of 5 cells (20%) have every point on a no-slip wall" in log, log
    assert ("their velocity is zero throughout" if split == "6" else "only their centre values are free") in log, log


# --- solved cases -----------------------------------------------------------
#
# The cases below are OpenFOAM solutions too large to check in: an obstacle
# array in 2D (rsa2d, 61,801 cells, cyclic, the obstacles no-slip, a Stokes
# flow under an oblique body force), four spheres in a channel (spheres3d,
# snappyHexMesh), the jittered and graded channels (Poiseuille), a Couette
# channel with a moving wall. PARTRAC_OPENFOAM_CASES lists the folders that
# hold them (os.pathsep between); each test skips when its case is absent.

CASE_DIRS = [d for d in os.environ.get("PARTRAC_OPENFOAM_CASES", "").split(os.pathsep) if d]


def solved(name):
    """The folder of the solved case `name`, or a skip."""
    for d in CASE_DIRS:
        if os.path.isdir(os.path.join(d, name, "constant")):
            return os.path.join(d, name)
    pytest.skip("the solved case %s is not in PARTRAC_OPENFOAM_CASES" % name)


def linked(src, time, d, lines, extra=()):
    """d with the case's constant, system and time linked in and a parameter file of `lines`."""
    d.mkdir(parents=True)
    for sub in ("constant", "system", time) + tuple(extra):
        os.symlink(os.path.realpath(os.path.join(src, sub)), d / sub)
    (d / "partrac_params.dat").write_text("".join(l + "\n" for l in lines))
    return d / "partrac_params.dat"


def probe_rows(params, time, nrw=2000):
    """interpol's own points (seed 1, one thread) in the case at `time`: X, U, J, log."""
    return probe(params, "t0=" + time, "Nrw=%d" % nrw, gradient=True)


@needs_openfoam
@pytest.mark.parametrize("name,time,split,cyclic", [("channel3d_jitter/n32", "1329", "12", 8192),
                                                    ("channel3d_jitter/n32", "1329", "6", 8192),
                                                    ("rsa2d", "189", "12", 844), ("spheres3d", "230", "12", 4608)])
def test_every_cyclic_facet_is_paired_on_the_solved_cases(tmp_path, name, time, split, cyclic):
    """A cyclic facet without its image is refused at load, so a load that
    completes has them all paired; the log counts them, two triangles a
    cyclic face in 3D (the channel's four 32 x 32 sides) and an edge in 2D."""
    params = linked(solved(name), time, tmp_path / "c", ["velocity_field=U", "ignore_pressure=true", "split=" + split])
    X, U, J, log = probe_rows(params, time, 200)
    assert ", %d cyclic (" % cyclic in log, log


def closed_form_errors(X, U, J, u, du):
    """RMS of |u - u_exact| and of |grad u - grad u_exact|_F for a flow u(z) along x."""
    z = X[:, 2]
    ue = np.zeros_like(U)
    ue[:, 0] = u(z)
    Je = np.zeros_like(J)
    Je[:, 0, 2] = du(z)
    e = np.linalg.norm(U - ue, axis=1)
    eg = np.linalg.norm((J - Je).reshape(len(z), 9), axis=1)
    return np.sqrt(np.mean(e ** 2)), np.sqrt(np.mean(eg ** 2))


@needs_openfoam
def test_poiseuille_on_the_graded_channel(tmp_path):
    """The graded channel (grading 4 toward both walls, n = 32) solved under a
    body force, against u = 4 z (1 - z): the error is the solver's and the
    split's, pinned at the values the construction gives, 10% either way;
    cellPoint's own (the inverse-distance nodes) is larger in value and
    smaller in gradient on this mesh."""
    src = solved("channel3d_graded/n32")
    errors = {}
    for nodes in ("least_squares", "inverse_distance"):
        params = linked(src, "3000", tmp_path / nodes, ["velocity_field=U", "ignore_pressure=true", "nodes=" + nodes])
        X, U, J, log = probe_rows(params, "3000")
        assert len(X) == 2000
        errors[nodes] = closed_form_errors(X, U, J, lambda z: 4 * z * (1 - z), lambda z: 4 - 8 * z)
    value, gradient = errors["least_squares"]
    assert 0.9 * 7.2e-4 < value < 1.1 * 7.2e-4, value
    assert 0.9 * 0.204 < gradient < 1.1 * 0.204, gradient
    assert errors["inverse_distance"][0] > value and errors["inverse_distance"][1] < gradient


@needs_openfoam
def test_couette_with_a_moving_wall_is_exact(tmp_path):
    """The same graded channel with no force and its top moving
    (movingWallVelocity): the solution u = z is linear, which the finite
    volumes reproduce to the solver's tolerance and the least-squares nodes
    exactly; the top is classed a moving wall, the bottom a no-slip one."""
    params = linked(solved("couette"), "3000", tmp_path / "c", ["velocity_field=U", "ignore_pressure=true"])
    X, U, J, log = probe_rows(params, "3000")
    assert "2048 no-slip wall (bottom), 2048 moving wall (top), 8192 cyclic" in log, log
    value, gradient = closed_form_errors(X, U, J, lambda z: z, np.ones_like)
    assert value < 1e-12 and gradient < 1e-11, (value, gradient)


def probes_output(case, func, field):
    """OpenFOAM's probes at the case's probe_points.txt: the last time's values
    of `field` (a row a point), NaN where a point is outside the mesh."""
    [f] = [os.path.join(dp, field) for dp, _, fs in os.walk(os.path.join(case, "postProcessing", func)) if field in fs]
    row = [l for l in open(f).read().splitlines() if l.strip() and not l.startswith("#")][-1]
    v = np.array(row.split(None, 1)[1].replace("(", " ").replace(")", " ").split(), dtype=float)
    v = v.reshape(len(np.loadtxt(os.path.join(case, "probe_points.txt"))), -1)
    v[np.abs(v) > 1e200] = np.nan
    return v


def at_probe_points(case, X, dim):
    """The row of interpol's output at each of the case's probe points."""
    from scipy.spatial import cKDTree
    d, i = cKDTree(X[:, :dim]).query(np.loadtxt(os.path.join(case, "probe_points.txt"))[:, :dim])
    assert d.max() < 1e-12
    return i


@needs_openfoam
def test_the_pressure_is_cellpoints(tmp_path):
    """The pressure through the inverse-distance nodes is OpenFOAM's cellPoint
    p at every probed point of the cavity, its zeroGradient walls included."""
    import h5py
    src = solved("cavity_p")
    params = linked(src, "10", tmp_path / "c", ["velocity_field=U", "pressure_field=p", "wall_p2=none",
                                                 "nodes=inverse_distance"])
    X, U, log = probe(params, "t0=10", "Nrw=2000")
    [out] = list(params.parent.rglob("interpolation.h5part"))
    with h5py.File(out, "r") as h:
        p = np.array(h[list(h.keys())[0]]["p"])
    ref = probes_output(src, "probesCPp", "p")[:, 0]
    i = at_probe_points(src, X, 2)
    assert np.isfinite(ref).all() and np.abs(ref).max() > 1e-2
    assert np.abs(p[i] - ref).max() < 1e-13, np.abs(p[i] - ref).max()


@needs_openfoam
def test_a_symmetry_plane_is_cellpoints_and_mirrors_exactly(tmp_path):
    """The lower half of the graded channel under the same force, its top a
    symmetryPlane: the inverse-distance nodes, whose points on the plane lose
    their normal component, are OpenFOAM's cellPoint; the least-squares nodes,
    which fit over the cells and their mirror images, are exact for a linear
    field with the plane's symmetry, u = (1, -0.5, 2 (z - 0.5))."""
    src = solved("halfchannel")
    params = linked(src, "3000", tmp_path / "a", ["velocity_field=U", "ignore_pressure=true", "wall_p2=none",
                                                   "nodes=inverse_distance"])
    X, U, log = probe(params, "t0=3000", "Nrw=2000")
    assert "2048 other (top)" in log, log
    ref = probes_output(src, "probesCP", "U")
    i = at_probe_points(src, X, 3)
    assert np.isfinite(ref).all()
    assert np.abs(U[i] - ref).max() < 1e-13, np.abs(U[i] - ref).max()
    params = linked(os.path.join(src, "mf_symlinear"), "0", tmp_path / "b",
                    ["velocity_field=U", "ignore_pressure=true", "wall_p2=none"])
    X, U, J, log = probe(params, "t0=0", "Nrw=2000", gradient=True)
    assert "1024 mirrored" in log, log
    ue = np.c_[np.ones(len(X)), -0.5 * np.ones(len(X)), 2 * (X[:, 2] - 0.5)]
    Je = np.zeros((3, 3))
    Je[2, 2] = 2.0
    assert np.abs(U - ue).max() < 1e-13 and np.abs(J - Je).max() < 1e-11


@needs_openfoam
def test_frozen_fields_on_the_steady_graded_channel(tmp_path):
    """A steady case's last iteration traced with frozen_fields, through
    tracers and partrac, beside an earlier iteration (here the field halved):
    the run sees the t_frozen field throughout, as on the case holding it alone."""
    src = solved("channel3d_graded/n32")
    ends = {}
    for name, init in (("tracers", ["init_mode=points_xyz", "Nrw=500", "Nrw_max=500"]),
                       ("partrac", ["init_mode=uniform_x", "Nrw=100", "Nrw_max=2000", "ds_max=0.05",
                                    "ds_min=0.001", "La=0.8"])):
        for tag, extra in (("two", ("1500",)), ("one", ())):
            d = tmp_path / (name + tag)
            params = linked(src, "3000", d, ["velocity_field=U", "ignore_pressure=true"])
            if extra:
                (d / "1500").mkdir()
                u = open(os.path.join(src, "3000", "U")).read()
                (d / "1500" / "U").write_text(VECTOR.sub(
                    lambda m: "(%r %r %r)" % tuple(0.5 * float(m.group(i)) for i in (1, 2, 3)), u))
            run_app(app(name), params, init, ["mode=openfoam", "int_order=1", "Dm=0", "dt=0.01",
                                              "t0=%s" % (extra[0] if extra else "3000"),
                                              "T=%g" % ((1500 if extra else 3000) + 0.5), "x0=0.5", "y0=0.5",
                                              "z0=0.3", "random=false", "seed=1", "dump_intv=0.5",
                                              "stat_intv=0.5", "checkpoint_intv=1e9", "frozen_fields=true",
                                              "t_frozen=3000", "num_threads=2"], timeout=600)
            D = all_dumps(d)
            ends[tag] = D[max(D)]
        for key in ("points", "u"):
            assert np.array_equal(ends["two"][key], ends["one"][key]), (name, key)
        assert np.abs(ends["one"]["u"][:, 0]).max() > 0.5


# The near-wall runs: tracers on a lattice through a checkpoint, Dm = 0, RK4,
# dt one mean cell at the fastest point, five throughflows, at rest where
# |u| < 1e-3 of the mean speed along the force

NEAR_WALL = {"rsa2d": dict(time="189", dt=0.6698121145098501, steps=7898, rest=9.452259612219838e-06, lattice=70),
             "spheres3d": dict(time="230", dt=0.6301918904131909, steps=368, rest=2.158362824587462e-05, lattice=16)}
OFFSETS = [(0.5, 0.5, 0.5), (0.25, 0.75, 0.4), (0.75, 0.3, 0.65)]


def foam_list(path):
    """The list of an ascii OpenFOAM mesh or field file, as text rows."""
    text = open(path).read()
    head = text.index("}", text.index("FoamFile")) + 1
    body = text[head:]
    m = re.search(r"(\d+)\s*\(", body)
    n = int(m.group(1))
    rows = re.findall(r"\(([^()]*)\)", body[m.end() - 1:])
    return rows[:n] if len(rows) >= n else rows


def front_plane(src, time):
    """rsa2d's fluid as triangles in the plane (each cell's front face fanned
    from the cell centre) and its obstacle facets (midpoint, unit normal into
    the fluid, length): the wall faces' edges on the front plane."""
    poly = os.path.join(src, "constant", "polyMesh")
    P = np.array([r.split() for r in foam_list(os.path.join(poly, "points"))], dtype=float)
    text = open(os.path.join(poly, "faces")).read()
    body = text[text.index("}", text.index("FoamFile")) + 1:]
    body = body[re.search(r"\d+\s*\(", body).end():]
    F = [np.array(m.split(), dtype=int) for m in re.findall(r"\d+\s*\(([^)]*)\)", body)]
    own = open(os.path.join(poly, "owner")).read()
    own = np.array(re.findall(r"\(([^)]*)\)", own[own.index("}", own.index("FoamFile")) + 1:], re.S)[0].split(), dtype=int)
    C = np.array([r.split() for r in foam_list(os.path.join(src, time, "C"))], dtype=float)
    bnd = open(os.path.join(poly, "boundary")).read()
    patches = {m.group(1): dict(re.findall(r"(\w+)\s+([^;]+);", m.group(2)))
               for m in re.finditer(r"(\w+)\s*\{([^}]*)\}", bnd[bnd.index("}", bnd.index("FoamFile")) + 1:])}
    z0 = P[:, 2].min()
    pe = patches["frontAndBack"]
    tris, nodes = [], [P[:, :2]]
    centre = len(P)
    for f in range(int(pe["startFace"]), int(pe["startFace"]) + int(pe["nFaces"])):
        if np.abs(P[F[f], 2] - z0).max() > 1e-12:
            continue
        cell = own[f]
        k = len(F[f])
        tris += [(centre, F[f][i], F[f][(i + 1) % k]) for i in range(k)]
        nodes.append(C[cell:cell + 1, :2])
        centre += 1
    X = np.concatenate(nodes)
    pw = patches["walls"]
    mids, normals, sizes = [], [], []
    for f in range(int(pw["startFace"]), int(pw["startFace"]) + int(pw["nFaces"])):
        front = [p for p in F[f] if abs(P[p, 2] - z0) < 1e-12]
        a, b = P[front[0], :2], P[front[1], :2]
        e = b - a
        n = np.array([-e[1], e[0]]) / np.linalg.norm(e)
        c = 0.5 * (a + b)
        if np.dot(C[own[f], :2] - c, n) < 0:
            n = -n
        mids.append(c)
        normals.append(n)
        sizes.append(np.linalg.norm(e))
    return X, np.array(tris), np.array(mids), np.array(normals), np.array(sizes), P[:, :2].min(0), P[:, :2].max(0)


def in_fluid(tri, X, lo, hi):
    """Whether each point, wrapped into the periodic box, lies in a fluid triangle."""
    Y = lo + np.mod(X[:, :2] - lo, hi - lo)
    return tri.get_trifinder()(Y[:, 0], Y[:, 1]) >= 0


def from_checkpoint(params, points, dim, c, run, args, threads=4):
    """tracers restarted from a checkpoint whose positions are `points`: a
    two-step run writes the checkpoint, its positions are replaced, the
    restart continues; returns the restart's dumps (in id order) and the id
    of each point."""
    t0 = float(c["time"])
    seed = ["x0=%r" % points[:, 0].mean(), "y0=%r" % points[:, 1].mean(),
            "z0=%r" % (points[:, 2].mean() if dim == 3 else 0.0)]
    base = run + seed + ["init_mode=points_" + "xyz"[:dim], "int_order=1", "dt=%r" % c["dt"],
                         "Nrw=%d" % len(points), "Nrw_max=%d" % len(points), "random=false", "seed=1",
                         "frozen_fields=true", "t_frozen=" + c["time"], "t0=" + c["time"]]
    run_app(TRACERS, params, base, ["T=%r" % (t0 + c["dt"]), "dump_intv=1e9", "stat_intv=1e9",
                                    "checkpoint_intv=%r" % c["dt"], "num_threads=%d" % threads], timeout=1200)
    [pos] = list(params.parent.rglob("Checkpoints/positions.pos"))
    np.savetxt(pos, np.c_[points, np.zeros((len(points), 3 - dim))], fmt="%.17g")
    ids = np.loadtxt(pos.parent / "id.list", dtype=int)
    run_app(TRACERS, params, base, args, ["checkpoint_intv=1e9", "restart_folder=" + str(pos.parent.parent),
                                          "num_threads=%d" % threads], timeout=3600)
    return all_dumps(pos.parent.parent), ids


def lattice(n, off, lo, hi, dim):
    axes = [lo[a] + (np.arange(n) + off[a]) / n * (hi[a] - lo[a]) for a in range(dim)]
    return np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, dim)


def at_rest(params, points, dim, c, run):
    """The fraction of the tracers from `points` at rest after five throughflows."""
    T = float(c["time"]) + c["steps"] * c["dt"]
    D, _ = from_checkpoint(params, points, dim, c, run,
                           ["T=%r" % T, "dump_intv=%r" % (c["steps"] * c["dt"]), "stat_intv=1e9"])
    last = D[max(D)]
    return float(np.mean(np.linalg.norm(last["u"][:, :dim], axis=1) < c["rest"])), last


@needs_openfoam
@pytest.mark.slow
@pytest.mark.parametrize("nodes", ["least_squares", "inverse_distance"])
@pytest.mark.parametrize("wall", ["none", "edge"])
def test_tracers_at_rest_near_the_obstacles_of_rsa2d(tmp_path, nodes, wall):
    """Between obstacles met obliquely, a P1 field's u.n grows linearly above
    each wall facet and tracers reach the walls and stop there (three
    quarters of them in five throughflows); with the near-wall rule u.n
    grows quadratically -- |u.n|/|u.t| ten times larger a decade further
    out -- and a tenth as many stop, those at the obstacles' polygon
    vertices. The nodes do not matter here: a wall facet's cells hold only
    wall points, zero under either rule, and the cell centre. Averaged over
    three shifted lattices, which is what separates the effect from the
    sampling."""
    pytest.importorskip("matplotlib")
    import matplotlib.tri as mtri
    c = NEAR_WALL["rsa2d"]
    src = solved("rsa2d")
    X, tris, mids, normals, sizes, lo, hi = front_plane(src, c["time"])
    tri = mtri.Triangulation(X[:, 0], X[:, 1], tris)
    params = linked(src, c["time"], tmp_path / "c", ["velocity_field=U", "ignore_pressure=true",
                                                      "nodes=" + nodes, "wall_p2=" + wall])
    run = ["mode=openfoam", "Dm=0", "scheme=RK4"]
    # |u.n|/|u.t| at 1e-3 and 1e-2 of each facet's length along its inward normal
    median = {}
    for rel in (1e-3, 1e-2):
        pts = mids + (rel * sizes)[:, None] * normals
        t = float(c["time"]) + 2 * c["dt"]
        D, ids = from_checkpoint(params, pts, 2, c, run, ["T=%r" % t, "dump_intv=%r" % (2 * c["dt"]),
                                                          "stat_intv=1e9"])
        g = D[min(D, key=lambda k: abs(k - t))]
        assert np.abs(g["points"][ids][:, :2] - pts).max() < 1e-12
        u = g["u"][ids][:, :2]
        un = np.einsum("ij,ij->i", u, normals)
        ut = np.linalg.norm(u - un[:, None] * normals, axis=1)
        median[rel] = float(np.median(np.abs(un) / ut))
        shutil.rmtree(tmp_path / "c" / "Tracers")
    print("rsa2d %s %s: |u.n|/|u.t| median %.4g at 1e-3, %.4g at 1e-2" % (nodes, wall, median[1e-3], median[1e-2]))
    if wall == "edge":
        assert abs(median[1e-3] / 2.2e-4 - 1) < 0.02 and abs(median[1e-2] / 2.2e-3 - 1) < 0.02, median
        assert abs(median[1e-2] / median[1e-3] - 10) < 0.1, median
    else:
        assert abs(median[1e-3] / 9.0e-2 - 1) < 0.02 and abs(median[1e-2] / 9.0e-2 - 1) < 0.02, median
    rest = []
    for off in OFFSETS:
        P = lattice(c["lattice"], off, lo, hi, 2)
        P = P[in_fluid(tri, P, lo, hi)]
        assert 3600 < len(P) < 3700
        rest.append(at_rest(params, P, 2, c, run)[0])
        shutil.rmtree(tmp_path / "c" / "Tracers")
    print("rsa2d %s %s: at rest %s, mean %.4f" % (nodes, wall, rest, np.mean(rest)))
    if wall == "edge":
        assert abs(np.mean(rest) - 0.075) < 0.012, rest
    else:
        assert abs(np.mean(rest) - 0.76) < 0.02, rest


@needs_openfoam
@pytest.mark.slow
@pytest.mark.parametrize("wall", ["none", "edge"])
def test_tracers_at_rest_near_the_spheres(tmp_path, wall):
    """The 3D obstacle case: with P1 about one tracer in sixty stops at a
    sphere in five throughflows, with the near-wall rule none."""
    c = NEAR_WALL["spheres3d"]
    params = linked(solved("spheres3d"), c["time"], tmp_path / "c", ["velocity_field=U", "ignore_pressure=true",
                                                                      "wall_p2=" + wall])
    centres = np.array([(0.30, 0.30, 0.30), (0.70, 0.35, 0.65), (0.35, 0.70, 0.70), (0.68, 0.68, 0.32)])
    rest = []
    for off in OFFSETS:
        P = lattice(c["lattice"], off, np.zeros(3), np.ones(3), 3)
        P = P[np.min(np.linalg.norm(P[:, None] - centres[None], axis=2), axis=1) > 0.13 + 1e-9]
        rest.append(at_rest(params, P, 3, c, ["mode=openfoam", "Dm=0", "scheme=RK4"])[0])
        shutil.rmtree(tmp_path / "c" / "Tracers")
    print("spheres3d %s: at rest %s, mean %.4f" % (wall, rest, np.mean(rest)))
    if wall == "edge":
        assert max(rest) <= 0.0005, rest
    else:
        assert abs(np.mean(rest) - 0.016) < 0.003, rest


@needs_openfoam
@pytest.mark.slow
def test_diffusing_tracers_stay_out_of_the_obstacles_of_rsa2d(tmp_path):
    """With Dm > 0 every step is walked from its start and mirrored at the
    walls: after five throughflows no tracer is inside an obstacle (in no
    fluid triangle, the box wrapped), and the noise carries tracers off the
    walls, so no more are at rest than without it."""
    pytest.importorskip("matplotlib")
    import matplotlib.tri as mtri
    c = NEAR_WALL["rsa2d"]
    src = solved("rsa2d")
    X, tris, _, _, _, lo, hi = front_plane(src, c["time"])
    tri = mtri.Triangulation(X[:, 0], X[:, 1], tris)
    params = linked(src, c["time"], tmp_path / "c", ["velocity_field=U", "ignore_pressure=true"])
    P = lattice(c["lattice"], OFFSETS[0], lo, hi, 2)
    P = P[in_fluid(tri, P, lo, hi)]
    # the noise a step about a tenth of a mean cell
    Dm = (0.1 * 0.024655944680018402) ** 2 / (2 * c["dt"])
    still, _ = at_rest(params, P, 2, c, ["mode=openfoam", "Dm=0", "scheme=RK4"])
    shutil.rmtree(tmp_path / "c" / "Tracers")
    moving, last = at_rest(params, P, 2, c, ["mode=openfoam", "Dm=%r" % Dm, "scheme=explicit"])
    print("rsa2d Dm = %.3g: %d tracers, at rest %.4f, %.4f with Dm = 0; %d outside the fluid"
          % (Dm, len(P), moving, still, (~in_fluid(tri, last["points"], lo, hi)).sum()))
    assert in_fluid(tri, last["points"], lo, hi).all()
    assert moving <= still, (moving, still)
