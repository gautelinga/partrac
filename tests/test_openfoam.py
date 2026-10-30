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
left unpaired). The cavity loads with either node rule, zero winning at the
lid's corners;
a stamp whose conditions differ from the first's is refused. The load
classes each boundary patch by its velocity condition (no-slip wall, moving
wall, cyclic, other) and says so; a condition built only from its dictionary
is read too. A 2D case may have its empty pair normal to any axis: the
cavity turned onto x and onto y is the same field and traces the same. A
steady case is traced frozen at one stamp: frozen_fields must give that
stamp's field at every step, through partrac. The refusals are
checked by their messages. A phase_field is P1 on the split, its nodes
cellPoint's, each mixed fanned cell's centre moved within [0, 1] so the
volume above 1/2 is the cells' volume fractions', with a gradient of its
own; a scalar condition OpenFOAM cannot build is not data, a velocity one is
refused. A case with cells whose every point is on a no-slip wall is warned
of, mode=openfoam with no reader to load is refused, saying what is missing,
and, marked slow, an installed tree reads the cavity with its own reader.

The tests after those run on cases foam_write builds at test time, each
field a closed form at OpenFOAM's cell centres and the reference: the
cyclic pairing complete on a jittered channel and on a 2D box whose
obstacles cut its sides; Poiseuille second order under refinement; Couette
with a moving wall, a linear pressure and a symmetry plane exact; and the
refusals that need a case of their own (an unwritten value that depends on
the patch, an inverted cell, a zero patch that changes between stamps, a
time without the pressure, times only a processor directory holds) and an
empty cyclic pair, which loads. And on
two small solved fixtures, a cylinder in 2D and a snapped sphere in 3D,
cyclic all round under an oblique body force: |u.n|/|u.t| at the cylinder's
wall the same a decade further out with a P1 field and ten times larger
with the near-wall rule; tracers started at the walls stopping there
without the rule and a tenth as many with it; diffusing tracers kept out of
the cylinder by the reflecting walk.

The reader is built with PARTRAC_ENABLE_OPENFOAM=ON; without it
(build_features.txt says openfoam=off) these tests skip, and fail instead
under PARTRAC_REQUIRE_OPENFOAM=1, as a job that has OpenFOAM runs them."""

import ast
import gzip
import hashlib
import os
import re
import shutil
import subprocess

import numpy as np
import pytest

import foam_write
import openfoam_expected
from cases import once
from dumps import all_dumps, deformation_gradient, dump_at
from paths import BIN, REPO, app, needs
from runs import cells_counts, from_points, interpolation, read_checkpoint, run_app

INTERPOL = app("interpol")
TRACERS = app("tracers")
CAVITY = os.path.join(REPO, "data_example", "openfoam_cavity")
CHANNEL = os.path.join(REPO, "data_example", "openfoam_channel3d")


def built_with_openfoam():
    features = os.path.join(BIN, "build_features.txt")
    return os.path.exists(features) and "openfoam=on" in open(features).read().split()


def _needs_openfoam():
    if built_with_openfoam():
        return needs(INTERPOL)
    if os.environ.get("PARTRAC_REQUIRE_OPENFOAM"):
        return pytest.mark.skipif(False, reason="")
    return pytest.mark.skip(reason="built without PARTRAC_ENABLE_OPENFOAM")


needs_openfoam = _needs_openfoam()
PROBE = ["mode=openfoam", "Nrw=200", "int_order=2", "random=false", "seed=1", "num_threads=1"]


def params_in(d, lines):
    (d / "partrac_params.dat").write_text("".join(l + "\n" for l in lines))
    return d / "partrac_params.dat"


def case_copy(src, d, params=None):
    """The case copied to d, its parameter file replaced by `params` lines if given."""
    shutil.copytree(src, d, ignore=shutil.ignore_patterns("generate_case.sh", "*.py", "split_expected.h5"))
    return d / "partrac_params.dat" if params is None else params_in(d, params)


def probe(params, *args, gradient=False):
    """interpol's points (x, y, z) and velocities, and its output, stdout
    then stderr; with gradient, the velocity gradients J[:, i, j] = du_i/dx_j
    after them."""
    r = run_app(INTERPOL, params, PROBE, list(args), env={"OMP_NUM_THREADS": "1"}, timeout=300)
    g = interpolation(params.parent)
    X, U = np.c_[g["x"], g["y"], g["z"]], np.c_[g["ux"], g["uy"], g["uz"]]
    if gradient:
        J = np.stack([np.c_[g["u%sx" % a], g["u%sy" % a], g["u%sz" % a]] for a in "xyz"], 1)
        return X, U, J, r.stdout + r.stderr
    return X, U, r.stdout + r.stderr


def against_cellpoint(src, d, dim, *args):
    """The largest difference between interpol and OpenFOAM's cellPoint over
    the points both find inside, and how many they are; args after t0 go to
    the parameter file."""
    params = case_copy(src, d, ["velocity_field=U", "pressure_field=p", "wall_upgrade=false",
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
def test_rk4cells_on_the_channel_is_rk4_on_its_linear_field(tmp_path):
    """The channel's field is linear (exact under the least-squares nodes), so
    landing on facets costs RK4cells nothing: its positions and F are RK4's to
    RK4's own error, and it crosses facets of the W12 split without relocating."""
    out = {}
    for scheme in ("RK4cells", "RK4"):
        params = case_copy(CHANNEL, tmp_path / scheme)
        for t in ("1", "2"):
            shutil.copytree(tmp_path / scheme / "0", tmp_path / scheme / t)
        r = run_app(app("tracertensors"), params,
                    ["mode=openfoam", "Nrw=200", "Nrw_max=200", "int_order=1", "Dm=0", "dt=0.02", "t0=1",
                     "T=1.4", "x0=0.5", "y0=0.5", "z0=0.5", "init_mode=points_xyz", "random=false", "seed=1",
                     "dump_intv=0.4", "stat_intv=0.4", "checkpoint_intv=1e9", "scheme=" + scheme], timeout=300)
        g = dump_at(tmp_path / scheme, 1.4)
        out[scheme] = (g["points"], deformation_gradient(g), r.stdout)
    (xc, Fc, log), (x4, F4, _) = out["RK4cells"], out["RK4"]
    c = cells_counts(log)
    assert c and c["crossings"] > 1 and c["relocations"] == 0, log[-2000:]
    assert np.abs(xc - x4).max() < 1e-8 and np.abs(Fc - F4).max() < 1e-8, (np.abs(xc - x4).max(), np.abs(Fc - F4).max())


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


PARTRAC = app("partrac")
FROZEN = ["mode=openfoam", "init_mode=uniform_x", "Nrw=100", "Nrw_max=2000", "ds_max=0.01", "ds_min=0.001",
          "La=0.08", "int_order=1", "Dm=0", "dt=0.002", "x0=0.05", "y0=0.07", "random=false", "seed=1",
          "dump_intv=0.1", "stat_intv=0.1", "checkpoint_intv=1e9", "frozen_fields=true", "num_threads=2"]


@needs_openfoam
def test_frozen_fields_hold_the_stamp_at_t_frozen_whatever_the_time(tmp_path):
    """A steady case is traced as a series frozen at one stamp: run from t = 5
    on stamps at t = 5 and 10 with the fields frozen at either, the tracers
    see that stamp's field at every step, bit for bit what they see on a case
    holding that stamp alone, and not the time weight's blend of the two. A
    time past the stamps freezes the last one."""
    ends = {}
    for tf in (5, 10):
        two = two_stamp_cavity(tmp_path / ("two%d" % tf))
        run_app(PARTRAC, two, FROZEN, ["t0=5", "T=5.2", "t_frozen=%d" % tf], timeout=300)
        one = two_stamp_cavity(tmp_path / ("one%d" % tf))
        shutil.rmtree(tmp_path / ("one%d" % tf) / str(15 - tf))
        run_app(PARTRAC, one, FROZEN, ["t0=%d" % tf, "T=%g" % (tf + 0.2), "t_frozen=%d" % tf], timeout=300)
        a, b = all_dumps(tmp_path / ("two%d" % tf))[5.2], all_dumps(tmp_path / ("one%d" % tf))[tf + 0.2]
        assert len(a["points"]) >= 100
        for key in ("points", "u"):
            assert np.array_equal(a[key], b[key]), (tf, key)
        ends[tf] = a["points"]
    assert np.abs(ends[5] - ends[10]).max() > 1e-4
    # between the stamps: their blend, here three quarters of the t = 10 field
    run_app(PARTRAC, two_stamp_cavity(tmp_path / "mid"), FROZEN, ["t0=5", "T=5.2", "t_frozen=7.5"], timeout=300)
    one = two_stamp_cavity(tmp_path / "quarter")
    shutil.rmtree(tmp_path / "quarter" / "5")
    u = tmp_path / "quarter" / "10" / "U"
    u.write_text(VECTOR.sub(lambda m: "(%r %r %r)" % tuple(0.75 * float(m.group(i)) for i in (1, 2, 3)),
                            u.read_text()))
    run_app(PARTRAC, one, FROZEN, ["t0=10", "T=10.2", "t_frozen=10"], timeout=300)
    a, b = all_dumps(tmp_path / "mid")[5.2], all_dumps(tmp_path / "quarter")[10.2]
    assert np.abs(a["points"] - b["points"]).max() < 1e-12
    run_app(PARTRAC, two_stamp_cavity(tmp_path / "past"), FROZEN, ["t0=5", "T=5.2", "t_frozen=20"], timeout=300)
    assert np.array_equal(all_dumps(tmp_path / "past")[5.2]["points"], ends[10])


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
    m = openfoam_expected.read_mesh(str(case))
    openfoam_expected.geometry(m)
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
    r = run_app(INTERPOL, params, PROBE, ["Nrw=%d" % nrw] + list(args), env={"OMP_NUM_THREADS": "1"}, timeout=300)
    g = interpolation(params.parent)
    G = np.c_[g["dphi_dx"], g["dphi_dy"], g["dphi_dz"]] if "dphi_dx" in g else None
    return np.c_[g["x"], g["y"], g["z"]], g["phi"], G, r.stdout


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


def digest(a):
    """A short name for the array's values."""
    return hashlib.sha1(np.ascontiguousarray(a).tobytes()).hexdigest()[:16]


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
def test_the_phase_fields_nodes_are_cellpoints_whatever_the_node_rule(tmp_path, tmp_path_factory, nodes):
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
    # the same field and points under either rule: probed once
    ref = once(tmp_path_factory, "foam_cavity_alpha_" + digest(P),
               lambda: cellpoint(tmp_path / "c", "alpha.water", "10", P))
    assert np.abs(phi - ref).max() < 1e-12


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
    params = case_copy(CAVITY, tmp_path / "c")
    write_alpha(tmp_path / "c", "10")
    _, _, log = probe(params, "t0=10")
    assert "alpha.water" not in log
    assert "phi" not in interpolation(params.parent)


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
    g = interpolation(params.parent)
    assert np.array_equal(np.c_[g["x"], g["y"], g["z"]], X0)
    assert np.array_equal(np.c_[g["ux"], g["uy"], g["uz"]], U0)


@needs(INTERPOL)
def test_the_mode_without_its_reader_is_refused(tmp_path):
    """mode=openfoam with no reader to load is refused saying what is missing:
    in a build without it, the build option; in one with it, the reader named
    by PARTRAC_OPENFOAM_PLUGIN, here a file that is not there."""
    params = case_copy(CAVITY, tmp_path / "c")
    if not built_with_openfoam():
        refused(params, "mode openfoam needs a build with PARTRAC_ENABLE_OPENFOAM=ON.")
        return
    r = run_app(INTERPOL, params, PROBE, ["t0=10"], env={"PARTRAC_OPENFOAM_PLUGIN": str(tmp_path / "none.so")},
                check=False, timeout=120)
    assert r.returncode == 2, r.stdout[-1000:] + r.stderr
    assert "cannot load the OpenFOAM reader: " in r.stderr and "none.so" in r.stderr, r.stderr


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
    shutil.rmtree(d / "constant" / "polyMesh")
    pid = lambda i, j, k: i + (nx + 1) * (j + 2 * k)
    X = np.array([(i, j, 0.1 * k) for k in range(2) for j in range(2) for i in range(nx + 1)], dtype=float)
    xface = lambda i: [pid(i, 0, 0), pid(i, 1, 0), pid(i, 1, 1), pid(i, 0, 1)]   # normal +x
    yface = lambda i, j: [pid(i, j, 0), pid(i, j, 1), pid(i + 1, j, 1), pid(i + 1, j, 0)]   # normal +y
    zface = lambda i, k: [pid(i, 0, k), pid(i + 1, 0, k), pid(i + 1, 1, k), pid(i, 1, k)]   # normal +z
    rev = lambda f: [f[0]] + f[:0:-1]
    faces = [xface(i + 1) for i in range(nx - 1)]
    owner = list(range(nx - 1))
    neighbour = list(range(1, nx))
    sides = [("left", "wall", [rev(xface(0))], [0]),
             ("bottom", "wall", [rev(yface(i, 0)) for i in range(nx)], list(range(nx))),
             ("lid", "wall", [yface(0, 1)], [0]),
             ("top", "patch", [yface(i, 1) for i in range(1, nx)], list(range(1, nx))),
             ("right", "patch", [xface(nx)], [nx - 1]),
             ("frontAndBack", "empty", [rev(zface(i, 0)) for i in range(nx)] + [zface(i, 1) for i in range(nx)],
              list(range(nx)) * 2)]
    patches = []
    for name, kind, fs, cells in sides:
        patches.append({"name": name, "type": kind, "start": len(faces), "size": len(fs)})
        faces += fs
        owner += cells
    foam_write.write_mesh(d, X, faces, owner, neighbour, patches)
    # the walls noSlip, the patches zeroGradient
    foam_write.write_field(d, "0", "U", {"patches": patches}, np.tile([1., 0., 0.], (nx, 1)))
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


# --- generated cases --------------------------------------------------------
#
# Cases foam_write builds at test time: lattices of hexes, graded toward the
# walls, their interior points jittered, cells removed as obstacles; each
# field a closed form at OpenFOAM's cell centres, which is the reference.


def graded(n, ratio=4.):
    """n + 1 lines on [0, 1], the cells at both ends ratio times smaller than in the middle."""
    assert n % 2 == 0 and n >= 4, "graded takes an even n of at least 4, not %d" % n
    h = ratio ** (np.arange(n // 2) / (n // 2 - 1.))[::-1]
    half = np.r_[0., np.cumsum(0.5 * h / h.sum())]
    return np.r_[half, 1. - half[-2::-1]]


def channel(d, n, x="cyclic", top="wall", height=1., jitter=0.):
    """The channel [0,1]^2 x [0,height] in n^3 hexes, graded toward its bottom
    and top, cyclic in y and, unless x names its two sides, in x."""
    return foam_write.lattice_case(d, np.linspace(0, 1, n + 1), np.linspace(0, 1, n + 1),
                                   height * graded(n), {"x": x, "y": "cyclic", "z": ("wall", top)},
                                   jitter=jitter)


def along_x(m, u):
    """The flow u(z) along x at the cell centres."""
    z = m["C"][:, 2]
    return np.c_[u(z), np.zeros_like(z), np.zeros_like(z)]


def probe_rows(params, time, nrw=2000):
    """interpol's own points (seed 1, one thread) in the case at `time`: X, U, J, log."""
    return probe(params, "t0=" + time, "Nrw=%d" % nrw, gradient=True)


def cyclic_facets(case, dim):
    """The cyclic facets the log counts: a face of k points k - 2 triangles in
    3D, a face an edge in 2D."""
    m = openfoam_expected.read_mesh(case)
    sizes = [len(m["faces"][f]) for p in m["patches"] if p["type"] == "cyclic"
             for f in range(p["start"], p["start"] + p["size"])]
    return sum(k - 2 for k in sizes) if dim == 3 else len(sizes)


@needs_openfoam
@pytest.mark.parametrize("case,split", [("channel", "12"), ("channel", "6"), ("box", "12"), ("sphere", "12")])
def test_every_cyclic_facet_is_paired(tmp_path, case, split):
    """A cyclic facet without its image is refused at load, so a load that
    completes has them all paired; the log counts them, on a jittered channel,
    a 2D box whose obstacles cut its cyclic sides, and the snapped sphere."""
    d = tmp_path / "c"
    lines = ["velocity_field=U", "ignore_pressure=true", "split=" + split]
    if case == "sphere":
        dim, time = 3, NEAR_WALL[SPHERE]["time"]
        params = linked(SPHERE, time, d, lines)
    else:
        time = "0"
        if case == "channel":
            m, dim = channel(d, 8, jitter=0.2), 3
        else:
            solid = np.zeros((12, 12, 1), dtype=bool)
            solid[3:6, 4:7] = solid[0:2, 8:11] = solid[10:12, 1:3] = solid[6:9, 0:2] = True
            m, dim = foam_write.lattice_case(d, np.linspace(0, 1, 13), np.linspace(0, 1, 13), [0, 0.1],
                                             {"x": "cyclic", "y": "cyclic", "z": "empty"}, solid=solid), 2
        foam_write.write_field(d, "0", "U", m, np.tile([1., 0.3, 0.], (len(m["C"]), 1)))
        params = params_in(d, lines)
    X, U, J, log = probe_rows(params, time, 200)
    assert ", %d cyclic (" % cyclic_facets(d, dim) in log, log


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
@pytest.mark.parametrize("nodes", ["least_squares", "inverse_distance"])
def test_poiseuille_on_the_graded_channel_is_second_order(tmp_path, nodes):
    """u = 4 z (1 - z) on the graded channel: refined twice, the velocity's
    error falls by about four and its gradient's by about two, the orders of
    a piecewise linear field."""
    errors = []
    for n in (8, 16):
        d = tmp_path / str(n)
        m = channel(d, n)
        foam_write.write_field(d, "0", "U", m, along_x(m, lambda z: 4 * z * (1 - z)))
        X, U, J, log = probe_rows(params_in(d, ["velocity_field=U", "ignore_pressure=true", "nodes=" + nodes]), "0")
        errors.append(closed_form_errors(X, U, J, lambda z: 4 * z * (1 - z), lambda z: 4 - 8 * z))
    value, gradient = errors[0][0] / errors[1][0], errors[0][1] / errors[1][1]
    assert 3.2 < value < 4.8 and 1.6 < gradient < 2.4, errors


@needs_openfoam
def test_couette_with_a_moving_wall_is_exact(tmp_path):
    """The jittered graded channel, its top moving (movingWallVelocity) and
    u = z: the field is linear, which the least-squares nodes reproduce
    exactly; the top is classed a moving wall, the bottom a no-slip one."""
    d = tmp_path / "c"
    m = channel(d, 8, jitter=0.2)
    foam_write.write_field(d, "0", "U", m, along_x(m, lambda z: z), {"top": ("movingWallVelocity", [1., 0., 0.])})
    X, U, J, log = probe_rows(params_in(d, ["velocity_field=U", "ignore_pressure=true"]), "0")
    assert "128 no-slip wall (bottom), 128 moving wall (top), 512 cyclic" in log, log
    value, gradient = closed_form_errors(X, U, J, lambda z: z, np.ones_like)
    assert value < 1e-12 and gradient < 1e-11, (value, gradient)


@needs_openfoam
def test_a_linear_pressure_is_exact(tmp_path):
    """p = x on the jittered channel, its walls zeroGradient and its x sides
    fixedValue: the least-squares nodes of the pressure reproduce it exactly."""
    d = tmp_path / "c"
    m = channel(d, 8, x=("patch", "patch"), jitter=0.2)
    foam_write.write_field(d, "0", "U", m, along_x(m, lambda z: 4 * z * (1 - z)))
    foam_write.write_field(d, "0", "p", m, m["C"][:, 0], {"left": ("fixedValue", 0.), "right": ("fixedValue", 1.)})
    params = params_in(d, ["velocity_field=U", "pressure_field=p", "nodes=least_squares"])
    X, U, log = probe(params, "t0=0", "Nrw=2000")
    p = interpolation(params.parent)["p"]
    assert np.abs(p - X[:, 0]).max() < 1e-13, np.abs(p - X[:, 0]).max()


@needs_openfoam
def test_a_symmetry_plane_mirrors_exactly(tmp_path):
    """The lower half of the jittered channel, its top a symmetryPlane and
    u = (1, -0.5, 2 (z - 0.5)), which has the plane's symmetry: the
    least-squares nodes, which fit over the cells and their mirror images,
    reproduce it exactly; the plane and the bottom's fixed values are
    classed other, the plane's nodes mirrored."""
    d = tmp_path / "c"
    m = channel(d, 8, top="symmetryPlane", height=0.5, jitter=0.2)
    u = lambda X: np.c_[np.ones(len(X)), -0.5 * np.ones(len(X)), 2 * (X[:, 2] - 0.5)]
    [bottom] = [p for p in m["patches"] if p["name"] == "bottom"]
    foam_write.write_field(d, "0", "U", m, u(m["C"]),
                           {"bottom": ("fixedValue", u(m["fc"][bottom["start"]:bottom["start"] + bottom["size"]]))})
    X, U, J, log = probe_rows(params_in(d, ["velocity_field=U", "ignore_pressure=true", "wall_upgrade=false"]), "0")
    assert "256 other (bottom top)" in log, log
    assert int(re.search(r"(\d+) mirrored", log).group(1)) > 0, log
    Je = np.zeros((3, 3))
    Je[2, 2] = 2.0
    assert np.abs(U - u(X)).max() < 1e-13 and np.abs(J - Je).max() < 1e-11


def small_channel(d, x=("patch", "patch"), times=("0",), bcs=None):
    """The channel in 4^3 hexes, its x sides plain patches, a uniform (1 0 0)
    inside at each of times, the conditions bcs[time] where given; its
    parameter file reads U alone."""
    m = channel(d, 4, x=x)
    for t in times:
        foam_write.write_field(d, t, "U", m, np.tile([1., 0., 0.], (len(m["C"]), 1)), (bcs or {}).get(t))
    return m, params_in(d, ["velocity_field=U", "ignore_pressure=true"])


def refused_at(params, t0, message):
    """interpol's output, refused with message at t0."""
    r = run_app(INTERPOL, params, PROBE, ["t0=" + t0], check=False, timeout=120)
    assert r.returncode == 2, r.stdout[-1000:] + r.stderr
    assert message in r.stderr, r.stderr
    return r.stdout


@needs_openfoam
@pytest.mark.parametrize("condition,entries", [
    ("surfaceNormalFixedValue", {"refValue": "uniform -1"}),
    ("flowRateInletVelocity", {"volumetricFlowRate": "0.5"})])
def test_an_unwritten_value_that_depends_on_the_patch_is_refused(tmp_path, condition, entries):
    """A condition that sets its value without writing it is evaluated on a
    one-cell probe: one whose value follows the patch's normals or areas
    would give the probe's value on every face, so it is refused, asking for
    the value written; noSlip, the same on any patch, is read."""
    _, params = small_channel(tmp_path / "a", bcs={"0": {"left": "noSlip"}})
    X, U, log = probe(params, "t0=0")
    assert len(X) == 200 and "no-slip wall (left bottom top)" in log, log
    _, params = small_channel(tmp_path / "b", bcs={"0": {"left": (condition, None, entries)}})
    refused_at(params, "0", "0/U: patch left is %s, whose value is not written and depends on the patch's faces; "
                            "write its value in the time directory" % condition)


@needs_openfoam
def test_a_cell_the_split_inverts_is_refused(tmp_path):
    """A 2D box of 2 x 2 cells, its middle point pulled towards the far
    corner: that corner's cell is a dart whose centre lies outside it, so the
    centre's fan has an inverted triangle and the load is refused, naming
    the cell and checkMesh."""
    d = tmp_path / "c"
    m = foam_write.lattice_case(d, [0, 0.5, 1], [0, 0.5, 1], [0, 0.1],
                                {"x": ("wall", "wall"), "y": ("wall", "wall"), "z": "empty"})
    pts = d / "constant" / "polyMesh" / "points"
    text = pts.read_text()
    for z in ("0.0", "0.1"):
        assert text.count("(0.5 0.5 %s)" % z) == 1
        text = text.replace("(0.5 0.5 %s)" % z, "(0.95 0.95 %s)" % z)
    pts.write_text(text)
    foam_write.write_field(d, "0", "U", m, np.zeros((len(m["C"]), 3)))
    refused_at(params_in(d, ["velocity_field=U", "ignore_pressure=true", "split=12"]), "0",
               "simplices of the split flat or inverted (under 1e-12 of their cell) in 1 cells, the first in "
               "cell 3; the mesh has warped, flat or inverted cells, which checkMesh -case . reports")


@needs_openfoam
@pytest.mark.parametrize("written", [True, False])
def test_a_patch_zero_at_the_first_stamp_and_not_at_a_later_is_refused(tmp_path, written):
    """Where the inlet meets the no-slip walls, least squares picks the zero
    patches by the first stamp's values: an inlet at rest at t = 1 and moving
    at t = 2 is refused at t = 2, naming the patch. Its value written, or
    unwritten from a table in time, evaluated at each stamp's own instant."""
    if written:
        bcs = {t: {"left": ("fixedValue", [v, 0., 0.])} for t, v in (("1", 0.), ("2", 1.))}
    else:
        table = {"uniformValue": "table ((1 (0 0 0)) (2 (1 0 0)))"}
        bcs = {t: {"left": ("uniformFixedValue", None, table)} for t in ("1", "2")}
    _, params = small_channel(tmp_path / "c", times=("1", "2"), bcs=bcs)
    refused_at(params, "2", "U at 2: patch left is not zero at a point it shares with another patch that fixes "
                            "the value, and was at the first stamp")


@needs_openfoam
def test_a_time_without_the_pressure_is_refused_at_load(tmp_path):
    """p at t = 1 but not at t = 2: refused before the node values are built,
    naming the time."""
    m, _ = small_channel(tmp_path / "c", times=("1", "2"))
    foam_write.write_field(tmp_path / "c", "1", "p", m, m["C"][:, 0])
    params = params_in(tmp_path / "c", ["velocity_field=U", "pressure_field=p"])
    out = refused_at(params, "1", "/2/p: no such file; set pressure_field, or ignore_pressure=true")
    assert "Nodes by" not in out, out


@needs_openfoam
def test_an_empty_cyclic_pair_loads(tmp_path):
    """A cyclic pair with no faces says nothing of the periodic axes: the
    channel with one more, empty, loads as without it."""
    out = []
    for ghosts in (False, True):
        d = tmp_path / str(ghosts)
        m = channel(d, 4)
        end = m["patches"][-1]["start"] + m["patches"][-1]["size"]
        for name, other in (("ghostA", "ghostB"), ("ghostB", "ghostA")) if ghosts else ():
            m["patches"].append({"name": name, "type": "cyclic", "start": end, "size": 0,
                                 "extra": {"neighbourPatch": other, "transformType": "translational",
                                           "separation": "(0 0 %s)" % ("1" if name == "ghostA" else "-1")}})
        foam_write.write_boundary(d, m["patches"])
        foam_write.write_field(d, "0", "U", m, along_x(m, lambda z: 4 * z * (1 - z)))
        out.append(probe(params_in(d, ["velocity_field=U", "ignore_pressure=true"]), "t0=0"))
    (X0, U0, _), (X, U, log) = out
    assert len(X) == 200 and np.array_equal(X, X0) and np.array_equal(U, U0), log


@needs_openfoam
def test_times_only_a_processor_directory_holds_are_refused(tmp_path):
    """A run decomposed after t = 1 and not reconstructed: processor0 holds
    U at t = 2 and 3, the root only at t = 1; refused naming the times."""
    d = tmp_path / "c"
    _, params = small_channel(d, times=("1",))
    for t in ("0", "2", "3"):
        shutil.copytree(d / "1", d / "processor0" / t)
    refused_at(params, "1", "processor0 holds U at times the case's root lacks (2, 3); reconstruct them "
                            "(reconstructPar -case . -newTimes) or remove the processor directories")


@needs_openfoam
@pytest.mark.parametrize("variable", ["WM_PROJECT_DIR", "WM_PROJECT_VERSION"])
def test_another_openfoam_in_the_environment_is_refused(tmp_path, variable):
    """The reader runs the OpenFOAM it was built against: another one sourced
    is refused, naming both, rather than mixing their files."""
    params = case_copy(CAVITY, tmp_path / "c")
    other = str(tmp_path) if variable == "WM_PROJECT_DIR" else "0"
    r = run_app(INTERPOL, params, PROBE, ["t0=10"], env={variable: other}, check=False, timeout=120)
    assert r.returncode == 2, r.stdout[-1000:] + r.stderr
    assert "%s is %s, the build's " % (variable, other) in r.stderr, r.stderr


# --- solved fixtures --------------------------------------------------------
#
# Steady Stokes flows about one obstacle, cyclic on every side, under a body
# force at 20 degrees to x, meshed and solved by OpenFOAM (their
# generate_case.sh): data_example/openfoam_cylinder2d, an O-grid, and
# openfoam_sphere3d, snappyHexMesh's. Curved no-slip walls met obliquely, where
# u.n grows as the square of the distance. Tracers start next to the walls;
# Dm = 0, RK4, dt one mean cell at the fastest point; at rest where |u| is
# under 1e-3 of the mean speed along the force.

CYLINDER = os.path.join(REPO, "data_example", "openfoam_cylinder2d")
SPHERE = os.path.join(REPO, "data_example", "openfoam_sphere3d")


def fixture_constants(path):
    """The literal constants a fixture's generating script assigns, read without running it."""
    out = {}
    with open(path) as f:
        body = ast.parse(f.read()).body
    for node in body:
        if not (isinstance(node, ast.Assign) and len(node.targets) == 1):
            continue
        try:
            value = ast.literal_eval(node.value)
        except ValueError:
            continue
        target = node.targets[0]
        if isinstance(target, ast.Name):
            out[target.id] = value
        elif isinstance(target, ast.Tuple):
            out.update(zip((e.id for e in target.elts), value))
    return out


SPHERE_CASE = fixture_constants(os.path.join(SPHERE, "sphere_case.py"))
SPHERE_CENTRE, SPHERE_RADIUS = np.array(SPHERE_CASE["C"]), SPHERE_CASE["R"]


def near_wall(time, dt, along, throughflows, **more):
    """A fixture's run from what generate_case.sh's measure prints (dt, along the
    mean speed along the force, h the mean cell size): at rest under 1e-3 of
    along, as many steps as the throughflows take."""
    return dict(time=time, dt=dt, rest=1e-3 * along, steps=int(np.ceil(throughflows / along / dt)), **more)


NEAR_WALL = {CYLINDER: near_wall("548", 0.37981, 0.024649, 40, h=0.01981),
             SPHERE: near_wall("219", 0.361877, 0.10918, 5)}


def linked(src, time, d, lines):
    """d with the case's constant, system and time linked in and a parameter file of `lines`."""
    d.mkdir(parents=True)
    for sub in ("constant", "system", time):
        os.symlink(os.path.realpath(os.path.join(src, sub)), d / sub)
    return params_in(d, lines)


def cylinder_wall():
    """The cylinder's wall facets on the front plane (midpoints, unit normals
    into the fluid, lengths) and the polygon they bound, its corners in order."""
    m = openfoam_expected.read_mesh(CYLINDER)
    X = m["X"]
    z0 = X[:, 2].min()
    [wall] = [p for p in m["patches"] if p["name"] == "cylinder"]
    ends = np.array([[X[q, :2] for q in m["faces"][f] if abs(X[q, 2] - z0) < 1e-12]
                     for f in range(wall["start"], wall["start"] + wall["size"])])
    a, b = ends[:, 0], ends[:, 1]
    mids, sizes = 0.5 * (a + b), np.linalg.norm(b - a, axis=1)
    normals = np.c_[a[:, 1] - b[:, 1], b[:, 0] - a[:, 0]] / sizes[:, None]
    centre = a.mean(0)
    normals *= np.sign(np.einsum("ij,ij->i", mids - centre, normals))[:, None]
    corners = np.concatenate([a, b])
    angle = np.arctan2(*(corners - centre).T[::-1])
    return mids, normals, sizes, corners[np.argsort(angle)]


def inside(polygon, X):
    """Whether each point, wrapped into the unit box, is inside the polygon (an odd crossing count)."""
    Y = np.mod(X[:, :2], 1.)[:, None]
    a, b = polygon[None], np.roll(polygon, -1, axis=0)[None]
    crosses = (a[..., 1] > Y[..., 1]) != (b[..., 1] > Y[..., 1])
    at = a[..., 0] + (Y[..., 1] - a[..., 1]) * (b[..., 0] - a[..., 0]) / np.where(crosses, b[..., 1] - a[..., 1], 1.)
    return (crosses & (Y[..., 0] < at)).sum(1) % 2 == 1


def near_cylinder(mids, normals, sizes):
    """Starts above every wall facet at a twentieth to two fifths of its length."""
    return np.concatenate([mids + (r * sizes)[:, None] * normals for r in (0.05, 0.1, 0.2, 0.4)])


def near_sphere(n=400):
    """Starts above the sphere at 0.002 to 0.02, n directions a shell (a
    Fibonacci lattice); the snapped facets lie inside the sphere."""
    k = np.arange(n) + 0.5
    theta, phi = np.arccos(1 - 2 * k / n), np.pi * (1 + 5 ** 0.5) * k
    dirs = np.c_[np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
    return np.concatenate([SPHERE_CENTRE + (SPHERE_RADIUS + d) * dirs for d in (0.002, 0.005, 0.01, 0.02)])


def from_checkpoint(params, points, c, run, args, threads=4):
    """tracers restarted at the fixture's time, frozen there, from a
    checkpoint whose positions are `points` (runs.from_points); returns the
    restart's dumps (in id order) and the id of each point."""
    folder, _ = from_points(TRACERS, params, points, c["dt"],
                            run + ["frozen_fields=true", "t_frozen=" + c["time"], "num_threads=%d" % threads], args,
                            t0=float(c["time"]), timeout=3600)
    return all_dumps(folder), read_checkpoint(folder)["id"][:, 0].astype(int)


def at_rest(params, points, dim, c, run):
    """The fraction of the tracers from `points` at rest after the run, and its last dump."""
    T = float(c["time"]) + c["steps"] * c["dt"]
    D, _ = from_checkpoint(params, points, c, run,
                           ["T=%r" % T, "dump_intv=%r" % (c["steps"] * c["dt"]), "stat_intv=1e9"])
    last = D[max(D)]
    return float(np.mean(np.linalg.norm(last["u"][:, :dim], axis=1) < c["rest"])), last


def rest_once(tmp_path_factory, d, src, points, dim, wall):
    """The fraction at rest with wall_upgrade=wall and Dm = 0, run in d once per session."""
    def run():
        params = linked(src, NEAR_WALL[src]["time"], d, ["velocity_field=U", "ignore_pressure=true",
                                                         "wall_upgrade=" + wall])
        return at_rest(params, points, dim, NEAR_WALL[src], ["mode=openfoam", "Dm=0", "scheme=RK4"])[0]
    return float(once(tmp_path_factory, "foam_rest_%s_%s_%s" % (os.path.basename(src), wall, digest(points)), run))


def rest_with_and_without(tmp_path, tmp_path_factory, src, points, dim):
    """The fraction at rest without the near-wall rule and with it."""
    return {wall: rest_once(tmp_path_factory, tmp_path / wall, src, points, dim, wall) for wall in ("false", "true")}


@needs_openfoam
@pytest.mark.parametrize("nodes", ["least_squares", "inverse_distance"])
@pytest.mark.parametrize("wall", ["false", "true"])
def test_the_near_wall_rule_makes_u_n_grow_as_the_square_at_the_cylinder(tmp_path, nodes, wall):
    """Above each wall facet of the cylinder, at 1e-3 and 1e-2 of its length:
    a P1 field's u.n grows linearly, as u.t does, so |u.n|/|u.t| is the same
    at both; with the near-wall rule it grows as the square, ten times larger
    a decade further out. The nodes do not matter: a wall facet's cells hold
    only wall points, zero under either rule, and the cell centre."""
    c = NEAR_WALL[CYLINDER]
    mids, normals, sizes, _ = cylinder_wall()
    params = linked(CYLINDER, c["time"], tmp_path / "c", ["velocity_field=U", "ignore_pressure=true",
                                                          "nodes=" + nodes, "wall_upgrade=" + wall])
    median = {}
    for rel in (1e-3, 1e-2):
        pts = mids + (rel * sizes)[:, None] * normals
        t = float(c["time"]) + 2 * c["dt"]
        D, ids = from_checkpoint(params, pts, c, ["mode=openfoam", "Dm=0", "scheme=RK4"],
                                 ["T=%r" % t, "dump_intv=%r" % (2 * c["dt"]), "stat_intv=1e9"])
        g = D[min(D, key=lambda k: abs(k - t))]
        assert np.abs(g["points"][ids][:, :2] - pts).max() < 1e-12
        u = g["u"][ids][:, :2]
        un = np.einsum("ij,ij->i", u, normals)
        ut = np.linalg.norm(u - un[:, None] * normals, axis=1)
        median[rel] = float(np.median(np.abs(un) / ut))
        shutil.rmtree(tmp_path / "c" / "Tracers")
    growth = median[1e-2] / median[1e-3]
    assert abs(growth - (10. if wall == "true" else 1.)) < 0.05 * growth, median


@needs_openfoam
def test_tracers_stop_at_the_cylinder_without_the_near_wall_rule(tmp_path, tmp_path_factory):
    """Tracers started next to the cylinder, forty throughflows: with a P1
    field u.n is linear at the wall and a fifth of them reach it and stop;
    with the near-wall rule a tenth as many stop, at the stagnation points."""
    mids, normals, sizes, _ = cylinder_wall()
    rest = rest_with_and_without(tmp_path, tmp_path_factory, CYLINDER, near_cylinder(mids, normals, sizes), 2)
    assert rest["false"] > 0.1 and rest["true"] < 0.1 * rest["false"], rest


@needs_openfoam
def test_tracers_stop_at_the_sphere_without_the_near_wall_rule(tmp_path, tmp_path_factory):
    """The same about the snapped sphere, five throughflows: with P1 some
    stop at the wall, with the near-wall rule a tenth as many."""
    rest = rest_with_and_without(tmp_path, tmp_path_factory, SPHERE, near_sphere(), 3)
    assert rest["false"] > 0.02 and rest["true"] < 0.1 * rest["false"], rest


@needs_openfoam
def test_diffusing_tracers_stay_out_of_the_cylinder(tmp_path, tmp_path_factory):
    """With Dm > 0 every step is walked from its start and mirrored at the
    walls: after forty throughflows no tracer is inside the cylinder (the
    box wrapped), and the noise carries tracers off the wall, so fewer are at
    rest than without it."""
    c = NEAR_WALL[CYLINDER]
    mids, normals, sizes, polygon = cylinder_wall()
    P = near_cylinder(mids, normals, sizes)
    params = linked(CYLINDER, c["time"], tmp_path / "c", ["velocity_field=U", "ignore_pressure=true",
                                                          "wall_upgrade=false"])
    # the noise a step about a tenth of a mean cell
    Dm = (0.1 * c["h"]) ** 2 / (2 * c["dt"])
    moving, last = at_rest(params, P, 2, c, ["mode=openfoam", "Dm=%r" % Dm, "scheme=explicit"])
    still = rest_once(tmp_path_factory, tmp_path / "still", CYLINDER, P, 2, "false")
    assert not inside(polygon, last["points"]).any()
    assert still > 0.1 and moving < still, (moving, still)
