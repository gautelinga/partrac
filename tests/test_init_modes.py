"""The initial states the init_modes place, and the ones they refuse.

pair_<dirs> is one pair of points about (x0, y0, z0): a single edge of length
ds_init centred on the given point, its direction drawn along the named axes
only, whatever Nrw asks for; a centre or an end outside the fluid is refused
rather than moved. points_<dirs> joins its points in order when it samples
one axis, a line, and is a cloud when it samples more; only the line needs
ds_init, ds_max and ds_min. The other modes are refused when they cannot
place what was asked for: a point outside the fluid, too few points inside
it, weights that are zero everywhere, a surface that leaves it, an axis the
domain is flat in, or too few points to step from end to end. Runs are T=0,
so the dump at t=0 is the initial state as the initializer placed it. The
fluid is plane Poiseuille flow, inside where |x| <= 1 and flowing along z, so
a shape is put across the wall at x = 1 by its x0 alone.
"""

import os
import re

import numpy as np
import pytest

from dumps import dump_at
from paths import REPO, app
from runs import copy_example, run_app

PARTRAC = app("partrac")
EXAMPLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")
LINEAR = os.path.join(REPO, "data_example", "linear_flow", "expr_params.dat")
BASE = ("mode=analytic Nrw=10 Nrw_max=100 ds_max=1 ds_min=0.01 ds_init=0.1 Dm=0 dt=0.01 T=0 "
        "int_order=1 dump_intv=0.01 stat_intv=0 checkpoint_intv=1e9 random=false seed=1").split()
CENTRE = np.array([0.2, -0.3, 0.4])

needs = pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")


def refused(r, message):
    """Exit code 2 and message on stderr."""
    assert r.returncode == 2, r.stdout[-1000:] + r.stderr[-1000:]
    assert message in r.stderr, r.stderr


@needs
@pytest.mark.parametrize("dirs", ["x", "y", "z", "xy", "xyz"])
def test_a_pair_is_one_edge_of_length_ds_init_about_the_centre(tmp_path, dirs):
    """Nrw=10, yet two points and one edge between them: ds_init apart, their
    midpoint the centre, and every coordinate off the named axes the centre's.
    The pair is what a separation study starts from, so its length and
    orientation are the experiment's initial condition."""
    params = copy_example(EXAMPLE, tmp_path / "c")
    run_app(PARTRAC, params, BASE, "init_mode=pair_%s" % dirs,
            "x0=%r y0=%r z0=%r" % tuple(CENTRE))
    g = dump_at(tmp_path / "c", 0.0)
    p = g["points"]
    assert p.shape == (2, 3)
    assert sorted(g["edges"].ravel()) == [0, 1]
    assert np.isclose(np.linalg.norm(p[0] - p[1]), 0.1, rtol=1e-12)
    assert np.allclose(p.mean(axis=0), CENTRE, rtol=0, atol=1e-15)
    for k, axis in enumerate("xyz"):
        if axis in dirs:
            assert p[0, k] != CENTRE[k]
        else:
            assert np.array_equal(p[:, k], [CENTRE[k]] * 2)


@needs
@pytest.mark.parametrize("args,message", [
    ("init_mode=pair_xy x0=5", "init_mode pair_xy: the pair centre is not inside the domain"),
    # the ends at x = 0.93 and 1.03, whichever way the pair points
    ("init_mode=pair_x x0=0.98", "init_mode pair_x: the pair is not inside the domain"),
])
def test_a_pair_outside_the_fluid_is_refused(tmp_path, args, message):
    """A centre outside the fluid, or an end outside it: the pair is not
    redrawn elsewhere (that is pairs_<dirs>_<dirs>), so the run stops and says
    which."""
    params = copy_example(EXAMPLE, tmp_path / "c")
    r = run_app(PARTRAC, params, BASE, args, check=False)
    refused(r, message)


@needs
def test_a_point_outside_the_fluid_is_refused(tmp_path):
    """point puts all Nrw particles at (x0, y0, z0); outside the fluid there
    is nothing to place, and a run with no particles would report statistics
    of nothing."""
    params = copy_example(EXAMPLE, tmp_path / "c")
    r = run_app(PARTRAC, params, BASE, "init_mode=point x0=5", check=False)
    refused(r, "init_mode point: (x0, y0, z0) is not inside the domain")


@needs
@pytest.mark.parametrize("args", [
    # every grid cell at x = 5
    "init_mode=points_y init_weight=uniform x0=5",
    # one draw in 500000 inside, at the end of the strip
    "init_mode=randomgaussianstrip_x_y La=1 Lb=0.1 x0=1.499998",
    # the disc's rim just across the wall
    "init_mode=randomgaussiancircle_z_z La=1 Lb=0.1 x0=1.49989",
])
def test_a_sampling_mode_that_cannot_place_nrw_points_is_refused(tmp_path, args):
    """Nrw=100 random points, few or none of whose draws land in the fluid:
    after a million misses in a row the run stops and says how many it placed,
    rather than drawing forever or going on with fewer points than asked for,
    which would change every statistic normalised by Nrw."""
    params = copy_example(EXAMPLE, tmp_path / "c")
    r = run_app(PARTRAC, params, BASE, "Nrw=100", args, check=False, timeout=120)
    mode = args.split()[0].split("=")[1]
    refused(r, "init_mode %s: no points inside the domain in 1000000 draws in a row, with " % mode)
    assert re.search(r"with \d+ of 100 placed", r.stderr), r.stderr


@needs
def test_points_weighted_by_a_velocity_component_that_is_zero_everywhere_are_refused(tmp_path):
    """The flow is along z, so init_weight=ux weighs every cell 0: there is
    nothing to sample by, and the weighted draw would put every point in the
    first cell."""
    params = copy_example(EXAMPLE, tmp_path / "c")
    r = run_app(PARTRAC, params, BASE, "init_mode=points_xy init_weight=ux", check=False)
    refused(r, "init_mode points_xy: init_weight ux has no finite positive total")


@needs
def test_an_ellipsoid_whose_surface_leaves_the_fluid_is_refused(tmp_path):
    """A flat ellipsoid 1.8 across in x centred at x = 0.5: its starting
    tetrahedron is inside the fluid, its surface reaches x = 1.4. A surface
    partly in the solid would be advected by a field that is not there."""
    params = copy_example(EXAMPLE, tmp_path / "c")
    r = run_app(PARTRAC, params, BASE, "Nrw_max=100000 ds_max=0.15 ds_min=0.03",
                "init_mode=ellipsoid_xy La=0.1 Lb=0.9 x0=0.5", check=False)
    refused(r, "init_mode ellipsoid_xy: the ellipsoid is not inside the domain")


@needs
def test_an_ellipsoid_whose_remeshing_cannot_settle_still_starts(tmp_path):
    """ds_min=0.09 and ds_max=0.1: each refinement makes edges that coarsening
    collapses again. The remeshing stops after a bounded number of passes, so
    the run starts rather than hanging."""
    params = copy_example(EXAMPLE, tmp_path / "c")
    r = run_app(PARTRAC, params, BASE, "Nrw_max=100000 ds_max=0.1 ds_min=0.09",
                "init_mode=ellipsoid_xy La=0.5 Lb=0.3", timeout=120)
    assert "Note: the ellipsoid's remeshing stopped after 100 passes" in r.stdout
    g = dump_at(tmp_path / "c", 0.0)
    assert len(g["points"]) >= 4
    assert np.all(np.abs(g["points"][:, 0]) <= 1)


@needs
def test_uniform_along_an_axis_the_domain_is_flat_in_is_refused(tmp_path):
    """A box of zero height in z: uniform_z would put Nrw points on one spot,
    joined by edges of length zero."""
    params = copy_example(LINEAR, tmp_path / "c")
    text = re.sub(r"^z_min=.*$", "z_min=0.0", params.read_text(), flags=re.M)
    text = re.sub(r"^z_max=.*$", "z_max=0.0", text, flags=re.M)
    params.write_text(re.sub(r"^Lz=.*$", "Lz=0.0", text, flags=re.M))
    r = run_app(PARTRAC, params, BASE, "init_mode=uniform_z", check=False)
    refused(r, "init_mode uniform_z: the domain has no extent along z")


@needs
@pytest.mark.parametrize("args,message", [
    ("init_mode=uniform_x", "init_mode uniform needs Nrw of 2 or more"),
    ("init_mode=strip_x La=1", "init_mode strip needs Nrw of 2 or more"),
])
def test_a_line_of_one_point_is_refused(tmp_path, args, message):
    """uniform and strip step from one end to the other in Nrw-1 steps; with
    Nrw=1 there is no step, and the one point's position would be 0/0."""
    params = copy_example(EXAMPLE, tmp_path / "c")
    r = run_app(PARTRAC, params, BASE, "Nrw=1", args, check=False)
    refused(r, message)


@needs
def test_points_along_one_axis_are_joined_in_order(tmp_path):
    """points_x samples a line: its points sorted along x, each edge joining
    two that follow each other. That order is the line's, so the edges are a
    filament that can stretch."""
    params = copy_example(EXAMPLE, tmp_path / "c")
    run_app(PARTRAC, params, BASE, "init_mode=points_x init_weight=uniform")
    g = dump_at(tmp_path / "c", 0.0)
    x = g["points"][:, 0]
    assert np.all(np.diff(x) >= 0)
    edges = np.sort(g["edges"], axis=1)
    assert len(edges) > 0
    assert np.all(edges[:, 1] - edges[:, 0] == 1)


@needs
@pytest.mark.parametrize("dirs", ["xy", "xyz"])
def test_points_over_a_plane_or_a_volume_are_a_cloud(tmp_path, dirs):
    """Points sampled over more than one axis have no order a line could
    follow, so they start with no edges."""
    params = copy_example(EXAMPLE, tmp_path / "c")
    run_app(PARTRAC, params, BASE, "init_mode=points_%s init_weight=uniform" % dirs)
    g = dump_at(tmp_path / "c", 0.0)
    assert len(g["points"]) == 10
    assert "edges" not in g


@needs
@pytest.mark.parametrize("init_mode,missing", [
    ("points_x", ["ds_init", "ds_max", "ds_min"]),
    ("points_xy", []),
    ("points_xyz", []),
])
def test_only_points_along_one_axis_need_edge_lengths(tmp_path, init_mode, missing):
    """points_x joins its points into edges, which read ds_init when joined
    and ds_max and ds_min as they evolve, so the schema asks for all three; a
    cloud over a plane or a volume has no edges and runs without them, rather
    than demanding values nothing reads."""
    params = copy_example(EXAMPLE, tmp_path / "c")
    base = [a for a in BASE if a.split("=")[0] not in ("ds_init", "ds_max", "ds_min")]
    r = run_app(PARTRAC, params, base, "init_mode=%s init_weight=uniform" % init_mode, check=False)
    if missing:
        for key in missing:
            refused(r, "missing required parameter '%s'" % key)
    else:
        assert r.returncode == 0, r.stdout[-1000:] + r.stderr[-1000:]
