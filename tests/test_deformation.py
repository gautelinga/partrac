"""Stretching of material lines and sheets in partrac, against a closed form.

The flow is plane Poiseuille, u = (0, 0, 1.5 U (1 - (x/R)^2)). It is axial and
depends only on x, so a particle never moves in x and keeps its velocity along
its own path: z(T) = z0 + u_z(x0) T exactly, whatever the time step. A line laid
along x therefore becomes the parabola z = u_z(x) T, whose arc length is known,
and a sheet laid in the x-y plane is sheared into a surface of larger area.
"""

import os

import numpy as np
import pytest

from dumps import all_dumps, read_stats
from paths import REPO, app
from runs import copy_example, run_app

PARTRAC = app("partrac")
EXAMPLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")

# the example's u_inf and half-width, the line length La, and the end time
U_INF, R, LA, T = 1.0, 1.0, 1.0, 0.5

BASE = ("mode=analytic Nrw=200 Nrw_max=5000 ds_max=1e9 ds_min=1e-9 "
        "refine=false coarsen=false Dm=0 int_order=1 dt=0.01 "
        "dump_intv=1e9 checkpoint_intv=1e9 random=false seed=1").split()


def run(tmp_path, extra):
    """Run partrac to T with `extra` overriding BASE; return the last statistics row by column."""
    d = tmp_path / "case"
    run_app(PARTRAC, copy_example(EXAMPLE, d), BASE, ["T=%g" % T, "stat_intv=%g" % T], extra)
    return {k: v[-1] for k, v in read_stats(d).items()}


def analytic_line(n=200):
    """Length at T of an n-node polyline initially along x, and its initial length."""
    x = np.linspace(-LA / 2, LA / 2, n)
    z = 1.5 * U_INF * (1 - (x / R) ** 2) * T
    return np.sqrt(np.diff(x) ** 2 + np.diff(z) ** 2).sum(), np.abs(np.diff(x)).sum()


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_line_length_matches_the_analytic_value(tmp_path):
    """A 200-node line along x reports the initial length s0 and the stretched
    length s of the same polyline in the exact flow, under length columns. If
    this fails, the length and log-elongation statistics of every line run are
    wrong."""
    st = run(tmp_path, ["init_mode=strip_x", "La=%g" % LA])
    # the same node count as the run, so both measure the same polyline
    s_exp, s0_exp = analytic_line()
    assert st["s0"] == pytest.approx(s0_exp, rel=1e-6)
    assert st["s"] == pytest.approx(s_exp, rel=1e-4)
    assert st["s"] > st["s0"]
    assert np.isfinite(st["logelong_wmean"])
    assert st["logelong_wmean"] > 0
    # a line writes length columns, never area columns
    assert "A" not in st and "A0" not in st


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_sheet_area_grows_and_is_finite(tmp_path):
    """A flat sheet sheared by the flow reports an area A larger than its initial
    area A0 and a finite, positive mean log-elongation, under area columns, so
    sheet statistics measure the stretching rather than a zero or a NaN, and a
    postprocessing script reads the quantity its label names."""
    st = run(tmp_path, ["init_mode=sheet_xy", "La=%g" % LA, "Lb=%g" % LA,
                        "ds_init=0.1"])
    assert st["A0"] > 0
    assert st["A"] > st["A0"]
    assert np.isfinite(st["logelong_wmean"])
    assert st["logelong_wmean"] > 0
    # a sheet writes area columns, never length columns
    assert "s" not in st and "s0" not in st


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_refinement_tracks_the_curve_better(tmp_path):
    """Refining a coarse line brings its length closer to the arc length of the
    exact parabola. The split edges must keep their reference length consistent,
    otherwise refinement would add length that the flow never produced."""
    # a polyline understates the arc length, so a 20001-node polyline stands in
    # for the exact curve
    exact, _ = analytic_line(20001)
    coarse = ["init_mode=strip_x", "La=%g" % LA, "Nrw=40", "Nrw_max=20000"]
    plain = run(tmp_path / "p", coarse + ["ds_max=1e9", "refine=false"])
    fine = run(tmp_path / "r", coarse + ["ds_max=0.01", "refine=true",
                                         "refine_intv=0.05"])
    assert fine["Nrw"] > plain["Nrw"]          # it really refined
    e_plain = abs(plain["s"] - exact) / exact
    e_fine = abs(fine["s"] - exact) / exact
    assert e_fine < e_plain
    assert e_fine < 1e-4


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
@pytest.mark.parametrize("scheme,int_order", [("explicit", 1), ("explicit", 2),
                                              ("RK4", 1)])
def test_plane_poiseuille_advects_only_along_the_axis(tmp_path, scheme, int_order):
    """u is axial and depends only on x, so every scheme and order keeps each
    particle's x and y to round-off while it moves along z. An integrator that
    mixes up components would carry particles across streamlines."""
    d = tmp_path / "axis"
    run_app(PARTRAC, copy_example(EXAMPLE, d),
            "mode=analytic init_mode=uniform_x Nrw=100 Nrw_max=5000 ds_max=0.4 ds_min=0.1 "
            "Dm=0 dt=0.01 dump_intv=0.05 stat_intv=0.05 checkpoint_intv=0.05 random=false seed=3",
            ["T=0.3", "scheme=" + scheme, "int_order=%d" % int_order])
    pytest.importorskip("h5py")
    dumps = all_dumps(d, raw=True)
    first, last = dumps[min(dumps)]["points"], dumps[max(dumps)]["points"]
    assert np.abs(last[:, 0] - first[:, 0]).max() < 1e-12
    assert np.abs(last[:, 1] - first[:, 1]).max() < 1e-12
    assert np.abs(last[:, 2] - first[:, 2]).max() > 1e-3  # it did move, so the checks above are not vacuous
