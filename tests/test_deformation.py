"""Stretching of lines and sheets, against an analytic answer.

Plane Poiseuille is axial and depends only on x, so a particle never moves in x
and its velocity is constant along its own path: z(T) = u_z(x)*T exactly,
whatever the timestep. The deformation of a material line or sheet laid in that
flow therefore has a closed form to check against.
"""

import os
import shutil
import subprocess

import numpy as np
import pytest

from paths import REPO, app

PARTRAC = app("partrac")
EXAMPLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")

U_INF, R, LA, T = 1.0, 1.0, 1.0, 0.5

BASE = ("mode=analytic Nrw=200 Nrw_max=5000 ds_max=1e9 ds_min=1e-9 "
        "refine=false coarsen=false Dm=0 int_order=1 dt=0.01 "
        "dump_intv=1e9 checkpoint_intv=1e9 random=false seed=1").split()


def run(tmp_path, extra):
    d = tmp_path / "case"
    d.mkdir(parents=True)
    shutil.copy(EXAMPLE, d / "expr_params.dat")
    keys = {a.split("=")[0] for a in extra}
    base = [a for a in BASE if a.split("=")[0] not in keys]
    r = subprocess.run([PARTRAC, str(d / "expr_params.dat")] + base
                       + ["T=%g" % T, "stat_intv=%g" % T] + extra,
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr
    f = list(d.rglob("tdata_from_t*.dat"))
    assert len(f) == 1
    lines = [l for l in f[0].read_text().splitlines() if l.strip()]
    head = [h.strip() for h in lines[0].lstrip("# ").split("\t") if h.strip()]
    last = [v for v in lines[-1].split("\t") if v.strip()]
    assert len(head) == len(last), "tdata columns do not match its header"
    return dict(zip(head, last))


def analytic_line(n=200):
    """Length of a line initially along x, after time T."""
    x = np.linspace(-LA / 2, LA / 2, n)
    z = 1.5 * U_INF * (1 - (x / R) ** 2) * T
    return np.sqrt(np.diff(x) ** 2 + np.diff(z) ** 2).sum(), np.abs(np.diff(x)).sum()


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_line_length_matches_the_analytic_value(tmp_path):
    st = run(tmp_path, ["init_mode=strip_x", "La=%g" % LA])
    s_exp, s0_exp = analytic_line()
    assert float(st["s0"]) == pytest.approx(s0_exp, rel=1e-6)
    assert float(st["s"]) == pytest.approx(s_exp, rel=1e-4)
    # it really did stretch, and the elongation is finite
    assert float(st["s"]) > float(st["s0"])
    assert np.isfinite(float(st["logelong_wmean"]))
    assert float(st["logelong_wmean"]) > 0


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_sheet_area_grows_and_is_finite(tmp_path):
    st = run(tmp_path, ["init_mode=sheet_xy", "La=%g" % LA, "Lb=%g" % LA,
                        "ds_init=0.1"])
    assert float(st["A0"]) > 0
    assert float(st["A"]) > float(st["A0"])
    assert np.isfinite(float(st["logelong_wmean"]))
    assert float(st["logelong_wmean"]) > 0


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_a_line_reports_length_and_a_sheet_reports_area(tmp_path):
    # the two were labelled the other way round
    line = run(tmp_path / "l", ["init_mode=strip_x", "La=%g" % LA])
    sheet = run(tmp_path / "s", ["init_mode=sheet_xy", "La=%g" % LA,
                                 "Lb=%g" % LA, "ds_init=0.1"])
    assert "s" in line and "s0" in line and "A" not in line
    assert "A" in sheet and "A0" in sheet and "s" not in sheet


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
def test_refinement_tracks_the_curve_better(tmp_path):
    # a polyline understates the arc length; refining should close the gap, and
    # the split edges must keep their reference length consistent
    exact, _ = analytic_line(20001)
    coarse = ["init_mode=strip_x", "La=%g" % LA, "Nrw=40", "Nrw_max=20000"]
    plain = run(tmp_path / "p", coarse + ["ds_max=1e9", "refine=false"])
    fine = run(tmp_path / "r", coarse + ["ds_max=0.01", "refine=true",
                                         "refine_intv=0.05"])
    assert int(fine["Nrw"]) > int(plain["Nrw"])          # it really refined
    e_plain = abs(float(plain["s"]) - exact) / exact
    e_fine = abs(float(fine["s"]) - exact) / exact
    assert e_fine < e_plain
    assert e_fine < 1e-4
