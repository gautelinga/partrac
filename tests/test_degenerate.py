"""What partrac does when a run has nothing left to compute.

Each of these used to end in a crash, a NaN, or a silent success. A run whose
mesh has lost its dimension, or that has no particles at all, is stopped rather
than left to write rows nobody can interpret.
"""

import os
import shutil
import subprocess

import numpy as np
import pytest

from paths import REPO, app

PARTRAC = app("partrac")
POISEUILLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")
TAYLOR_COUETTE = os.path.join(REPO, "data_example", "taylor_couette", "expr_params.dat")

BASE = ("mode=analytic ds_max=1e9 ds_min=1e-9 refine=false coarsen=false Dm=0 "
        "int_order=1 dt=0.01 stat_intv=0.01 dump_intv=1e9 checkpoint_intv=1e9 "
        "random=false seed=1").split()

needs_partrac = pytest.mark.skipif(not os.path.exists(PARTRAC),
                                   reason="partrac is not built")


def run(tmp_path, example, extra):
    tmp_path.mkdir(parents=True, exist_ok=True)
    shutil.copy(example, tmp_path / "expr_params.dat")
    keys = {a.split("=")[0] for a in extra}
    argv = [a for a in BASE if a.split("=")[0] not in keys] + extra
    return subprocess.run([PARTRAC, str(tmp_path / "expr_params.dat")] + argv,
                          capture_output=True, text=True, timeout=600)


def last_row(tmp_path):
    f = list(tmp_path.rglob("tdata_from_t*.dat"))
    assert len(f) == 1
    rows = [l for l in f[0].read_text().splitlines() if l.strip()]
    head = [h.strip() for h in rows[0].lstrip("# ").split("\t") if h.strip()]
    values = [v for v in rows[-1].split("\t") if v.strip()]
    assert len(head) == len(values), "tdata columns do not match its header"
    return dict(zip(head, values))


@needs_partrac
def test_a_strip_coarsened_to_nothing_keeps_its_end_edges(tmp_path):
    # a ds_min larger than the whole strip collapses every interior edge. The
    # two end edges are never collapsed: a tip node carries a single edge, and
    # collapsing it would shorten the strip while its ds0 had nowhere to go.
    r = run(tmp_path, POISEUILLE,
            ["init_mode=strip_x", "La=0.5", "x0=0", "y0=0", "z0=0", "Nrw=10",
             "Nrw_max=2000", "ds_min=100", "coarsen=true", "coarsen_intv=0.01",
             "T=0.1"])
    assert r.returncode == 0, r.stdout + r.stderr
    st = last_row(tmp_path)
    assert int(st["Nrw"]) == 3                        # the two end edges survive
    assert float(st["s0"]) == pytest.approx(0.5)      # and keep the whole s0


@needs_partrac
@pytest.mark.parametrize("ds_min", ["0.02", "0.06", "0.2"])
def test_coarsening_conserves_the_reference_length(tmp_path, ds_min):
    # every collapse hands its ds0 to the neighbours it had, so the sum over
    # edges is invariant however many of them go
    d = tmp_path / ds_min
    r = run(d, POISEUILLE,
            ["init_mode=strip_x", "La=0.5", "x0=0", "y0=0", "z0=0", "Nrw=50",
             "Nrw_max=2000", "ds_min=" + ds_min, "coarsen=true",
             "coarsen_intv=0.01", "T=0.1"])
    assert r.returncode == 0, r.stdout + r.stderr
    st = last_row(d)
    assert float(st["s0"]) == pytest.approx(0.5, rel=1e-9)
    assert 3 <= int(st["Nrw"]) < 50            # it really coarsened


@needs_partrac
def test_a_run_that_loses_every_particle_stops(tmp_path):
    r = run(tmp_path, POISEUILLE,
            ["init_mode=strip_x", "La=0.5", "x0=0", "y0=0", "z0=0", "Nrw=50",
             "Nrw_max=2000", "T=2.0", "stat_intv=0.1", "exit_plane=z",
             "Ln=0.05", "filter_intv=0.05"])
    assert r.returncode != 0
    assert "no particles left" in r.stdout + r.stderr


@needs_partrac
def test_an_initializer_that_places_nothing_stops(tmp_path):
    # the axis is inside the solid inner cylinder, so every particle is rejected
    r = run(tmp_path, TAYLOR_COUETTE,
            ["init_mode=point", "x0=0", "y0=0", "z0=0", "Nrw=1", "Nrw_max=1",
             "T=0.02"])
    assert r.returncode != 0
    assert "no particles left" in r.stdout + r.stderr


@needs_partrac
def test_a_single_particle_reports_zero_spread_not_nan(tmp_path):
    # the sample variance divides by Nrw - 1
    r = run(tmp_path, POISEUILLE,
            ["init_mode=point", "x0=0", "y0=0.2", "z0=0", "Nrw=1", "Nrw_max=1",
             "T=0.02"])
    assert r.returncode == 0, r.stdout + r.stderr
    st = last_row(tmp_path)
    for key, value in st.items():
        assert np.isfinite(float(value)), key
    assert float(st["dx2_mean"]) == 0.0
