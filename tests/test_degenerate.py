"""What partrac does when a run has nothing left to compute.

Each of these used to end in a crash, a NaN, or a silent success. A run whose
mesh has lost its dimension, or that has no particles at all, is stopped rather
than left to write rows nobody can interpret.
"""

import os
import shutil
import subprocess

import h5py
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
@pytest.mark.parametrize("ds_init", ["0.1", "0.05", "0.025"])
def test_coarsening_a_flat_sheet_keeps_every_elongation(tmp_path, ds_init):
    # nothing moves, so every face must still read dA/dA0 = 1 however much the
    # mesh is coarsened. Sharing the removed dA0 out by clamped relative weights
    # left an rms error of 18% here that refining the mesh did not reduce.
    d = tmp_path / ds_init
    d.mkdir(parents=True)
    (d / "expr_params.dat").write_text(
        open(POISEUILLE).read().replace("u_inf=1.0", "u_inf=0.0"))
    ds = float(ds_init)

    def areas(ds_min):
        case = d / ds_min
        case.mkdir(parents=True)
        (case / "expr_params.dat").write_text((d / "expr_params.dat").read_text())
        r = subprocess.run(
            [PARTRAC, str(case / "expr_params.dat")]
            + [a for a in BASE
               if a.split("=")[0] not in {"ds_min", "dump_intv", "coarsen"}]
            + ["init_mode=sheet_xy", "La=0.5", "Lb=0.5", "ds_init=" + ds_init,
               "x0=0", "y0=0", "z0=0", "Nrw=100", "Nrw_max=200000",
               "ds_min=" + ds_min, "coarsen=true", "coarsen_intv=1e9",
               "T=0.01", "dump_intv=0.01"],
            capture_output=True, text=True, timeout=600)
        assert r.returncode == 0, r.stdout + r.stderr
        dump = list(case.rglob("data_from_t*.h5"))
        assert len(dump) == 1
        h = h5py.File(dump[0], "r")
        key = sorted(h.keys(), key=float)[0]
        return (np.array(h[key + "/dA"]).ravel(),
                np.array(h[key + "/dA0"]).ravel())

    dA_ref, _ = areas("1e-12")                      # nothing collapses
    dA, dA0 = areas("%g" % (1.2 * ds))

    assert len(dA) < 0.6 * len(dA_ref)              # it really coarsened
    assert dA0.min() > 0                            # and never through zero
    assert np.max(np.abs(dA / dA0 - 1)) < 1e-12     # every elongation kept
    assert dA0.sum() == pytest.approx(dA_ref.sum(), rel=1e-9)   # mass conserved


@needs_partrac
def test_coarsening_a_sheet_of_uniform_tau_keeps_it_uniform(tmp_path):
    # Nothing moves, so tau grows identically on every face, and coarsening
    # must not separate them. On a flat sheet the dA0 share-out is exact to
    # the bit, so the taus stay identical and a collapse leaves them alone --
    # including the first one, at t = 0, where every tau is zero and sharing
    # dA0/sqrt(tau) would divide by it. This pins that path, and the property.
    d = tmp_path
    d.mkdir(parents=True, exist_ok=True)
    (d / "expr_params.dat").write_text(
        open(POISEUILLE).read().replace("u_inf=1.0", "u_inf=0.0"))

    def taus(ds_min):
        case = d / ds_min
        case.mkdir(parents=True)
        (case / "expr_params.dat").write_text((d / "expr_params.dat").read_text())
        r = subprocess.run(
            [PARTRAC, str(case / "expr_params.dat")]
            + [a for a in BASE if a.split("=")[0]
               not in {"ds_min", "dump_intv", "coarsen", "dt"}]
            + ["init_mode=sheet_xy", "La=0.5", "Lb=0.5", "ds_init=0.05",
               "x0=0", "y0=0", "z0=0", "Nrw=100", "Nrw_max=200000",
               "integrate_tau=true", "tau_intv=0.001", "tau_max=0",
               "dt=0.001", "T=0.01", "dump_intv=0.01",
               "ds_min=" + ds_min, "coarsen=true", "coarsen_intv=0.005"],
            capture_output=True, text=True, timeout=600)
        assert r.returncode == 0, r.stdout + r.stderr
        dump = list(case.rglob("data_from_t*.h5"))
        assert len(dump) == 1
        h = h5py.File(dump[0], "r")
        key = sorted(h.keys(), key=float)[-1]
        return np.array(h[key + "/tau"]).ravel()

    tau_ref = taus("1e-12")                         # nothing collapses
    tau = taus("0.06")

    assert len(tau) < 0.6 * len(tau_ref)            # it really coarsened
    assert tau_ref.min() > 0                        # and tau really ran
    assert tau_ref == pytest.approx(tau_ref[0], rel=1e-12)
    assert tau == pytest.approx(tau_ref[0], rel=1e-12)


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
