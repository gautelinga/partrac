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
    # Nothing moves, so every face must keep dA/dA0 = 1 through coarsening.
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
    # Nothing moves, so tau is the same on every face and must stay so through
    # coarsening, including at t = 0 where every tau is zero.
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
def test_an_exit_plane_on_a_point_cloud_removes_only_what_crosses(tmp_path):
    # A cloud has no edges, so "a node no surviving edge reaches is dead" must
    # not be applied to it: nothing is beyond the plane, so nothing may go.
    r = run(tmp_path, POISEUILLE,
            ["init_mode=uniform_x", "clear_initial_edges=true", "La=0.5",
             "x0=0", "y0=0", "z0=0", "Nrw=50", "Nrw_max=2000", "T=0.05",
             "exit_plane=z", "Ln=100", "filter_intv=0.01"])
    assert r.returncode == 0, r.stdout + r.stderr
    assert int(last_row(tmp_path)["Nrw"]) == 50


@needs_partrac
def test_coarsening_a_point_cloud_keeps_it(tmp_path):
    # a cloud has no edges to collapse, so a pass over it must leave it alone.
    # Deriving "dead" from the edges instead deleted every particle.
    r = run(tmp_path, POISEUILLE,
            ["init_mode=point", "x0=0", "y0=0.1", "z0=0", "Nrw=20",
             "Nrw_max=2000", "coarsen=true", "coarsen_intv=0.01", "T=0.05"])
    assert r.returncode == 0, r.stdout + r.stderr
    assert int(last_row(tmp_path)["Nrw"]) == 20


@needs_partrac
def test_an_exit_plane_on_a_strip_removes_exactly_the_nodes_beyond(tmp_path):
    # and on a strip it is the nodes beyond plus what they orphan, which for a
    # strip is nothing: the node at the cut keeps its inward edge
    def points(case, Ln):
        r = run(tmp_path / case, POISEUILLE,
                ["init_mode=strip_x", "La=0.5", "x0=0", "y0=0", "z0=0",
                 "Nrw=51", "Nrw_max=2000", "T=0", "exit_plane=x",
                 "Ln=%g" % Ln, "filter_intv=0.01"])
        assert r.returncode == 0, r.stdout + r.stderr
        dump = list((tmp_path / case).rglob("data_from_t*.h5"))
        assert len(dump) == 1
        h = h5py.File(dump[0], "r")
        key = sorted(h.keys(), key=float)[0]
        return np.array(h[key + "/points"]), np.array(h[key + "/edges"])

    whole, _ = points("whole", 1e9)               # strip_x is centred on x0
    kept, edges = points("cut", 0.105)
    assert len(kept) < len(whole)                     # it really cut
    assert np.array_equal(kept, whole[whole[:, 0] <= 0.105])
    assert len(edges) == len(kept) - 1                # still one open strip
    assert edges.max() == len(kept) - 1               # and no dangling index


@needs_partrac
def test_an_exit_plane_on_a_sheet_cuts_it(tmp_path):
    # a sheet used to refuse this outright; the nodes beyond go, and with them
    # every face and edge that reached one
    def state(case, Ln):
        r = run(tmp_path / case, POISEUILLE,
                ["init_mode=sheet_xy", "La=0.5", "Lb=0.5", "ds_init=0.05",
                 "x0=0", "y0=0", "z0=0", "Nrw=100", "Nrw_max=200000", "T=0",
                 "exit_plane=x", "Ln=%g" % Ln, "filter_intv=0.01"])
        assert r.returncode == 0, r.stdout + r.stderr
        dump = list((tmp_path / case).rglob("data_from_t*.h5"))
        assert len(dump) == 1
        h = h5py.File(dump[0], "r")
        key = sorted(h.keys(), key=float)[0]
        return tuple(np.array(h[key + "/" + n]) for n in ("points", "faces", "dA0"))

    whole, faces_whole, dA0_whole = state("whole", 1e9)
    kept, faces, dA0 = state("cut", 0.1)
    assert (kept[:, 0] <= 0.1).all()                  # nothing beyond survived
    assert 0 < len(kept) < len(whole)                 # and it is not empty
    assert 0 < len(faces) < len(faces_whole)
    assert faces.max() < len(kept) * 3                # no index left dangling
    assert dA0.min() > 0                              # no face lost its area
    assert dA0.sum() < dA0_whole.sum()


INJECT = ["init_mode=uniform_x", "Nrw=20", "Nrw_max=200000", "inject=true",
          "inject_edges=true", "inject_intv=0.05", "T=0.3"]


@needs_partrac
def test_edge_injection_sweeps_the_inlet_into_a_sheet(tmp_path):
    # inject_edges advects a 1-D inlet and stitches each generation to the last,
    # so the mesh is 2-D from the first injection on. check_dim refused that.
    # The inlet ends sit on the no-slip walls and never leave them, so those two
    # nodes are reused rather than injected again: 18 nodes and 36 faces a
    # generation, not 20 and 38, and the same swept area either way.
    r = run(tmp_path, POISEUILLE, INJECT + ["dump_intv=0.1", "stat_intv=1e9"])
    assert r.returncode == 0, r.stdout + r.stderr
    dump = list(tmp_path.rglob("data_from_t*.h5"))
    assert len(dump) == 1
    h = h5py.File(dump[0], "r")
    keys = sorted(h.keys(), key=float)
    faces = [len(np.array(h[k + "/faces"])) for k in keys if "faces" in h[k]]
    mass = [np.array(h[k + "/dA0"]).sum() for k in keys if "dA0" in h[k]]
    assert len(faces) == 3                            # one generation per dump
    assert faces == [72, 144, 216]                    # 36 per injection
    assert mass[1] == pytest.approx(2 * mass[0], rel=1e-9)   # a strip each time
    assert mass[2] == pytest.approx(3 * mass[0], rel=1e-9)
    dA0 = np.array(h[keys[-1] + "/dA0"]).ravel()
    assert dA0.min() > 0                              # nothing collapsed to a line
    tri = np.array(h[keys[-1] + "/faces"])
    assert all(len(set(row.tolist())) == 3 for row in tri)   # no repeated node


@needs_partrac
def test_an_injected_sheet_carries_a_finite_compressed_time(tmp_path):
    # a zero-area face makes rho = dA/dA0 a 0/0, which spread NaN through tau
    # and into every statistics row after it
    r = run(tmp_path, POISEUILLE,
            INJECT + ["ds_max=0.1", "ds_min=0.02", "refine=true",
                      "refine_intv=0.05", "coarsen=true", "coarsen_intv=0.05",
                      "integrate_tau=true", "tau_intv=0.01", "tau_max=0",
                      "dump_intv=0.3", "stat_intv=0.05"])
    assert r.returncode == 0, r.stdout + r.stderr
    dump = list(tmp_path.rglob("data_from_t*.h5"))
    h = h5py.File(dump[0], "r")
    key = sorted(h.keys(), key=float)[-1]
    for name in ("dA", "dA0", "tau"):
        a = np.array(h[key + "/" + name]).astype(float)
        assert np.isfinite(a).all(), name
    assert np.array(h[key + "/dA0"]).min() > 0
    stats = list(tmp_path.rglob("tdata_from_t*.dat"))
    assert len(stats) == 1
    text = stats[0].read_text().lower()
    assert "nan" not in text and "inf" not in text


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
