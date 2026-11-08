"""partrac's crossing records (cross_dn) in plane Poiseuille, u_z = 1.5 (1 - x^2):
a strip along x moves each node straight up z, so the node at x first passes
z = k dn at t = k dn / (1.5 (1 - x^2)), exactly (the motion is linear in a step),
and every plane below the slowest node's last z gets every half edge once."""

import os

import numpy as np

from dumps import all_dumps
from paths import REPO, app, needs
from runs import copy_example, run_app, same_dumps

PARTRAC = app("partrac")
POISEUILLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")
BASE = ("mode=analytic ds_max=1e9 ds_min=1e-9 refine=false coarsen=false Dm=0 int_order=1 "
        "dt=0.01 stat_intv=0.01 dump_intv=1e9 checkpoint_intv=1e9 random=false seed=1 "
        "init_mode=strip_x La=0.5 x0=0 y0=0 z0=0 Nrw=51 Nrw_max=2000").split()


def run(tmp_path, extra):
    return run_app(PARTRAC, copy_example(POISEUILLE, tmp_path), BASE, extra, check=False, timeout=600)


def records(tmp_path, columns=6):
    """c, t, ln rho, weight, length, k; for particles id, t, k."""
    body = next(tmp_path.rglob("crossings.bin")).read_bytes().split(b"\n", 1)[1]
    return np.frombuffer(body, np.float64).reshape(-1, columns).T


@needs(PARTRAC)
def test_every_node_records_every_plane_once_at_its_time(tmp_path):
    """Planes z = 0.1 k up to 0.7: 100 half edges each, weighing the strip's 0.5;
    a filter of 50 edges to 20 is logged as such."""
    r = run(tmp_path, ["T=0.5", "cross_n=z", "cross_dn=0.1"])
    assert r.returncode == 0, r.stdout + r.stderr
    c, t, _, weight, _, k = records(tmp_path)
    g = all_dumps(tmp_path, raw=True)[0.0]
    x = dict(zip(g["c"][:, 0], g["points"][:, 0]))
    xs = np.array([x[ci] for ci in c])
    assert np.allclose(t, 0.1 * k / (1.5 * (1 - xs ** 2)), rtol=0, atol=1e-12)
    assert k.max() == 7
    for kk in range(1, 8):
        assert (k == kk).sum() == 100 and np.isclose(weight[k == kk].sum(), 0.5, rtol=0, atol=1e-12)

    r = run(tmp_path / "filter", ["T=0.1", "cross_n=z", "cross_dn=0.1", "filter=true",
                                  "filter_target=20", "filter_intv=0.05"])
    assert r.returncode == 0, r.stdout + r.stderr
    log = np.loadtxt(next((tmp_path / "filter").rglob("crossings.bin.filter")), ndmin=2)
    assert log[0, 1:].tolist() == [50, 20] and (log[1:, 1:] == 20).all()


@needs(PARTRAC)
def test_nodes_added_on_the_way_record_from_their_first_step(tmp_path):
    """Refined as it shears, a strip of 160 edges of 0.01 doubles before z = 0.25:
    each plane still weighs 1.6 to within one half edge (a split across a plane);
    new nodes skipped in their first step lost 2 % of a plane."""
    r = run(tmp_path, ["T=0.5", "La=1.6", "Nrw=161", "cross_n=z", "cross_dn=0.05",
                       "refine=true", "refine_intv=0.01", "ds_max=0.0103"])
    assert r.returncode == 0, r.stdout + r.stderr
    _, _, _, weight, _, k = records(tmp_path)
    for kk in range(1, 6):
        assert abs(weight[k == kk].sum() - 1.6) <= 0.006, kk


@needs(PARTRAC)
def test_an_axis_exit_plane_and_its_normal_remove_alike(tmp_path):
    """exit_plane=z and exit_plane=n with cross_n=0,0,2 remove the same middle of
    a strip whose ends stay below z = 0.2: the dumps agree bit for bit."""
    common = ["La=1.6", "T=0.2", "Ln=0.2", "filter_intv=0.05", "dump_intv=0.05"]
    for case, extra in (("axis", ["exit_plane=z"]), ("normal", ["exit_plane=n", "cross_n=0,0,2"])):
        r = run(tmp_path / case, common + extra)
        assert r.returncode == 0, r.stdout + r.stderr
    same_dumps(tmp_path / "axis", tmp_path / "normal")
    last = all_dumps(tmp_path / "axis", raw=True)
    assert 0 < len(last[max(last)]["id"]) < 51


@needs(PARTRAC)
def test_particles_record_every_plane_once_at_their_time(tmp_path):
    """Particles without edges, from z = 0: one record each per plane z = 0.1 k they pass,
    up to T + dt (the last step, after the last dump, ends at the final checkpoint)."""
    r = run(tmp_path, ["T=0.5", "init_mode=points_xy", "init_weight=uniform", "cross_n=z", "cross_dn=0.1"])
    assert r.returncode == 0, r.stdout + r.stderr
    assert b"float64 id t k;" in next(tmp_path.rglob("crossings.bin")).read_bytes()
    ids, t, k = records(tmp_path, 3)
    g = all_dumps(tmp_path, raw=True)[0.0]
    x = dict(zip(g["id"][:, 0], g["points"][:, 0]))
    assert (g["points"][:, 2] == 0).all()
    xs = np.array([x[i] for i in ids])
    assert np.allclose(t, 0.1 * k / (1.5 * (1 - xs ** 2)), rtol=0, atol=1e-12)
    T_end = 0.5 + 0.01
    passed = np.floor(1.5 * (1 - g["points"][:, 0] ** 2) * T_end / 0.1)
    assert passed.sum() > 0
    assert len(ids) == passed.sum() and len(set(zip(ids, k))) == len(ids)


@needs(PARTRAC)
def test_two_normals_are_refused(tmp_path):
    r = run(tmp_path, ["T=0.1", "exit_plane=x", "cross_n=z", "Ln=1", "filter_intv=0.1"])
    assert r.returncode != 0 and "names the normal" in r.stdout + r.stderr


@needs(PARTRAC)
def test_a_sheet_is_refused(tmp_path):
    r = run(tmp_path, ["T=0.1", "init_mode=sheet_xy", "Lb=0.5", "ds_init=0.05", "cross_n=z", "cross_dn=0.1"])
    assert r.returncode != 0 and "a sheet is not supported" in r.stdout + r.stderr
