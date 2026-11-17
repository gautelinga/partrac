"""Crossing records (cross_intv) of partrac and the tracer apps in plane Poiseuille,
u_z = 1.5 (1 - x^2): a strip along x moves each node straight up z, so the node
at x first passes z = k dn at t = k dn / (1.5 (1 - x^2)), exactly (the motion is
linear in a step), and every plane below the slowest node's last z gets every
half edge once. A line element there, l(0) = (a, b, 0), becomes
l(t) = (a, b, -3 x a t), exactly in a step too (J^2 = 0), so its log stretch and
direction at a crossing are its closed form's, interpolated linearly in the step."""

import os

import numpy as np
import pytest

from dumps import all_dumps
from paths import REPO, app, needs
from runs import (checkpoint_folder, continuous_and_resumed, copy_example, kill_after_checkpoint, run_app,
                  same_dumps)

PARTRAC = app("partrac")
TRACERS = app("tracers")
VECTORS = app("tracervectors")
TENSORS = app("tracertensors")
POISEUILLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")
BASE = ("mode=analytic ds_max=1e9 ds_min=1e-9 refine=false coarsen=false Dm=0 int_order=1 "
        "dt=0.01 stat_intv=0.01 dump_intv=1e9 checkpoint_intv=1e9 random=false seed=1 "
        "init_mode=strip_x La=0.5 x0=0 y0=0 z0=0 Nrw=51 Nrw_max=2000").split()
TRACER_BASE = ("mode=analytic init_mode=points_xy x0=0 y0=0 z0=0 Nrw=40 Nrw_max=40 Dm=0 "
               "int_order=1 dt=0.01 stat_intv=1e9 dump_intv=1e9 checkpoint_intv=1e9 "
               "random=false seed=1").split()
POINT_COLUMNS = "id t k x y z".split()
VECTOR_COLUMNS = POINT_COLUMNS + "lnrho nx ny nz".split()
DT = 0.01


def run(binary, d, base, extra, check=True):
    """binary on a copy of the Poiseuille example in the folder d; base and
    extra its arguments. With check, the run must exit 0."""
    return run_app(binary, copy_example(POISEUILLE, d), base, extra, check=check, timeout=600)


def read_records(d):
    """The one crossing file under d: its records by the column names its first line gives."""
    head, body = next(d.rglob("crossings.bin")).read_bytes().split(b"\n", 1)
    names = head.decode().split("float64 ")[1].split(";")[0].split()
    return dict(zip(names, np.frombuffer(body, np.float64).reshape(-1, len(names)).T))


@needs(PARTRAC)
def test_every_node_records_every_plane_once_at_its_time(tmp_path):
    """Planes z = 0.1 k up to 0.7: 100 half edges each, weighing the strip's 0.5;
    a filter of 50 edges to 20 is logged as such."""
    run(PARTRAC, tmp_path, BASE, ["T=0.5", "cross_n=z", "cross_intv=0.1"])
    r = read_records(tmp_path)
    c, t, weight, k = r["c"], r["t"], r["weight"], r["k"]
    g = all_dumps(tmp_path, raw=True)[0.0]
    x = dict(zip(g["c"][:, 0], g["points"][:, 0]))
    xs = np.array([x[ci] for ci in c])
    assert np.allclose(t, 0.1 * k / (1.5 * (1 - xs ** 2)), rtol=0, atol=1e-12)
    assert k.max() == 7
    for kk in range(1, 8):
        assert (k == kk).sum() == 100 and np.isclose(weight[k == kk].sum(), 0.5, rtol=0, atol=1e-12)

    run(PARTRAC, tmp_path / "filter", BASE, ["T=0.1", "cross_n=z", "cross_intv=0.1", "filter=true",
                                             "filter_target=20", "filter_intv=0.05"])
    log = np.loadtxt(next((tmp_path / "filter").rglob("crossings.bin.filter")), ndmin=2)
    assert log[0, 1:].tolist() == [50, 20] and (log[1:, 1:] == 20).all()


@needs(PARTRAC)
def test_nodes_added_on_the_way_record_from_their_first_step(tmp_path):
    """Refined as it shears, a strip of 160 edges of 0.01 doubles before z = 0.25:
    each plane still weighs 1.6 to within one half edge (a split across a plane);
    skipping new nodes in their first step would lose 2 % of a plane."""
    run(PARTRAC, tmp_path, BASE, ["T=0.5", "La=1.6", "Nrw=161", "cross_n=z", "cross_intv=0.05",
                                  "refine=true", "refine_intv=0.01", "ds_max=0.0103"])
    r = read_records(tmp_path)
    weight, k = r["weight"], r["k"]
    for kk in range(1, 6):
        assert abs(weight[k == kk].sum() - 1.6) <= 0.006, kk


@needs(PARTRAC)
def test_an_axis_exit_plane_and_its_normal_remove_alike(tmp_path):
    """exit_plane=z and exit_plane=n with cross_n=0,0,2 remove the same middle of
    a strip whose ends stay below z = 0.2: the dumps agree bit for bit."""
    common = ["La=1.6", "T=0.2", "Ln=0.2", "filter_intv=0.05", "dump_intv=0.05"]
    for case, extra in (("axis", ["exit_plane=z"]), ("normal", ["exit_plane=n", "cross_n=0,0,2"])):
        run(PARTRAC, tmp_path / case, BASE, common + extra)
    same_dumps(tmp_path / "axis", tmp_path / "normal")
    last = all_dumps(tmp_path / "axis", raw=True)
    assert 0 < len(last[max(last)]["id"]) < 51


@needs(PARTRAC)
def test_particles_record_every_plane_once_at_their_time(tmp_path):
    """Particles without edges, from z = 0: one record each per plane z = 0.1 k they pass,
    up to T + dt (the last step, after the last dump, ends at the final checkpoint)."""
    run(PARTRAC, tmp_path, BASE, ["T=0.5", "init_mode=points_xy", "init_weight=uniform", "cross_n=z",
                                  "cross_intv=0.1"])
    r = read_records(tmp_path)
    assert list(r) == POINT_COLUMNS
    ids, t, k, xc, yc, zc = (r[c] for c in POINT_COLUMNS)
    g = all_dumps(tmp_path, raw=True)[0.0]
    x = dict(zip(g["id"][:, 0], g["points"][:, 0]))
    y = dict(zip(g["id"][:, 0], g["points"][:, 1]))
    assert (g["points"][:, 2] == 0).all()
    xs = np.array([x[i] for i in ids])
    assert np.allclose(t, 0.1 * k / (1.5 * (1 - xs ** 2)), rtol=0, atol=1e-12)
    assert np.array_equal(xc, xs) and np.array_equal(yc, [y[i] for i in ids])
    assert np.allclose(zc, 0.1 * k, rtol=0, atol=1e-12)
    T_end = 0.5 + 0.01
    passed = np.floor(1.5 * (1 - g["points"][:, 0] ** 2) * T_end / 0.1)
    assert passed.sum() > 0
    assert len(ids) == passed.sum() and len(set(zip(ids, k))) == len(ids)


@needs(PARTRAC)
def test_two_normals_are_refused(tmp_path):
    r = run(PARTRAC, tmp_path, BASE, ["T=0.1", "exit_plane=x", "cross_n=z", "Ln=1", "filter_intv=0.1"],
            check=False)
    assert r.returncode != 0 and "names the normal" in r.stdout + r.stderr


@needs(PARTRAC)
def test_a_sheet_is_refused(tmp_path):
    r = run(PARTRAC, tmp_path, BASE, ["T=0.1", "init_mode=sheet_xy", "Lb=0.5", "ds_init=0.05", "cross_n=z",
                                      "cross_intv=0.1"], check=False)
    assert r.returncode != 0 and "a sheet is not supported" in r.stdout + r.stderr


# --- the tracer apps -------------------------------------------------------------

def start(d):
    """The t = 0 dump's datasets, rows by particle id."""
    g = all_dumps(d, raw=True)[0.0]
    ids = g["id"][:, 0]
    out = {}
    for k, v in g.items():
        out[k] = np.empty_like(v)
        out[k][ids] = v
    return out


def line_element(x, l0, t):
    """l(t) = (a, b, c - 3 x a t) of l0 = (a, b, c) at x, rows by particle."""
    l = np.array(l0, dtype=float)
    l[:, 2] += -3.0 * x * l0[:, 0] * t
    return l


@needs(TRACERS, VECTORS, TENSORS)
@pytest.mark.parametrize("binary, columns", [(TRACERS, POINT_COLUMNS), (TENSORS, POINT_COLUMNS),
                                             (VECTORS, VECTOR_COLUMNS)],
                         ids=["tracers", "tracertensors", "tracervectors"])
def test_tracers_record_every_plane_once_where_and_when_they_pass(tmp_path, binary, columns):
    """Points from z = 0: one record per plane z = 0.1 k each passes up to T + dt, at
    its exact time and place; the line element's columns where the app carries one."""
    run(binary, tmp_path, TRACER_BASE, ["T=0.5", "cross_n=z", "cross_intv=0.1"])
    r = read_records(tmp_path)
    assert list(r) == columns
    g = start(tmp_path)
    ids = r["id"].astype(int)
    assert np.array_equal(ids, r["id"]) and (g["points"][:, 2] == 0).all()
    x0, y0 = g["points"][ids, 0], g["points"][ids, 1]
    assert np.allclose(r["t"], 0.1 * r["k"] / (1.5 * (1 - x0 ** 2)), rtol=0, atol=1e-12)
    assert np.array_equal(r["x"], x0) and np.array_equal(r["y"], y0)
    assert np.allclose(r["z"], 0.1 * r["k"], rtol=0, atol=1e-12)
    passed = np.floor(1.5 * (1 - g["points"][:, 0] ** 2) * (0.5 + DT) / 0.1)
    assert passed.sum() > 50
    assert len(ids) == passed.sum() and len(set(zip(ids, r["k"]))) == len(ids)


@needs(VECTORS)
def test_a_line_element_crossing_carries_its_stretch_and_direction(tmp_path):
    """lnrho and n at a crossing: w and rhohat at the step's ends, at the time's
    fraction (rhohat normalized), to round-off; within the step's error of the
    closed form at the crossing time; n a unit vector."""
    run(VECTORS, tmp_path, TRACER_BASE, ["T=0.5", "cross_n=z", "cross_intv=0.05"])
    r = read_records(tmp_path)
    g = start(tmp_path)
    ids = r["id"].astype(int)
    x, l0 = g["points"][ids, 0], g["n"][ids]
    assert np.allclose(np.linalg.norm(l0, axis=1), 1, rtol=0, atol=1e-14) and (g["w"] == 0).all()
    assert np.abs(l0[:, 2]).max() == 0 and np.abs(l0[:, 0]).min() > 0
    tc = r["t"]
    ta = (np.ceil(tc / DT - 1e-9) - 1) * DT
    a = (tc - ta) / DT
    la, lb = line_element(x, l0, ta), line_element(x, l0, ta + DT)
    na, nb = np.linalg.norm(la, axis=1), np.linalg.norm(lb, axis=1)
    w = np.log(na) + a * (np.log(nb) - np.log(na))
    n = la / na[:, None] + a[:, None] * (lb / nb[:, None] - la / na[:, None])
    n /= np.linalg.norm(n, axis=1)[:, None]
    nc = np.stack([r["nx"], r["ny"], r["nz"]], axis=1)
    assert np.allclose(r["lnrho"], w, rtol=0, atol=1e-12)
    assert np.allclose(nc, n, rtol=0, atol=1e-12)
    assert np.allclose(np.linalg.norm(nc, axis=1), 1, rtol=0, atol=1e-14)
    # The closed form at the crossing: w'' <= (3 x a)^2 <= 9, so within 9 dt^2 / 8
    lc = line_element(x, l0, tc)
    assert np.abs(r["lnrho"] - np.log(np.linalg.norm(lc, axis=1))).max() < 9 * DT ** 2 / 8
    assert np.abs(nc - lc / np.linalg.norm(lc, axis=1)[:, None]).max() < 1e-3
    assert r["lnrho"].max() > 0.1


@needs(TRACERS, VECTORS)
@pytest.mark.parametrize("binary", [TRACERS, VECTORS], ids=["tracers", "tracervectors"])
def test_recording_crossings_changes_no_dump(tmp_path, binary):
    """cross_intv = 0, the default, records nothing; recording changes no dump."""
    args = ["T=0.3", "dump_intv=0.1"]
    run(binary, tmp_path / "off", TRACER_BASE, args)
    run(binary, tmp_path / "on", TRACER_BASE, args + ["cross_n=0,0.2,1", "cross_intv=0.05"])
    assert not list((tmp_path / "off").rglob("crossings.bin*"))
    assert not list((tmp_path / "on").rglob("crossings.bin.filter"))   # a strip's log only
    assert len(read_records(tmp_path / "on")["t"]) > 0
    same_dumps(tmp_path / "off", tmp_path / "on")


@needs(TRACERS, VECTORS)
@pytest.mark.parametrize("binary", [TRACERS, VECTORS], ids=["tracers", "tracervectors"])
def test_a_resumed_run_appends_the_records_of_an_uninterrupted_one(tmp_path, binary):
    """Stopped at T = 0.25 and resumed to 0.5: one header, and the records of the
    run to 0.5, bit for bit."""
    cont, split = continuous_and_resumed(binary, POISEUILLE, tmp_path,
                                         [TRACER_BASE, "cross_n=z cross_intv=0.05"], "T=0.25", "T=0.5")
    files = list(split.rglob("crossings.bin"))
    assert len(files) == 1 and files[0].read_bytes().count(b"partrac crossings") == 1
    a, b = read_records(cont), read_records(split)
    assert list(a) == list(b)
    rows = lambda r: np.stack(list(r.values()), axis=1)
    sort = lambda m: m[np.lexsort((m[:, 2], m[:, 0]))]
    assert len(a["t"]) > 50 and np.array_equal(sort(rows(a)), sort(rows(b)))


@needs(TRACERS)
def test_a_resume_starts_a_crossing_file_or_continues_only_its_own(tmp_path):
    """A run without crossings resumed with them: a file with its text line. Resumed
    again with other planes: refused, the file as it was."""
    params = copy_example(POISEUILLE, tmp_path)
    run_app(TRACERS, params, TRACER_BASE, "T=0.15")
    run_app(TRACERS, params, TRACER_BASE, "T=0.3 checkpoint_intv=0.1 cross_n=z cross_intv=0.05",
            "restart_folder=%s" % checkpoint_folder(tmp_path))
    r = read_records(tmp_path)
    assert list(r) == POINT_COLUMNS and r["t"].min() > 0.15 and len(r["t"]) > 20
    f = next(tmp_path.rglob("crossings.bin"))
    before = f.read_bytes()
    p = run_app(TRACERS, params, TRACER_BASE, "T=0.4 cross_n=z cross_intv=0.1",
                "restart_folder=%s" % checkpoint_folder(tmp_path), check=False)
    assert p.returncode != 0 and "holds other records" in p.stdout + p.stderr
    assert f.read_bytes() == before


@needs(TRACERS)
def test_a_new_run_in_the_folder_of_another_starts_its_records_afresh(tmp_path):
    """The same parameters again write into the same folder: the records are
    the new run's alone."""
    params = copy_example(POISEUILLE, tmp_path)
    run_app(TRACERS, params, TRACER_BASE, "T=0.15 cross_n=z cross_intv=0.05")
    first = next(tmp_path.rglob("crossings.bin")).read_bytes()
    run_app(TRACERS, params, TRACER_BASE, "T=0.15 cross_n=z cross_intv=0.05")
    assert next(tmp_path.rglob("crossings.bin")).read_bytes() == first


@needs(TRACERS)
def test_a_resume_refuses_records_its_checkpoint_does_not_account_for(tmp_path):
    """A record file shorter than its checkpoint counted is refused; so is one
    the checkpoint holds no state of (written without it): the file as it was."""
    h5py = pytest.importorskip("h5py")
    params = copy_example(POISEUILLE, tmp_path)
    run_app(TRACERS, params, TRACER_BASE, "T=0.15 cross_n=z cross_intv=0.05")
    f = next(tmp_path.rglob("crossings.bin"))
    resume = ["T=0.3 cross_n=z cross_intv=0.05", "restart_folder=%s" % checkpoint_folder(tmp_path)]
    whole = f.read_bytes()
    f.write_bytes(whole[:-8])
    p = run_app(TRACERS, params, TRACER_BASE, resume, check=False)
    assert p.returncode != 0 and "fewer than the %d its checkpoint counted" % len(whole) in p.stdout + p.stderr
    assert f.read_bytes() == whole[:-8]
    f.write_bytes(whole)
    with h5py.File(checkpoint_folder(tmp_path) / "Checkpoints" / "checkpoint.h5", "r+") as h:
        del h.attrs["crossings"]
    p = run_app(TRACERS, params, TRACER_BASE, resume, check=False)
    assert p.returncode != 0 and "holds no state of" in p.stdout + p.stderr
    assert f.read_bytes() == whole


@needs(TRACERS, PARTRAC)
@pytest.mark.parametrize("binary, base", [(TRACERS, TRACER_BASE), (PARTRAC, BASE)], ids=["tracers", "partrac"])
@pytest.mark.parametrize("args, message", [
    (["cross_intv=0.1"], "needs cross_n"),
    (["cross_n=z", "cross_intv=-0.1"], "cross_intv: '-0.1': intervals cannot be negative"),
    (["cross_n=z", "cross_intv=0.1:1,log0:5"], "cross_intv: 'log0:5': N in logN"),
    (["cross_n=z", "cross_intv=0.1", "cross_origin=0,1"], "cross_origin must be start or x,y,z"),
    (["cross_n=0,0,0"], "cross_n must be x, y, z or nx,ny,nz, not zero"),
    (["cross_n=w"], "cross_n must be x, y, z or nx,ny,nz, not zero"),
    (["cross_n=none", "cross_intv=0.1"], "cross_n must be x, y, z or nx,ny,nz, not zero"),
], ids=["no normal", "negative", "bad schedule", "bad origin", "zero", "no axis", "none"])
def test_planes_are_checked_with_the_parameters(tmp_path, binary, base, args, message):
    """--check reports a bad plane, as a run does, and makes no file."""
    r = run(binary, tmp_path, base, ["T=0.1", "--check"] + args, check=False)
    assert r.returncode != 0 and message in r.stdout + r.stderr
    assert not list(tmp_path.rglob("crossings.bin"))


# --- levels and origins ------------------------------------------------------------

@needs(TRACERS)
def test_the_header_names_the_columns_the_normal_the_levels_and_the_origin(tmp_path):
    run(TRACERS, tmp_path, TRACER_BASE, ["T=0.1", "cross_n=z", "cross_intv=0.05"])
    head = next(tmp_path.rglob("crossings.bin")).read_bytes().split(b"\n", 1)[0].decode()
    assert head == "partrac crossings: float64 id t k x y z; n = (0, 0, 1), cross_intv = 0.05, cross_origin = 0,0,0"
    assert not list(tmp_path.rglob("crossings.bin.levels"))     # a bare number's are k cross_intv


@needs(TRACERS)
def test_a_schedules_levels_are_the_values_it_dumps_at(tmp_path):
    """cross_intv and dump_intv the same schedule: the levels listed are the
    dump times above 0 (snapped to steps of 1e-4), one for one, k = 1, 2, ..."""
    schedule = "0.05:0.2,log5:2,0.3"
    args = ["T=3", "dt=1e-4", "Nrw=2", "Nrw_max=2", "cross_n=z", "cross_intv=" + schedule,
            "dump_intv=" + schedule]
    run(TRACERS, tmp_path, TRACER_BASE, args)
    levels = np.loadtxt(next(tmp_path.rglob("crossings.bin.levels")), ndmin=2)
    dumped = np.array(sorted(t for t in all_dumps(tmp_path, raw=True) if t > 0))
    listed = levels[levels[:, 1] <= 3 + 1e-9]
    assert np.array_equal(listed[:, 0], np.arange(1, len(listed) + 1))
    assert len(listed) == len(dumped) == 12
    assert np.abs(listed[:, 1] - dumped).max() <= 0.5e-4 + 1e-12


@needs(TRACERS, VECTORS)
@pytest.mark.parametrize("binary", [TRACERS, VECTORS], ids=["tracers", "tracervectors"])
def test_from_its_start_each_point_records_its_own_downstream_distance(tmp_path, binary):
    """cross_origin=start, points spread in z: X = z - z0 passes k dn at
    t = k dn / (1.5 (1 - x^2)), whatever z0, at x0 + (0, 0, k dn); a line
    element's lnrho and n are those of the closed form at the step's ends,
    interpolated."""
    run(binary, tmp_path, TRACER_BASE, ["T=0.5", "init_mode=points_xyz", "cross_n=z", "cross_intv=0.05",
                                        "cross_origin=start"])
    r = read_records(tmp_path)
    g = start(tmp_path)
    ids = r["id"].astype(int)
    x0 = g["points"][ids]
    assert np.ptp(g["points"][:, 2]) > 0.5
    assert np.allclose(r["t"], 0.05 * r["k"] / (1.5 * (1 - x0[:, 0] ** 2)), rtol=0, atol=1e-12)
    assert np.allclose(r["z"], x0[:, 2] + 0.05 * r["k"], rtol=0, atol=1e-12)
    assert (r["k"] >= 1).all() and len(set(zip(ids, r["k"]))) == len(ids) > 50
    if binary == VECTORS:
        l0 = g["n"][ids]
        ta = (np.ceil(r["t"] / DT - 1e-9) - 1) * DT
        a = (r["t"] - ta) / DT
        la, lb = line_element(x0[:, 0], l0, ta), line_element(x0[:, 0], l0, ta + DT)
        w = np.log(np.linalg.norm(la, axis=1)) + a * (np.log(np.linalg.norm(lb, axis=1))
                                                      - np.log(np.linalg.norm(la, axis=1)))
        assert np.allclose(r["lnrho"], w, rtol=0, atol=1e-12)


@needs(TRACERS)
def test_an_origin_moves_the_planes(tmp_path):
    """cross_origin=0,0,0.03 from z = 0: the planes at z = 0.03 + 0.1 k, k
    from 0 (X starts at -0.03, below the level 0)."""
    run(TRACERS, tmp_path, TRACER_BASE, ["T=0.5", "cross_n=z", "cross_intv=0.1", "cross_origin=0,0,0.03"])
    r = read_records(tmp_path)
    x0 = start(tmp_path)["points"][r["id"].astype(int)]
    assert r["k"].min() == 0
    assert np.allclose(r["z"], 0.03 + 0.1 * r["k"], rtol=0, atol=1e-12)
    assert np.allclose(r["t"], (0.03 + 0.1 * r["k"]) / (1.5 * (1 - x0[:, 0] ** 2)), rtol=0, atol=1e-12)


@needs(PARTRAC)
def test_a_strip_has_no_start_of_its_own(tmp_path):
    """cross_origin=start is refused for a strip, by --check too."""
    for check in (["--check"], []):
        r = run(PARTRAC, tmp_path / str(len(check)), BASE,
                check + ["T=0.1", "cross_n=z", "cross_intv=0.1", "cross_origin=start"], check=False)
        assert r.returncode != 0 and "a strip's new nodes have no start of their own" in r.stdout + r.stderr


# --- killed and resumed -----------------------------------------------------------

ABC = os.path.join(REPO, "data_example", "abc_flow_unsteady", "expr_params.dat")
KILL_BASE = "Dm=0 int_order=1 scheme=RK4 dt=0.01 stat_intv=1e9 dump_intv=1e9 random=false seed=1"
# In the ABC flow x.n goes back and forth across the planes, and a refining
# strip filters again and again, at random; the Poiseuille strip filters once,
# then its filtering at every checkpoint changes nothing
KILLED = {
    "tracers": (TRACERS, ABC, "mode=analytic init_mode=points_xyz x0=5 y0=5 z0=5 Nrw=300 Nrw_max=300 "
                "T=12 checkpoint_intv=0.3 cross_n=0.3,0.2,1 cross_intv=0.05", 1.0),
    "tracers schedule": (TRACERS, ABC, "mode=analytic init_mode=points_xyz x0=5 y0=5 z0=5 Nrw=300 Nrw_max=300 "
                         "T=12 checkpoint_intv=0.3 cross_n=0.3,0.2,1 cross_intv=0.05:0.5,log10:3,0.2", 1.0),
    "tracervectors start": (VECTORS, ABC, "mode=analytic init_mode=points_xyz x0=5 y0=5 z0=5 Nrw=300 "
                            "Nrw_max=300 T=12 checkpoint_intv=0.3 cross_n=0.3,0.2,1 cross_intv=0.05 "
                            "cross_origin=start", 1.0),
    "tracervectors": (VECTORS, ABC, "mode=analytic init_mode=points_xyz x0=5 y0=5 z0=5 Nrw=300 Nrw_max=300 "
                      "T=12 checkpoint_intv=0.3 cross_n=0.3,0.2,1 cross_intv=0.05", 1.0),
    "strip": (PARTRAC, ABC, "mode=analytic init_mode=strip_x La=1 x0=5 y0=5 z0=5 Nrw=101 Nrw_max=200000 "
              "ds_max=0.01 ds_min=0.002 refine=true coarsen=true refine_intv=0.02 coarsen_intv=0.02 "
              "T=4 checkpoint_intv=0.1 cross_n=0.2,0.1,-1 cross_intv=0.02", 0.8),
    "refiltered strip": (PARTRAC, ABC, "mode=analytic init_mode=strip_x La=1 x0=5 y0=5 z0=5 Nrw=101 "
                         "Nrw_max=200000 ds_max=0.01 ds_min=1e-9 refine=true refine_intv=0.02 filter=true "
                         "filter_target=150 filter_intv=0.1 T=4 checkpoint_intv=0.1 cross_n=0.2,0.1,-1 "
                         "cross_intv=0.02", 1.0),
    "filtered strip": (PARTRAC, POISEUILLE, "mode=analytic init_mode=strip_x La=0.5 x0=0 y0=0 z0=0 Nrw=51 "
                       "Nrw_max=2000 ds_max=1e9 ds_min=1e-9 filter=true filter_target=20 filter_intv=0.05 "
                       "T=40 checkpoint_intv=0.1 cross_n=z cross_intv=0.001", 1.0),
}


def kill_past(binary, params, args, t_kill):
    """kill_after_checkpoint once the record file runs past the bytes the
    checkpoint counted; then the file cut within a record, as a kill mid-write
    leaves it. Returns the checkpoint's crossings_bytes."""
    def past(attrs):
        files = list(params.parent.rglob("crossings.bin"))
        return bool(files) and files[0].stat().st_size > attrs["crossings_bytes"][0]
    attrs = kill_after_checkpoint(binary, params, args, t_kill, past)
    with open(next(params.parent.rglob("crossings.bin")), "ab") as f:
        f.write(b"torn record")
    return attrs["crossings_bytes"]


@needs(TRACERS, VECTORS, PARTRAC)
@pytest.mark.parametrize("case", list(KILLED))
def test_a_run_killed_and_resumed_records_what_the_run_never_stopped_does(tmp_path, case):
    """Killed after a mid-run checkpoint, its record file past the checkpoint
    and cut within a record, and resumed from it: the records and a strip's
    filter log are those of the run never stopped, byte for byte. The resume
    cuts the files back to the checkpoint's lengths and restores the running
    max: in the ABC flow nodes stand back from theirs at the checkpoint, so a
    max afresh would record their planes again. A strip filtering after the
    restart drops the edges it would have, its generator checkpointed."""
    binary, example, args, t_kill = KILLED[case]
    args = [KILL_BASE, args]
    cont = copy_example(example, tmp_path / "cont")
    run_app(binary, cont, args)
    killed = copy_example(example, tmp_path / "killed")
    lengths = kill_past(binary, killed, args, t_kill)
    resume = checkpoint_folder(killed.parent)
    if "ABC" in example:
        h5py = pytest.importorskip("h5py")
        with h5py.File(resume / "Checkpoints" / "checkpoint.h5", "r") as h:
            n = np.array([float(v) for v in args[1].split("cross_n=")[1].split()[0].split(",")])
            o = h["crossings_start"][()] if "start" in case else 0.
            behind = h["crossings_max"][:, 0] - (h["points"][()] - o) @ (n / np.linalg.norm(n))
        assert (behind > 0.02).sum() > 10
    if case == "filtered strip":
        assert lengths[2] == 1      # checkpointed after its step's filtering
    if case == "refiltered strip":
        log = np.loadtxt(next(cont.parent.rglob("crossings.bin.filter")), ndmin=2)
        assert ((log[:, 0] > 2 * t_kill) & (log[:, 1] > log[:, 2])).sum() > 5
    run_app(binary, killed, args, "restart_folder=%s" % resume)
    for name in ("crossings.bin", "crossings.bin.filter", "crossings.bin.levels"):
        a, b = list(cont.parent.rglob(name)), list(killed.parent.rglob(name))
        made = {"crossings.bin": True, "crossings.bin.filter": binary == PARTRAC,
                "crossings.bin.levels": "schedule" in case}[name]
        assert len(a) == len(b) == int(made)
        if a:
            assert a[0].read_bytes() == b[0].read_bytes(), name
    assert len(read_records(cont.parent)["t"]) > 1000
