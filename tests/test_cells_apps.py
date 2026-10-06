"""scheme=RK4cells in partrac and filaments: strips, sheets and pairs.

These apps create and move particles outside the step: refinement puts a
node between two, injection adds the inlet's nodes, coarsening merges two,
resizing pulls a pair's node in, outside=reinject moves a stuck node or edge.
Each operation leaves the particle in the cell it located, or in none (-1),
and RK4cells starts it there or locates it afresh. A cell left over from
before the move would hold the particle no longer, and RK4cells would relocate
it, which the relocations it prints count: here they stay at zero.

Every operation runs under RK4cells as under RK4 at a fine dt (same events,
same particles, positions within 1e-5, mesh invariants kept); an analytic
field is RK4 bit for bit; restarts and the thread count change no bit; Dm > 0
and a field without cells are refused.
"""

import os
import re

import numpy as np
import pytest

from cases import ABC, args_for
from dumps import all_dumps
from paths import app
from runs import (cells_counts, checkpoint_file, checkpoint_folder, copy_case, copy_example, read_checkpoint,
                  run_app, same_checkpoint, same_dumps, same_stats)

PARTRAC = app("partrac")
FILAMENTS = app("filaments")

needs_apps = pytest.mark.skipif(not all(os.path.exists(a) for a in (PARTRAC, FILAMENTS)),
                                reason="partrac or filaments is not built")

# dt and the intervals exact in binary; T and an interval for each case
COMMON = "Dm=0 int_order=1 dt=0.00390625 stat_intv=0.125 checkpoint_intv=1e9 random=false seed=1 verbose=true"

# On the periodic P2 tet case, u = (sin 2 pi y, sin 2 pi z, sin 2 pi x)
TET_CASES = {
    # an inlet sweeping a sheet, refined and coarsened
    "sheet": (PARTRAC, "mode=tet init_mode=strip_x La=0.5 x0=0.5 y0=0.53 z0=0.47 Nrw=11 Nrw_max=20000 "
              "inject=true inject_edges=true inject_intv=0.0625 refine=true refine_intv=0.0625 ds_max=0.05 "
              "coarsen=true coarsen_intv=0.0625 ds_min=0.01 T=0.5 dump_intv=0.125"),
    # a strip refined, coarsened, filtered and cut at an exit plane
    "strip": (PARTRAC, "mode=tet init_mode=strip_y La=0.6 x0=0.5 y0=0.5 z0=0.47 Nrw=31 Nrw_max=20000 "
              "refine=true refine_intv=0.0625 ds_max=0.04 coarsen=true coarsen_intv=0.0625 ds_min=0.0199 "
              "filter=true filter_intv=0.125 filter_target=80 exit_plane=x Ln=0.7 T=1.0 dump_intv=0.125"),
    # pairs halved back to ds_init
    "doublings": (FILAMENTS, "mode=tet init_mode=pairs_xyz x0=0.5 y0=0.5 z0=0.5 Nrw=400 Nrw_max=400 "
                  "ds_init=0.01 resize=doublings resize_target=ds_init resize_intv=0.0625 T=1.0 dump_intv=0.25"),
    # pairs rescaled to ds_max
    "rescale": (FILAMENTS, "mode=tet init_mode=pairs_xyz x0=0.5 y0=0.5 z0=0.5 Nrw=400 Nrw_max=400 "
                "ds_init=0.01 ds_max=0.02 resize=rescale resize_intv=0.0625 T=1.0 dump_intv=0.25"),
}

# What each case's run must have done (verbose lines), and the longest edge at a dump
CHECKS = {
    "sheet": ([r"Added [1-9]\d* nodes", r"Added [1-9]\d* edges", r"Removed [1-9]\d* edges"], 0.05),
    "strip": ([r"Added [1-9]\d* edges", r"Removed [1-9]\d* edges", r"Filtered edges",
               r"Removed [1-9]\d* nodes that were beyond"], 0.04),
    "doublings": ([r"Resized edges"], 0.01),
    "rescale": ([r"Resized edges"], 0.02),
}
EVENT = r"(Added|Removed|Filtered|Resized)"


def tet_run(mesh_dir, root, name, case, *args):
    """A case on a copy of the tet case under root/name; (folder, stdout)."""
    binary, own = TET_CASES[case]
    d = copy_case(mesh_dir("tet"), root / name)
    r = run_app(binary, d / "dolfin_params.dat", COMMON, own, *args)
    return d, r.stdout


def matched(x, y):
    """The largest distance from a point of x or y to the nearest of the other's: refinement
    may number two new nodes the other way round where two edges are as long to round-off."""
    from scipy.spatial import cKDTree
    return max(cKDTree(y).query(x)[0].max(), cKDTree(x).query(y)[0].max())


def check_mesh(g, ds_max):
    """A dump's curve or sheet: edges and faces on existing distinct nodes, none longer than ds_max."""
    x = g["points"]
    n = len(x)
    if "faces" in g:
        f = g["faces"][:, :3].astype(np.int64)
        assert (f < n).all() and (f[:, 0] != f[:, 1]).all() and (f[:, 1] != f[:, 2]).all() \
            and (f[:, 0] != f[:, 2]).all()
        assert (g["dA"] > 0).all()
        e = np.unique(np.sort(np.vstack([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]]), axis=1), axis=0)
    else:
        e = g["edges"].astype(np.int64)
        assert (e < n).all() and (e[:, 0] != e[:, 1]).all()
    assert np.isfinite(x).all()
    assert np.linalg.norm(x[e[:, 0]] - x[e[:, 1]], axis=1).max() <= ds_max * (1 + 1e-9)


@needs_apps
@pytest.mark.parametrize("case", list(TET_CASES))
def test_remeshing_and_resizing_under_rk4cells_as_under_rk4(mesh_dir, tmp_path, case):
    """At dt |u| = h/25 RK4cells and RK4 remesh or resize alike: the same
    events, the same particles (by id) at every dump, positions within 1e-5,
    the mesh kept (edges no longer than ds_max after refinement, faces with
    area), and no particle relocated: every node made or moved has its cell."""
    out = {s: tet_run(mesh_dir, tmp_path, s, case, "scheme=" + s) for s in ("RK4cells", "RK4")}
    (dc, sc), (d4, s4) = out["RK4cells"], out["RK4"]
    patterns, longest = CHECKS[case]
    for pattern in patterns:
        assert re.search(pattern, sc), pattern
    assert [l for l in sc.splitlines() if re.match(EVENT, l)] == [l for l in s4.splitlines() if re.match(EVENT, l)]
    assert cells_counts(sc).get("relocations") == 0.0, sc
    a, b = all_dumps(dc), all_dumps(d4)
    raw, raw4 = all_dumps(dc, raw=True), all_dumps(d4, raw=True)
    assert len(a) >= 4 and set(a) == set(b)
    for t in a:
        assert np.array_equal(np.sort(raw[t]["id"].ravel()), np.sort(raw4[t]["id"].ravel())), t
        assert matched(a[t]["points"], b[t]["points"]) < 1e-5, t
        if t > 0:
            check_mesh(raw[t], longest)
        if "doublings" in a[t]:
            assert np.array_equal(a[t]["doublings"], b[t]["doublings"]), t


@needs_apps
@pytest.mark.parametrize("app_name", ["partrac", "filaments"])
def test_reinjection_under_rk4cells_as_under_rk4(stamp_mesh, tmp_path, app_name):
    """u = (x, -y) on P1 triangles carries a strip (refined as it stretches)
    or pairs (halved back) out of the square; outside=reinject moves the
    stuck node, or the whole edge, back in at a random offset. Under RK4cells
    the same nodes get stuck, are moved by the same draws, and are not
    relocated: each takes the cell its position was located in."""
    if app_name == "partrac":
        binary = PARTRAC
        args = ("mode=triangle init_mode=strip_x La=2 x0=0 y0=0.3 z0=0 Nrw=21 Nrw_max=20000 "
                "refine=true refine_intv=0.1 ds_max=0.2 coarsen=false ds_min=1e-9")
        stuck = "nodes could not move"
    else:
        binary = FILAMENTS
        args = ("mode=triangle init_mode=pairs_xy_xy x0=0 y0=0 z0=0 Nrw=200 Nrw_max=200 ds_init=0.02 "
                "resize=doublings resize_target=ds_init resize_intv=0.05")
        stuck = "Some nodes are outside"
    out = {}
    for s in ("RK4cells", "RK4"):
        d = copy_case(stamp_mesh, tmp_path / s)
        (d / "timestamps.dat").write_text("0.0 out.h5\n3.0 out.h5\n")
        r = run_app(binary, d / "dolfin_params.dat", COMMON, args, "outside=reinject dt=0.05 T=1.0 dump_intv=0.5",
                    "scheme=" + s)
        out[s] = (d, r.stdout)
    (dc, sc), (d4, s4) = out["RK4cells"], out["RK4"]
    assert sc.count(stuck) >= 5 and sc.count(stuck) == s4.count(stuck)
    assert cells_counts(sc).get("relocations") == 0.0, sc
    a, b = all_dumps(dc), all_dumps(d4)
    assert set(a) == set(b)
    for t in a:
        assert len(a[t]["points"]) == len(b[t]["points"])
        assert np.abs(a[t]["points"] - b[t]["points"]).max() < 1e-5, t
        assert (np.abs(a[t]["points"][:, :2]) <= 2).all()


def next_id(d, drop=False):
    """The next node id saved in the checkpoint under d, None if it has none; with drop, deleted from it."""
    import h5py
    with h5py.File(checkpoint_file(d), "r+" if drop else "r") as h:
        n = int(h.attrs["next_id"]) if "next_id" in h.attrs else None
        if drop and n is not None:
            del h.attrs["next_id"]
    return n


def stopped_and_resumed(mesh_dir, tmp_path, case, scheme, drop=False):
    """(cont, split, checkpoint): a run to T = 0.5 and one stopped one step
    before the dump at 0.25 and resumed from its final checkpoint; with drop,
    the checkpoint loses its next id first, as a checkpoint without one."""
    binary = TET_CASES[case][0]
    end, stop = "T=0.5", "T=%r" % (0.25 - 0.00390625)
    cont, _ = tet_run(mesh_dir, tmp_path, "cont", case, "scheme=" + scheme, end)
    split, _ = tet_run(mesh_dir, tmp_path, "split", case, "scheme=" + scheme, stop)
    ck = read_checkpoint(split)
    ck["next_id"] = next_id(split, drop)
    run_app(binary, split / "dolfin_params.dat", COMMON, TET_CASES[case][1], end,
            "scheme=%s restart_folder=%s" % (scheme, checkpoint_folder(split)))
    return cont, split, ck


@needs_apps
@pytest.mark.parametrize("case,scheme", [("sheet", "RK4cells"), ("sheet", "RK4"), ("doublings", "RK4cells")])
def test_a_resumed_run_is_the_run_never_stopped(mesh_dir, tmp_path, case, scheme):
    """Stopped one step before a dump and resumed from the final checkpoint,
    which holds the cells and the next node id: every later dump bit for bit,
    ids too, with injection, refinement and coarsening (partrac) or resizing
    (filaments) on the way. In the sheet, coarsening has removed the newest
    node before the stop, so the ids present no longer give the next one."""
    cont, split, ck = stopped_and_resumed(mesh_dir, tmp_path, case, scheme)
    if scheme == "RK4cells":
        assert (ck["cell_id"] >= 0).all()
    assert ck["next_id"] > ck["id"].max()
    if case == "sheet":
        assert ck["next_id"] > ck["id"].max() + 1
    same_dumps(cont, split, after=0.25)


@needs_apps
def test_a_checkpoint_without_the_next_id_still_resumes(mesh_dir, tmp_path):
    """A checkpoint without the next id: the resumed run numbers new nodes on
    from the largest id present, and all else is the run never stopped's."""
    cont, split, ck = stopped_and_resumed(mesh_dir, tmp_path, "sheet", "RK4cells", drop=True)
    assert ck["next_id"] > ck["id"].max() + 1
    same_dumps(cont, split, after=0.25, skip=("id",))
    g = all_dumps(split, raw=True)
    t = min(t for t in g if t > 0.25)
    new = np.setdiff1d(g[t]["id"].ravel(), ck["id"].ravel())
    assert new.size and new.min() == ck["id"].max() + 1


@needs_apps
@pytest.mark.parametrize("case", ["sheet", "doublings"])
def test_one_and_eight_threads_agree_bit_for_bit(mesh_dir, tmp_path, case):
    """The steps are taken by the threads dynamically and the particle
    operations serially: every dump, the cells in the final checkpoint and
    the counts printed are the same at 1 and 8 threads, in partrac's sheet
    (injected, refined and coarsened) and filaments' pairs (resized)."""
    out = {}
    for n in (1, 8):
        d, stdout = tet_run(mesh_dir, tmp_path, "t%d" % n, case, "scheme=RK4cells num_threads=%d" % n,
                            "T=0.25")
        out[n] = (d, [l for l in stdout.splitlines() if "RK4cells:" in l])
    (a, sa), (b, sb) = out[1], out[8]
    assert sa and sa == sb
    same_dumps(a, b)
    assert "cell_id" in same_checkpoint(a, b)


@needs_apps
@pytest.mark.parametrize("app_name", ["partrac", "filaments"])
def test_an_analytic_field_takes_rk4s_loop(tmp_path, app_name):
    """On the unsteady ABC flow scheme=RK4cells is RK4: every dump, the
    statistics and the checkpoint bit for bit, with injection and refinement
    (partrac) or resizing (filaments)."""
    if app_name == "partrac":
        binary = PARTRAC
        args = ("init_mode=strip_x La=1 x0=3 y0=3 z0=3 Nrw=11 Nrw_max=20000 inject=true inject_intv=0.05 "
                "refine=true refine_intv=0.05 ds_max=0.1 coarsen=true coarsen_intv=0.05 ds_min=0.01")
    else:
        binary = FILAMENTS
        args = ("init_mode=pairs_xyz x0=3 y0=3 z0=3 Nrw=200 Nrw_max=200 ds_init=0.01 "
                "resize=doublings resize_target=ds_init resize_intv=0.05")
    runs = {}
    for s in ("RK4", "RK4cells"):
        params = copy_example(ABC, tmp_path / s)
        run_app(binary, params, COMMON, args, "mode=analytic dt=0.05 T=0.5 dump_intv=0.25 stat_intv=0.1 "
                "checkpoint_intv=0.25 scheme=" + s)
        runs[s] = params.parent
    same_dumps(runs["RK4"], runs["RK4cells"])
    same_checkpoint(runs["RK4"], runs["RK4cells"])
    same_stats(runs["RK4"], runs["RK4cells"])


@needs_apps
@pytest.mark.parametrize("app_name", ["partrac", "filaments"])
def test_rk4cells_refuses_diffusion_and_fields_without_cells(mesh_dir, felbm_dir, tmp_path, app_name):
    """As in the tracer apps: Dm > 0 is refused, and so is the constant lattice."""
    binary = PARTRAC if app_name == "partrac" else FILAMENTS
    d = copy_case(mesh_dir("tet"), tmp_path / "dm")
    r = run_app(binary, d / "dolfin_params.dat", args_for(app_name, "tet"), "scheme=RK4cells Dm=1e-4",
                check=False)
    assert r.returncode == 2 and "Dm must be 0" in r.stderr, r.stderr
    felbm = copy_case(felbm_dir, tmp_path / "felbm")
    with open(felbm / "felbm_params.dat", "a") as f:
        f.write("interpolation=constant\n")
    r = run_app(binary, felbm / "felbm_params.dat", args_for(app_name, "felbm"), "scheme=RK4cells", check=False)
    assert r.returncode == 2 and "scheme=RK4cells steps in the cells of a mesh" in r.stderr, r.stderr
    assert "use RK4" in r.stderr, r.stderr
