"""scheme=RK4cells in the tracer apps: RK4 held in the cells of a mesh field.

On a mesh field the velocity gradient jumps at every facet, so an RK4 step
across one loses the scheme's order: F converges at first order and positions
at second. RK4cells cuts each major step dt where the path is predicted to
meet a facet, lands on it by a Newton correction and steps on in the cell
beyond, so every RK4 step it takes sees one polynomial and the order is
restored. dt is its only control. A field without cells (analytic) has no
kinks in space, and the scheme hands it to RK4's loop.

The kernel's own properties (a step meeting no facet is RK4's bit for bit; a
particle on a facet or beside a vertex) are in test_cells.cpp. Here: the
order on a P2 tet field, the linear flow against its exact flow map (F and a
line element), the analytic field, restarts with the cells in the
checkpoint, the thread count, and diffusion, which it refuses.
"""

import os

import numpy as np
import pytest

from cases import ABC, KINDS, LINEAR, stamp_case
from dumps import deformation_gradient, dump_at
from paths import app
from runs import (cells_counts, checkpoint_file, checkpoint_folder, copy_case, copy_example, halving_ratios,
                  read_checkpoint, run_app, same_checkpoint, same_dumps, same_stats)
from test_stamp_cuts import TENSOR_BASE, alternating, flow_map

TRACERS = app("tracers")
VECTORS = app("tracervectors")
TENSORS = app("tracertensors")

needs_apps = pytest.mark.skipif(not all(os.path.exists(a) for a in (TRACERS, VECTORS, TENSORS)),
                                reason="the tracer apps are not built")

# 200 particles over the periodic 10^3 P2 tet case, u = (sin 2 pi y, sin 2 pi z, sin 2 pi x)
TET = ("mode=tet init_mode=points_xyz x0=0.5 y0=0.5 z0=0.5 Nrw=200 Nrw_max=200 Dm=0 int_order=1 "
       "stat_intv=1e9 checkpoint_intv=1e9 random=false seed=1")


def tet_case(mesh_dir, root, name):
    """A copy of the P2 tet case under root/name; its parameter file."""
    return copy_case(mesh_dir("tet"), root / name) / "dolfin_params.dat"


def tet_run(mesh_dir, root, name, *args):
    """tracertensors on a copy of the tet case; (case folder, stdout)."""
    params = tet_case(mesh_dir, root, name)
    r = run_app(TENSORS, params, TET, *args)
    return params.parent, r.stdout


def x_and_F(d, t):
    """Positions and F at t, in id order."""
    g = dump_at(d, t)
    return g["points"], deformation_gradient(g)


@needs_apps
def test_rk4cells_converges_at_fourth_order_where_rk4_does_not(mesh_dir, tmp_path):
    """Halving dt divides the change of F and x by about 16 under RK4cells
    (p90 over particles, T = 0.3), from dt |u| near h/4 down; under RK4 F's
    change halves (first order: a step across a facet is wrong to O(dt^2))."""
    T = 0.3
    out = {}
    for scheme in ("RK4cells", "RK4"):
        runs = [tet_run(mesh_dir, tmp_path, "%s%d" % (scheme, k), "scheme=" + scheme,
                        "dt=%r T=%r dump_intv=%r" % (0.025 / 2 ** k, T, T * (1 + 1e-9)))[0] for k in range(4)]
        out[scheme] = halving_ratios([x_and_F(d, T) for d in runs])
    for rF, rx in out["RK4cells"]:
        assert rF > 9 and rx > 9, out
    assert out["RK4cells"][-1][0] > 12 and out["RK4cells"][-1][1] > 12, out
    assert all(rF < 3 for rF, _ in out["RK4"]), out


@needs_apps
@pytest.mark.parametrize("kind", ["trianglefreq", "tetfreq"])
def test_frequency_fields_converge_at_fourth_order(mesh_dir, tmp_path, kind):
    """The frequency loaders' cells are crossed as the stamped ones are, no
    step falling back, and halving dt divides the change of F and x by about
    16 (p90 over particles, T = 0.4)."""
    runs = []
    for dt in (0.04, 0.02, 0.01):
        d = copy_case(mesh_dir(kind), tmp_path / ("dt%g" % dt))
        r = run_app(TENSORS, d / "dolfin_params.dat", "mode=%s init_mode=points_xyz %s" % (kind, KINDS[kind][1]),
                    "Nrw=100 Nrw_max=100 Dm=0 int_order=1 stat_intv=1e9 checkpoint_intv=1e9 random=false seed=1",
                    "scheme=RK4cells dt=%r T=0.4 dump_intv=0.4" % dt)
        c = cells_counts(r.stdout)
        assert c and c["crossings"] > 0.1 and c["fallbacks"] == 0, r.stdout
        runs.append(x_and_F(d, 0.4))
    [(rF, rx)] = halving_ratios(runs)
    assert rF > 12 and rx > 12, (rF, rx)


def linear_errors(stamp_mesh, root, scheme, dt, stamps, exact, T=1.2):
    """(max |F - Phi|, max |x - Phi x0|) at T of tracertensors on the linear flow, against the flow map exact."""
    params = stamp_case(stamp_mesh, root, stamps, "%s%g" % (scheme, dt))
    r = run_app(TENSORS, params, TENSOR_BASE, "scheme=%s dt=%r T=%r dump_intv=%r" % (scheme, dt, T, T + dt / 2))
    g0, gT = dump_at(params.parent, 0.0), dump_at(params.parent, T)
    # far from the walls, where no particle leaves the square
    keep = np.abs(g0["points"][:, 0]) < 0.8
    assert keep.sum() > 5
    F = deformation_gradient(gT)[keep]
    eF = np.abs(F[:, :2, :2] - exact).max()
    ex = np.abs(gT["points"][keep, :2] - g0["points"][keep, :2] @ exact.T).max()
    return eF, ex, r.stdout


@needs_apps
@pytest.mark.parametrize("field", ["steady", "alternating"])
def test_the_linear_flow_at_fourth_order_and_no_less_accurate_than_rk4(stamp_mesh, tmp_path, field):
    """u = A x on P1 triangles is exact on the mesh, J = A in every cell, so
    RK4 keeps its order there and the facets cost RK4cells nothing in
    accuracy: F and x converge at fourth order against the flow map (the
    matrix exponential for a steady A; fine RK4 between the stamps for A(t)
    alternating at stamps off the dt grid), and RK4cells, whose steps a
    crossing only shortens, is as accurate as RK4 at every dt (to 0.1%)."""
    from scipy.linalg import expm
    T = 1.2
    if field == "steady":
        stamps = [(0.0, "lin_a"), (3.0, "lin_a")]
        exact = expm(np.array(LINEAR["lin_a"]) * T)
    else:
        stamps = alternating("lin")
        exact = flow_map(stamps, T)
    errs = {s: [linear_errors(stamp_mesh, tmp_path, s, dt, stamps, exact) for dt in (0.1, 0.05, 0.025)]
            for s in ("RK4cells", "RK4")}
    cells = errs["RK4cells"]
    assert cells_counts(cells[0][2]).get("crossings", 0) > 0, cells[0][2]
    for (eF0, ex0, _), (eF1, ex1, _) in zip(cells[:-1], cells[1:]):
        assert eF0 / eF1 > 12 and ex0 / ex1 > 12, cells
    for (eF, ex, _), (eF4, ex4, _) in zip(cells, errs["RK4"]):
        assert eF <= eF4 * (1 + 1e-3) and ex <= ex4 * (1 + 1e-3), (cells, errs["RK4"])


@needs_apps
def test_a_line_element_landed_on_facets_follows_the_flow_map(stamp_mesh, tmp_path):
    """tracervectors on the steady u = A x on P1 triangles: a line element is
    carried to each facet a step lands on, so its direction and its log
    stretch w converge at fourth order to the flow map's, exp(A T) n0."""
    from scipy.linalg import expm
    T = 1.2
    Phi = np.zeros((3, 3))
    Phi[:2, :2] = np.array(LINEAR["lin_a"]) * T
    Phi = expm(Phi)
    errs = []
    for dt in (0.05, 0.025, 0.0125):
        params = stamp_case(stamp_mesh, tmp_path, [(0.0, "lin_a"), (3.0, "lin_a")], "v%g" % dt)
        r = run_app(VECTORS, params, TENSOR_BASE, "scheme=RK4cells dt=%r T=%r dump_intv=%r" % (dt, T, T + dt / 2))
        assert cells_counts(r.stdout).get("crossings", 0) > 0, r.stdout
        g0, gT = dump_at(params.parent, 0.0), dump_at(params.parent, T)
        keep = np.abs(g0["points"][:, 0]) < 0.8
        v = (Phi @ g0["n"][keep].T).T
        en = np.abs(gT["n"][keep] - v / np.linalg.norm(v, axis=1)[:, None]).max()
        ew = np.abs(gT["w"][keep].ravel() - np.log(np.linalg.norm(v, axis=1))).max()
        errs.append((en, ew))
    for (n0, w0), (n1, w1) in zip(errs[:-1], errs[1:]):
        assert n0 / n1 > 12 and w0 / w1 > 12, errs


@needs_apps
@pytest.mark.parametrize("app_name", ["tracers", "tracervectors", "tracertensors"])
def test_an_analytic_field_takes_rk4s_loop(tmp_path, app_name):
    """On the unsteady ABC flow scheme=RK4cells is RK4: every dump, the
    statistics and the checkpoint bit for bit, so a run switches schemes
    without switching fields."""
    args = ("mode=analytic init_mode=points_xyz x0=3 y0=3 z0=3 Nrw=100 Nrw_max=100 Dm=0 int_order=1 "
            "dt=0.05 T=0.5 dump_intv=0.25 stat_intv=0.1 checkpoint_intv=0.25 random=false seed=1")
    runs = {}
    for scheme in ("RK4", "RK4cells"):
        params = copy_example(ABC, tmp_path / scheme)
        run_app(app(app_name), params, args, "scheme=" + scheme)
        runs[scheme] = params.parent
    assert len(same_dumps(runs["RK4"], runs["RK4cells"])) == 3
    same_checkpoint(runs["RK4"], runs["RK4cells"])
    same_stats(runs["RK4"], runs["RK4cells"])


@needs_apps
def test_a_resumed_run_continues_in_the_cell_entered(mesh_dir, tmp_path):
    """The checkpoint holds each particle's cell. A particle that landed on a
    facet sits on both of its cells; resumed in the one it entered, it steps
    on as the run never stopped would, bit for bit. A checkpoint without the
    cells is located afresh and resumes too."""
    common = ["scheme=RK4cells", "dt=0.0625", "dump_intv=0.25"]
    end = ["T=1.0", "checkpoint_intv=1e9"]
    cont, _ = tet_run(mesh_dir, tmp_path, "cont", common, end)
    split, _ = tet_run(mesh_dir, tmp_path, "split", common, ["T=0.4375", "checkpoint_intv=0.4375"])
    ck = read_checkpoint(split)
    assert ck["cell_id"].shape == (200, 1) and (ck["cell_id"] >= 0).all()
    resume = checkpoint_folder(split)
    run_app(TENSORS, split / "dolfin_params.dat", TET, common, end, ["restart_folder=%s" % resume])
    assert len(same_dumps(cont, split, after=0.5)) == 2
    # Without the cells
    old, _ = tet_run(mesh_dir, tmp_path, "old", common, ["T=0.4375", "checkpoint_intv=0.4375"])
    h5py = pytest.importorskip("h5py")
    with h5py.File(checkpoint_file(old), "r+") as h:
        del h["cell_id"]
    run_app(TENSORS, old / "dolfin_params.dat", TET, common, end, ["restart_folder=%s" % checkpoint_folder(old)])
    xa, Fa = x_and_F(cont, 1.0)
    xo, Fo = x_and_F(old, 1.0)
    assert np.abs(xa - xo).max() < 1e-9 and np.abs(Fa - Fo).max() < 1e-7


@needs_apps
def test_one_and_eight_threads_agree_bit_for_bit(mesh_dir, tmp_path):
    """Particles are taken by the threads dynamically, but each one's steps
    depend on nothing but its own state: every dump, the cells in the
    checkpoint and the counts printed are the same at 1 and 8 threads."""
    out = {}
    for n in (1, 8):
        d, stdout = tet_run(mesh_dir, tmp_path, "t%d" % n, "scheme=RK4cells", "dt=0.05", "T=0.5",
                            "dump_intv=0.25", "checkpoint_intv=0.25", "num_threads=%d" % n)
        out[n] = (d, [l for l in stdout.splitlines() if "RK4cells:" in l])
    (a, sa), (b, sb) = out[1], out[8]
    assert sa and sa == sb
    same_dumps(a, b)
    assert "cell_id" in same_checkpoint(a, b)


@needs_apps
def test_rk4cells_refuses_diffusion(mesh_dir, tmp_path):
    """RK4cells is deterministic: Dm > 0 is refused. The fields without cells
    it refuses are test_cells_lattice's."""
    r = run_app(TENSORS, tet_case(mesh_dir, tmp_path, "dm"), TET, "scheme=RK4cells dt=0.05 T=0.1 Dm=1e-4",
                check=False)
    assert r.returncode == 2 and "Dm must be 0" in r.stderr, r.stderr
