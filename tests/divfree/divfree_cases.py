"""The meshes, fields and cases the tests of python/divfree/divfree_clean.py build.

Importing this module imports the cleaner, which needs petsc4py and mpi4py:
without them it skips the importing test module, or fails it under
PARTRAC_REQUIRE_MPI."""

import functools
import glob
import itertools
import os
import re
import shutil
import subprocess
import sys

import numpy as np
import pytest

from paths import REPO, app

sys.path.insert(0, os.path.join(REPO, "python", "divfree"))

if not os.environ.get("PARTRAC_REQUIRE_MPI"):
    pytest.importorskip("petsc4py", reason="the cleaner solves with PETSc")
    pytest.importorskip("mpi4py", reason="the cleaner runs under MPI")

import divfree_clean as D          # noqa: E402

INTERPOL = app("interpol")


# ----------------------------------------------------------- meshes and fields


def channel2d(n=6, squash=None):
    """Unit square, periodic in x, walls at y = 0 and y = 1; squares split along
    the diagonal. With squash, the bottom row of cells is flattened by that
    factor, so those cells have that aspect ratio."""
    g = np.linspace(0, 1, n + 1)
    if squash is not None:
        g[1] = g[1] / squash
    X = np.array([(x, y) for x in g for y in g])
    nid = lambda i, j: i * (n + 1) + j
    cells = []
    for i in range(n):
        for j in range(n):
            a, b, c, d = nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)
            cells += [(a, b, c), (a, c, d)]
    return X, np.array(cells)


def channel3d(n=3, squash=None):
    """Unit cube, periodic in x and y, walls at z = 0 and z = 1; each box split
    into six Kuhn tets. With squash, the bottom layer is flattened."""
    g = np.linspace(0, 1, n + 1)
    if squash is not None:
        g[1] = g[1] / squash
    X = np.array([(x, y, z) for x in g for y in g for z in g])
    nid = lambda a: (a[0] * (n + 1) + a[1]) * (n + 1) + a[2]
    tets = []
    for base in itertools.product(range(n), repeat=3):
        for perm in itertools.permutations(range(3)):
            v = [np.array(base)]
            for a in perm:
                w = v[-1].copy()
                w[a] += 1
                v.append(w)
            tets.append([nid(w) for w in v])
    return X, np.array(tets)


def blocked2d(n=12):
    """channel2d with a slit wall across it at x = 0.5: the vertices there are
    duplicated, so the facets on the line are exterior on both sides and the
    throughput is zero for any divergence-free field."""
    X, cells = channel2d(n)
    on = np.nonzero(np.isclose(X[:, 0], 0.5))[0]
    right = X[cells].mean(axis=1)[:, 0] > 0.5
    X = np.vstack([X, X[on]])
    cells = cells.copy()
    for i, v in enumerate(on):
        sel = right & np.any(cells == v, axis=1)
        cells[sel] = np.where(cells[sel] == v, len(X) - len(on) + i, cells[sel])
    return X, cells


def blocked_field(p):
    """Divergence-free, at rest on the walls and on the whole line x = 0.5, so
    the nodes of the periodic seam at x = 0 are at rest and interior."""
    x, y = p[:, 0], p[:, 1]
    f, df = np.sin(2 * np.pi * x) ** 2, 2 * np.pi * np.sin(4 * np.pi * x)
    g, dg = (y * (1 - y)) ** 2, 2 * y * (1 - y) * (1 - 2 * y)
    return np.stack([f * dg, -df * g], axis=1)


def pinch(X, cells, rng, count=4, eps=0.02):
    """A few near-degenerate cells: interior vertices slid almost onto a
    neighbour, so the cells holding both are slivers."""
    X = X.copy()
    interior = np.nonzero(np.all((X > 1e-12) & (X < 1 - 1e-12), axis=1))[0]
    for v in rng.choice(interior, size=min(count, len(interior)), replace=False):
        w = cells[np.argmax(np.any(cells == v, axis=1))]
        other = [k for k in w if k != v][0]
        X[v] = X[other] + eps * (X[v] - X[other])
    return X


def poiseuille(X):
    """u = (y(1-y), 0[, 0]): divergence-free, quadratic, zero on the walls."""
    last = X[:, -1]
    U = np.zeros_like(X)
    U[:, 0] = last * (1 - last)
    return U


def smooth_noslip(p):
    """A smooth, non-polynomial, divergence-free field at rest on the walls of
    the channel: the curl of a stream function that vanishes there twice."""
    if p.shape[1] == 2:
        x, y = p[:, 0], p[:, 1]
        s = np.sin(2 * np.pi * x)
        c = 2 * np.pi * np.cos(2 * np.pi * x)
        f = (y * (1 - y)) ** 2
        df = 2 * y * (1 - y) * (1 - 2 * y)
        # psi = s f, u = (d psi/dy, -d psi/dx)
        return np.stack([s * df + y * (1 - y), -c * f], axis=1)
    x, y, z = p[:, 0], p[:, 1], p[:, 2]
    f = (z * (1 - z)) ** 2
    df = 2 * z * (1 - z) * (1 - 2 * z)
    a = np.sin(2 * np.pi * x) * f                    # vector potential (a, b, 0)
    b = np.sin(2 * np.pi * y) * f
    da = np.sin(2 * np.pi * x) * df
    db = np.sin(2 * np.pi * y) * df
    dbx = 2 * np.pi * np.cos(2 * np.pi * x) * f
    dby = 2 * np.pi * np.cos(2 * np.pi * y) * f
    # u = curl (a, b, 0), plus a plug flow that also vanishes on the walls
    return np.stack([-db + z * (1 - z), da, dby - dbx], axis=1)


def random_divfree(d, rng):
    """A divergence-free quadratic field: the curl of a random cubic."""
    if d == 2:
        c = rng.normal(size=10)

        def f(p):
            x, y = p[..., 0], p[..., 1]
            one = np.ones_like(x)
            dx = [3 * x ** 2, 0 * x, 2 * x * y, y ** 2, 2 * x, 0 * x, y, one, 0 * x, 0 * x]
            dy = [0 * x, 3 * y ** 2, x ** 2, 2 * x * y, 0 * x, 2 * y, x, 0 * x, one, 0 * x]
            return np.stack([sum(ci * a for ci, a in zip(c, dy)),
                             -sum(ci * a for ci, a in zip(c, dx))], -1)
        return f
    E = [e for e in itertools.product(range(4), repeat=3) if sum(e) <= 3]
    c = rng.normal(size=(3, len(E)))

    def dA(p, j):
        out = []
        for k in range(3):
            s = 0
            for i, e in enumerate(E):
                if e[j] == 0:
                    continue
                q = list(e)
                q[j] -= 1
                s = s + c[k][i] * e[j] * np.prod(p ** np.array(q), axis=-1)
            out.append(s)
        return np.stack(out, -1)

    def f(p):
        d0, d1, d2 = dA(p, 0), dA(p, 1), dA(p, 2)
        return np.stack([d1[..., 2] - d2[..., 1], d2[..., 0] - d0[..., 2],
                         d0[..., 1] - d1[..., 0]], -1)
    return f


# ------------------------------------------------------------------ topologies


def held_topo(X, cells, per, U):
    """A topology whose held set is the one this field gives: the nodes at rest
    that lie on an exterior facet the mesh does not pair."""
    rest = D.at_rest_nodes(np.linalg.norm(U, axis=1), float(np.abs(U).max()))
    return D.Topo(X, cells, per, at_rest=rest)


def walled(X, cells, dim, field=poiseuille):
    """(topology with the dataset's held set, the field as P2 nodes)."""
    per = [True] * (dim - 1) + [False]
    t = D.Topo(X, cells, per)
    U = field(t.node_x)
    return held_topo(X, cells, per, U), U


def facet_flux(t, U):
    """P2 flux out through every (cell, local facet), (ncells, nv)."""
    Um = U[t.nverts:]
    out = np.zeros((t.ncells, t.nv))
    for o in range(t.nv):
        f = [i for i in range(t.nv) if i != o]
        if t.dim == 2:
            v = (U[t.cells[:, f[0]]] + 4 * Um[t.cell_edge[:, o]] + U[t.cells[:, f[1]]]) / 6.0
        else:
            loc = [j for j, (a, b) in enumerate(D.LOC_EDGES[4]) if a != o and b != o]
            v = Um[t.cell_edge[:, loc]].sum(axis=1) / 3.0
        out[:, o] = np.einsum('ij,ij->i', v, t.facet_n[:, o])
    return out


def all_held_facets(t):
    """Exterior facets the mesh does not pair whose every P2 node is held."""
    out = np.zeros((t.ncells, t.nv), bool)
    for o in range(t.nv):
        loc = [j for j, (a, b) in enumerate(D.LOC_EDGES[t.nv]) if a != o and b != o]
        nodes = np.concatenate([t.facet_verts[:, o], t.nverts + t.cell_edge[:, loc]], axis=1)
        out[:, o] = t.held_node[nodes].all(axis=1)
    return out & t.facet_boundary


CASES = [(2, channel2d()), (3, channel3d())]
IDS = ["2d", "3d"]


def shuffled(X, cells, seed=5):
    """The same mesh with its vertices and cells numbered in random order and
    each cell's vertices rotated, so nothing may lean on a structured
    numbering."""
    rng = np.random.default_rng(seed)
    p = rng.permutation(len(X))
    new = np.empty_like(p)
    new[p] = np.arange(len(X))
    c = new[cells][rng.permutation(len(cells))]
    return X[p], np.roll(c, 1, axis=1)


def topo_cases():
    """The meshes the distributed tables are checked against the serial ones on,
    by name: (X, cells, periodic). Periodic in one to three directions, so the
    corners of the doubly and triply periodic ones reach their master through a
    chain of links; squashed, so a layer of slivers lies against a wall; and
    numbered at random."""
    return {
        "2d": channel2d() + ([True, False],),
        "3d": channel3d() + ([True, True, False],),
        "2d_xy": channel2d() + ([True, True],),
        "3d_x": channel3d(4) + ([True, False, False],),
        "3d_xyz": channel3d() + ([True, True, True],),
        "2d_squash": channel2d(6, 20) + ([True, False],),
        "3d_squash": channel3d(3, 20) + ([True, True, False],),
        "3d_xyz_shuffled": shuffled(*channel3d(4)) + ([True, True, True],),
    }

PATHS = [None, 0.0]
PATH_IDS = ["penalised", "smallest"]


# --------------------------------------------------------------- cases on disk


def write_case(folder, X, cells, fields, periodic, stamps=None, freq=None, degree=2,
               phi=None):
    """A dolfin HDF5 case written with the tool's own writer: the mesh, one file
    a stamp and the parameter file the loaders read. phi, one P2 scalar a stamp,
    adds a phase field."""
    folder.mkdir(parents=True, exist_ok=True)
    t = D.Topo(X, cells, periodic)
    D.write_mesh_h5(folder / "mesh.h5", X, cells)
    dim = X.shape[1]
    names = []
    for i, U in enumerate(fields):
        name = "up_%d.h5" % i
        names.append(name)
        vals = U if degree == 2 else U[:t.nverts]
        if i == 0:
            D.write_checkpoint(folder / name, "u", vals, t, degree)
        else:
            D.write_vector(folder / name, "u", vals, degree, dim, dim)
        p = np.zeros((t.nverts, 1))
        if i == 0:
            D.write_checkpoint(folder / name, "p", p, t, 1, mode="a")
        else:
            D.write_vector(folder / name, "p", p, 1, 1, dim, mode="a")
        if phi is not None:
            if i == 0:
                D.write_checkpoint(folder / name, "phi", phi[i], t, 2, mode="a")
            else:
                D.write_vector(folder / name, "phi", phi[i], 2, 1, dim, mode="a")
    per = "".join("periodic_%s=%s\n" % (a, "true" if p else "false")
                  for a, p in zip("xyz", list(periodic) + [False] * (3 - dim)))
    if freq is not None:
        (folder / "freqstamps.dat").write_text(
            "\n".join(" ".join(list(cols) + [n]) for cols, n in zip(freq, names)) + "\n")
        key = "freqstamps=freqstamps.dat\ntau=1.0\nt_min=0\nt_max=1e8\n"
    else:
        stamps = stamps or [str(i) for i in range(len(names))]
        (folder / "timestamps.dat").write_text(
            "\n".join("%s %s" % (s, n) for s, n in zip(stamps, names)) + "\n")
        key = "timestamps=timestamps.dat\n"
    (folder / "dolfin_params.dat").write_text(
        "velocity_space=P%d\npressure_space=P1\nmesh=mesh.h5\n" % degree + key + per + "rho=1.0\n"
        + ("include_phi=true\n" if phi is not None else ""))
    return folder / "dolfin_params.dat"


def relabel_by_global_id(folder, seed=3):
    """A case relabelled as a run on several ranks writes it: `mesh/cell_indices`
    is a permutation and not the identity, and every checkpoint's `cells` names
    the same global id for the same row, which is what dolfin writes on both
    sides. Nothing moves; only the labels of the rows change."""
    h5py = pytest.importorskip("h5py")
    with h5py.File(folder / "mesh.h5", "r+") as f:
        gid = np.random.default_rng(seed).permutation(len(f["mesh/topology"]))
        f["mesh/cell_indices"][...] = gid
    for name in sorted(os.listdir(folder)):
        if not name.endswith(".h5") or name == "mesh.h5":
            continue
        with h5py.File(folder / name, "r+") as f:
            for field in f:
                if "cells" in f[field]:
                    f[field + "/cells"][...] = gid
    return gid


def two_cells(folder):
    """channel2d(1): two cells, walls at rest, periodic in x, the diagonal's
    midpoint moved so that there is a step to take. Its two unknowns, the
    diagonal and the periodic pair of vertical edges, both lie in the first
    cell, so on two ranks the second cell's rank owns neither of them under
    the first-cell rule, and on three a rank has no cell at all."""
    X, cells = channel2d(1)
    t = D.Topo(X, cells, [True, False])
    U = smooth_noslip(t.node_x)
    mid = t.node_x[t.nverts:]
    U[t.nverts + np.nonzero(np.abs(mid[:, 0] - mid[:, 1]) < 1e-12)[0]] += [0.3, -0.1]
    return write_case(folder, X, cells, [U], [True, False])


def pocket_case(folder, dim, kind):
    """The 2D or 3D channel with a near-stagnant pocket across it: the field
    scaled by 1e-8 at 0.3 < x < 0.7, so the pocket's cells' own fluxes are far
    below FLUX_FLOOR of the stamp's largest and are measured against the floor.
    On three ranks under blocks the middle rank holds the pocket's cells and no
    others, so its own largest facet flux is 1e-8 of the stamp's.

    kind 'closed': walls at rest, so every boundary midpoint is held and the
    closed domain's flux rows leave the net flux untouched. In 2D the data's
    periodic seam carries a net flux of 1e-11 of the cells' scales summed,
    which is shared out over the cells and left in every one of them, below
    what a loader refuses. The 3D channel's rows are dependent beyond the
    constant (rank 159 of 162), so a net flux shared out by scale leaves the
    solve a part it cannot reach, and there the net is the data's round-off.
    kind 'moving wall': the top wall moves along x, so its midpoints are free,
    at the boundary weight. Two stamps, the second half the first."""
    X, cells = channel2d(6) if dim == 2 else channel3d(3)
    per = [True] * (dim - 1) + [False]
    t = D.Topo(X, cells, per)
    x = t.node_x
    U = smooth_noslip(x)
    if kind == "moving wall":
        U[:, 0] += x[:, -1] ** 2
    U[(x[:, 0] > 0.3) & (x[:, 0] < 0.7)] *= 1e-8
    if kind == "closed" and dim == 2:
        r, big = D.cell_flux(t, U)
        total = np.maximum(big, D.FLUX_FLOOR * big.max()).sum()
        seam = np.zeros_like(U)
        on = np.abs(x[:, 0] - 1.0) < 1e-12
        seam[on, 0] = U[on, 0]
        U = U + (1e-11 * total - r.sum()) / D.cell_flux(t, seam)[0].sum() * seam
    return write_case(folder, X, cells, [U, 0.5 * U], per, stamps=["0", "1"])


NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")
# what legitimately differs with the rank count: the rank count, times, the
# iteration path of the solve, the route and the output folder
RUN_ONLY = [(re.compile(r", \d+ ranks"), ""),
            (re.compile(r"\b\d+\.\d s\b"), "T s"),
            (re.compile(r"\d+ iterations, \d+ refinements"), "iterations, refinements"),
            (re.compile(r"(wrote \d+ stamps) to .*"), r"\1"),
            (re.compile(r"(refinement stopped): .*"), r"\1")]


def same_log(a, b, rtol=1e-10, atol=1e-11):
    """Two runs' printed reports, line for line: the same words, every integer
    the same and every other number within rtol of it or atol, what can differ
    with the rank count aside (RUN_ONLY). A number reduced over some ranks only
    moves by far more than round-off. atol is a decade below what the
    refinement asks of a cell's balance, which is what the solve's own
    round-off reaches; the fixtures' velocities and fluxes are of order one."""
    la = [l for l in a.split("\n") if l.strip()]
    lb = [l for l in b.split("\n") if l.strip()]
    assert len(la) == len(lb), (la, lb)
    for x, y in zip(la, lb):
        for pat, rep in RUN_ONLY:
            x, y = pat.sub(rep, x), pat.sub(rep, y)
        assert NUMBER.sub("#", x) == NUMBER.sub("#", y), (x, y)
        for p, q in zip(NUMBER.findall(x), NUMBER.findall(y)):
            if not any(c in p + q for c in ".eE"):
                assert p == q, (x, y)
                continue
            p, q = float(p), float(q)
            assert abs(p - q) <= rtol * max(abs(p), abs(q)) + atol, (x, y)


def read_back(out, names, field="u"):
    """The field of each named file of a cleaned case by node, read through the
    tool's own reader on the whole mesh: a file whose values sit at the wrong
    rows reads back as another field, or is refused, whatever wrote it."""
    case = D.read_case(out / "dolfin_params.dat")
    t = D.Topo(case["X"], case["cells"], case["periodic"])
    dof = D.DofTable(out / names[0], field, t, case["cell_indices"])
    return [dof.values(out / name) for name in names]


def same_output(a, b, scale, bound=1e-12):
    """Two cleaned cases' files against each other: the same datasets in every
    file, integers identical and floats within bound of scale, the same
    attributes, and the text files equal. The files are not compared as bytes,
    since HDF5 stores times in them."""
    h5py = pytest.importorskip("h5py")
    names = sorted(os.listdir(a))
    assert names == sorted(os.listdir(b))
    for name in names:
        if not name.endswith(".h5"):
            assert (a / name).read_text() == (b / name).read_text(), name
            continue
        with h5py.File(a / name, "r") as f, h5py.File(b / name, "r") as g:
            got = []
            f.visititems(lambda n, o: got.append(n) if isinstance(o, h5py.Dataset) else None)
            want = []
            g.visititems(lambda n, o: want.append(n) if isinstance(o, h5py.Dataset) else None)
            assert got == want, name
            for n in got:
                x, y = f[n][()], g[n][()]
                assert x.dtype == y.dtype and x.shape == y.shape, (name, n)
                if x.dtype.kind == "f":
                    assert np.abs(x - y).max(initial=0.0) <= bound * scale, (name, n)
                else:
                    assert np.array_equal(x, y), (name, n)
                assert sorted(f[n].attrs) == sorted(g[n].attrs), (name, n)


# ------------------------------------------------------------ the apps and MPI


@functools.lru_cache(maxsize=None)
def mpi_launcher():
    """The first launcher on the path that starts one job of two ranks for
    mpi4py, or None. A launcher from another MPI than mpi4py's starts
    independent one-rank processes instead, which each clean the whole case."""
    probe = "from mpi4py import MPI; print(MPI.COMM_WORLD.size)"
    for name in ("mpiexec", "mpirun", "mpiexec.mpich", "mpirun.mpich",
                 "mpiexec.openmpi", "mpirun.openmpi"):
        launcher = shutil.which(name)
        if launcher is None:
            continue
        try:
            r = subprocess.run([launcher, "-n", "2", sys.executable, "-c", probe],
                               capture_output=True, text=True, timeout=60,
                               env=dict(os.environ, OMP_NUM_THREADS="1"))
        except subprocess.TimeoutExpired:
            continue
        if r.returncode == 0 and r.stdout.split() == ["2", "2"]:
            return launcher
    return None


def launcher_or_skip():
    """mpi_launcher(), or the test skipped where there is none, and failed there
    with PARTRAC_REQUIRE_MPI set."""
    launcher = mpi_launcher()
    if launcher is None:
        if os.environ.get("PARTRAC_REQUIRE_MPI"):
            pytest.fail("no launcher here starts an MPI job for mpi4py")
        pytest.skip("no launcher here starts an MPI job for mpi4py")
    return launcher


def mpi_run(ranks, args, timeout=900, env=None):
    """python3 args as a job of that many ranks, one OpenMP thread a rank: the
    machine is shared out by rank here, and the solve is PETSc's, not
    OpenMP's. A job that has not ended by timeout fails the test: a rank that
    raised alone leaves the others waiting forever, and that is what it
    reports. The launcher is stopped with SIGTERM, which it passes on to the
    ranks; killed outright it would leave them running."""
    cmd = [launcher_or_skip(), "-n", str(ranks), sys.executable] + [str(a) for a in args]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                         env=dict(env or os.environ, OMP_NUM_THREADS="1"))
    try:
        out, err = p.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        p.terminate()
        try:
            out, err = p.communicate(timeout=30)
        except subprocess.TimeoutExpired:
            p.kill()
            out, err = p.communicate()
        pytest.fail("the job of %d ranks did not end in %d s:\n%s"
                    % (ranks, timeout, (out + err)[-3000:]))
    return subprocess.CompletedProcess(cmd, p.returncode, out, err)


def mpi_clean(cfg, out, ranks, extra=(), ok=True, timeout=900):
    """The tool as its own job on that many ranks. With ok=False the job is
    expected to fail and the caller reads the output."""
    tool = os.path.join(REPO, "python", "divfree", "divfree_clean.py")
    r = mpi_run(ranks, [tool, cfg, "--out", out] + list(extra), timeout)
    if ok:
        assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    return r


def interpol_probe(params, mode, npts, int_order=1):
    """What the interpol app gives at npts points of the case's box, a dict of
    its dumped arrays sorted by position. One thread: the app draws its points
    per thread, so two cases on one box are probed at the same points."""
    h5py = pytest.importorskip("h5py")
    r = subprocess.run([INTERPOL, str(params), "mode=" + mode, "Nrw=%d" % npts,
                        "int_order=%d" % int_order, "random=false", "seed=1", "t0=0"],
                       capture_output=True, text=True, timeout=600,
                       env=dict(os.environ, OMP_NUM_THREADS="1"))
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr
    [out] = glob.glob(os.path.join(os.path.dirname(str(params)), "Interpolation", "**",
                                   "interpolation.h5part"), recursive=True)
    with h5py.File(out, "r") as h:
        g = h[list(h.keys())[0]]
        d = {k: np.array(g[k]) for k in g}
    order = np.lexsort((d["z"], d["y"], d["x"]))
    return {k: v[order] for k, v in d.items()}
