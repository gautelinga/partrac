"""Running an app on a private copy of an example, and comparing runs.

Every app writes its output beside its parameter file, so a test runs on a
copy. The apps refuse a parameter given twice, and a test states a common set
of arguments and overrides a few of them, so the arguments are merged by key
here, the later value winning.
"""

import os
import pathlib
import re
import shutil
import subprocess

import numpy as np
import pytest

from dumps import all_dumps, deformation_gradient, read_stats


def merged(*argsets):
    """The arguments of every set, one per key, a later set's value winning.

    A set is a string of whitespace-separated key=value arguments, or a list
    of such strings and sets. A word with no = after a path continues that
    path (a restart folder may hold a space); anywhere else it stays an
    argument of its own, which the app then refuses.
    """
    argv = {}
    for args in argsets:
        words = []
        for w in args.split() if isinstance(args, str) else merged(*args):
            if words and "/" in words[-1] and "=" not in w and not w.startswith("--"):
                words[-1] += " " + w
            else:
                words.append(w)
        argv.update((w.split("=")[0], w) for w in words)
    return list(argv.values())


def copy_case(src, d):
    """Copy the input folder src, less the mesh generator script, to the new folder d; return d."""
    shutil.copytree(src, d, ignore=shutil.ignore_patterns("generate_up.py"))
    return d


def copy_example(example, d):
    """Copy the parameter file `example` into the folder d (made if missing); return the copy."""
    d.mkdir(parents=True, exist_ok=True)
    shutil.copy(example, d / os.path.basename(example))
    return d / os.path.basename(example)


def run_app(binary, params, *argsets, env=None, check=True, timeout=900):
    """Run binary on the parameter file params with the merged arguments; return the process.

    env adds to the environment. With check, the run must exit 0.
    """
    r = subprocess.run([str(binary), str(params)] + merged(*argsets),
                       capture_output=True, text=True, timeout=timeout,
                       env=dict(os.environ, **(env or {})))
    if check:
        assert r.returncode == 0, r.stdout + r.stderr
    return r


def checkpoint_file(d):
    """The one Checkpoints/checkpoint.h5 under d."""
    found = list(d.rglob("Checkpoints/checkpoint.h5"))
    assert len(found) == 1, found
    return found[0]


def checkpoint_folder(d):
    """The folder to resume from: the parent of the one Checkpoints folder under d."""
    return checkpoint_file(d).parent.parent


def read_checkpoint(d):
    """Every dataset of the one checkpoint under d, by name, as written."""
    h5py = pytest.importorskip("h5py")
    with h5py.File(checkpoint_file(d), "r") as h:
        return {k: h[k][()] for k in h}


def write_checkpoint(d, **datasets):
    """Replace datasets of the one checkpoint under d, each keeping its type and columns."""
    h5py = pytest.importorskip("h5py")
    with h5py.File(checkpoint_file(d), "r+") as h:
        for k, v in datasets.items():
            old = h[k]
            v = np.asarray(v, dtype=old.dtype).reshape(-1, old.shape[1])
            del h[k]
            h.create_dataset(k, data=v)


# One row a particle; points first
PER_PARTICLE = ("points", "id", "c", "t_loc", "rhohat", "w", "S", "Q", "logstretch", "U", "generation",
                "cell_id")


def put_points(d, points):
    """Put the particles of the checkpoint under d at points (n x 3 or n x dim).

    With as many points as particles the ids and fields are kept; otherwise
    the ids are 0..n-1, c runs from 0 to 1 and t_loc is 0, and a checkpoint
    carrying an element cannot be resized. The cells are unknown (-1): the
    resume locates the points.
    """
    points = np.asarray(points, dtype=float)
    points = np.c_[points, np.zeros((len(points), 3 - points.shape[1]))]
    ck = read_checkpoint(d)
    n = len(points)
    if n == len(ck["points"]):
        write_checkpoint(d, points=points, **({"cell_id": np.full(n, -1)} if "cell_id" in ck else {}))
        return
    new = {"points": points, "id": np.arange(n), "c": np.linspace(0., 1., n), "t_loc": np.zeros(n),
           "cell_id": np.full(n, -1)}
    extra = [k for k in PER_PARTICLE if k in ck and k not in new]
    assert not extra, "cannot resize %s" % extra
    write_checkpoint(d, **{k: v for k, v in new.items() if k in ck})


def write_text_checkpoint(folder, F_whole=False):
    """Rewrite the checkpoint.h5 in the Checkpoints folder `folder` as the text
    files older runs wrote, and remove it. With F_whole a tensor element is
    written whole, as F.ten, as before its factors were carried."""
    h5py = pytest.importorskip("h5py")
    folder = pathlib.Path(folder)
    with h5py.File(folder / "checkpoint.h5", "r") as h:
        ck = {k: h[k][()] for k in h}

    def rows(name, *columns):
        with open(folder / name, "w") as f:
            for r in zip(*columns):
                f.write(" ".join(repr(int(v)) if isinstance(v, np.integer) else "%.17g" % v
                                 for v in r) + "\n")

    def edges(name, suffix):
        e = ck["edges" + suffix]
        rows(name, e[:, 0], e[:, 1], ck["dl0" + suffix][:, 0], ck["edge_tau" + suffix][:, 0],
             ck["edge_rho_prev" + suffix][:, 0])

    def columns(a):
        return [a[:, j] for j in range(a.shape[1])]

    rows("positions.pos", *columns(ck["points"]))
    rows("id.list", ck["id"][:, 0])
    rows("colors.col", ck["c"][:, 0])
    edges("edges.edge", "")
    f = ck["face_edges"]
    rows("faces.face", f[:, 0], f[:, 1], f[:, 2], ck["dA0"][:, 0], ck["face_tau"][:, 0],
         ck["face_rho_prev"][:, 0])
    for name, stem in (("doublings", "doublings.list"), ("edges_inlet", "edges_inlet.list"),
                       ("nodes_inlet", "nodes_inlet.list"), ("t_loc", "t_loc.dat"),
                       ("w", "w.dat"), ("S", "S.dat"), ("generation", "generation.dat")):
        if name in ck:
            rows(stem, ck[name][:, 0])
    if "positions_inj" in ck:
        rows("positions_inj.pos", *columns(ck["positions_inj"]))
        edges("edges_inj.edge", "_inj")
    for name, stem in (("rhohat", "rhohat.vec"), ("logstretch", "logstretch.vec"), ("U", "U.vec")):
        if name in ck and not (F_whole and name != "rhohat"):
            rows(stem, *columns(ck[name]))
    if "Q" in ck:
        if F_whole:
            F = deformation_gradient({"Q": ck["Q"], "logstretch": ck["logstretch"], "U": ck["U"]})
            rows("F.ten", *columns(F.reshape(-1, 9)))
        else:
            rows("Q.ten", *columns(ck["Q"]))
    (folder / "checkpoint.h5").unlink()


def continuous_and_resumed(binary, example, root, args, stop, end, env=None, text=False):
    """(cont, split) folders: a run to `end`, and one to `stop` resumed from its checkpoint to `end`.

    args, stop and end are argument sets; stop and end override args, and end
    serves both the uninterrupted run and the resumed one. The final
    checkpoint is written one step past the stop, and the resumed run writes
    into the tree of the run it resumes. With text, the checkpoint is
    rewritten in the old text format before the resume.
    """
    cont = copy_example(example, root / "cont")
    split = copy_example(example, root / "split")
    run_app(binary, cont, args, end, env=env)
    run_app(binary, split, args, stop, env=env)
    resume = checkpoint_folder(split.parent)
    if text:
        write_text_checkpoint(resume / "Checkpoints")
    run_app(binary, split, args, end, ["restart_folder=%s" % resume], env=env)
    return cont.parent, split.parent


# --- comparing runs ------------------------------------------------------------

def cells_counts(stdout):
    """RK4cells' counts per particle-step, by name, from the first report in stdout; {} without one."""
    m = re.search(r"RK4cells: (\d+) particle-steps; per particle-step (.*)", stdout)
    if not m:
        return {}
    return {k: float(v) for v, k in re.findall(r"([0-9.e+-]+) ([a-z ]+?)(?:,|$)", m.group(2))}


def same_dumps(a, b, after=-1.0, skip=()):
    """The dumps of the runs in the folders a and b after a time: the same
    times and datasets, every dataset but those in skip bit for bit, as
    written; the times compared."""
    da, db = all_dumps(a, raw=True), all_dumps(b, raw=True)
    times = sorted(t for t in da if t > after)
    assert times and times == sorted(t for t in db if t > after), (sorted(da), sorted(db))
    for t in times:
        assert set(da[t]) == set(db[t]), t
        for k in set(da[t]) - set(skip):
            assert np.array_equal(da[t][k], db[t][k]), (t, k)
    return times


def same_checkpoint(a, b):
    """The checkpoints under the folders a and b: the same datasets, each bit
    for bit; a's checkpoint."""
    ca, cb = read_checkpoint(a), read_checkpoint(b)
    assert set(ca) == set(cb), (sorted(ca), sorted(cb))
    for k in ca:
        assert np.array_equal(ca[k], cb[k]), k
    return ca


def same_stats(a, b):
    """The statistics files under the folders a and b: the same columns, each bit for bit."""
    sa, sb = read_stats(a), read_stats(b)
    assert list(sa) == list(sb), (list(sa), list(sb))
    for k in sa:
        assert np.array_equal(sa[k], sb[k]), k


def halving_ratios(runs):
    """For runs [(x, F)] at dt halving from one to the next, the ratios of the
    successive changes of F and of x, p90 over the particles: [(rF, rx)], one
    for each three runs in a row."""
    ratios = []
    for (x0, F0), (x1, F1), (x2, F2) in zip(runs, runs[1:], runs[2:]):
        dF = [np.percentile(np.linalg.norm(a - b, axis=(1, 2)), 90) for a, b in ((F0, F1), (F1, F2))]
        dx = [np.percentile(np.linalg.norm(a - b, axis=1), 90) for a, b in ((x0, x1), (x1, x2))]
        ratios.append((dF[0] / dF[1], dx[0] / dx[1]))
    return ratios
