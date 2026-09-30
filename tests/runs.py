"""Running an app on a private copy of an example.

Every app writes its output beside its parameter file, so a test runs on a
copy. The apps refuse a parameter given twice, and a test states a common set
of arguments and overrides a few of them, so the arguments are merged by key
here, the later value winning.
"""

import os
import pathlib
import shutil
import subprocess

import numpy as np
import pytest


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
PER_PARTICLE = ("points", "id", "c", "t_loc", "rhohat", "w", "S", "Q", "logstretch", "U", "generation")


def put_points(d, points):
    """Put the particles of the checkpoint under d at points (n x 3 or n x dim).

    With as many points as particles the ids and fields are kept; otherwise
    the ids are 0..n-1, c runs from 0 to 1 and t_loc is 0, and a checkpoint
    carrying an element cannot be resized.
    """
    points = np.asarray(points, dtype=float)
    points = np.c_[points, np.zeros((len(points), 3 - points.shape[1]))]
    ck = read_checkpoint(d)
    n = len(points)
    if n == len(ck["points"]):
        write_checkpoint(d, points=points)
        return
    new = {"points": points, "id": np.arange(n), "c": np.linspace(0., 1., n), "t_loc": np.zeros(n)}
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
            from dumps import deformation_gradient
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
