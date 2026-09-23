"""Reading partrac's HDF5 dumps and statistics files.

partrac sorts its particles by cell as it runs, so a particle's slot moves
between dumps and between runs; only the `id` dataset names it. Anything that
compares dumps across runs that may sort differently has to go through the
ids, which by_id does; a comparison of runs that must agree slot for slot
reads the dump raw.
"""

import numpy as np
import pytest


def by_id(group):
    """A dump group's datasets in id order, keyed by name (without `id`).

    Per-edge datasets (edges, dl, ...) are not permuted: sorting the particles
    keeps the edge order and only renames the nodes the edges refer to. A
    dataset is taken as per-particle when it has one row per id, so on a
    closed curve, with as many edges as nodes, the per-edge datasets would be
    permuted too; no test here dumps one.
    """
    order = np.argsort(np.array(group["id"])[:, 0], kind="stable")
    return {k: np.array(group[k])[order] if group[k].shape[0] == len(order)
            else np.array(group[k]) for k in group if k != "id"}


def as_written(group):
    """A dump group's datasets as written, keyed by name, `id` included."""
    return {k: np.array(group[k]) for k in group}


def dump_at(case_dir, t, raw=False):
    """Every dataset dumped at time t, from whichever file holds it.

    In id order, or as written (with `id`) if raw.
    """
    h5py = pytest.importorskip("h5py")
    read = as_written if raw else by_id
    for f in sorted(case_dir.rglob("data_from_t*.h5")):
        with h5py.File(f, "r") as h:
            for key in h:
                # group names are the dump time printed as text
                if abs(float(key) - t) < 1e-9:
                    return read(h[key])
    raise AssertionError("no dump at t = %g under %s" % (t, case_dir))


def all_dumps(case_dir, pattern="data_from_t*.h5", raw=False):
    """time -> every dataset dumped then, over every file under case_dir matching pattern.

    In id order, or as written (with `id`) if raw. Times are rounded to 1e-9,
    so they can be looked up by the value they were dumped at.
    """
    h5py = pytest.importorskip("h5py")
    read = as_written if raw else by_id
    out = {}
    for f in sorted(case_dir.rglob(pattern)):
        with h5py.File(f, "r") as h:
            for key in h:
                out[round(float(key), 9)] = read(h[key])
    return out


def read_stats(path):
    """The statistics file at path, or the one under the folder path, as column -> array.

    Every row must have one value per column name: the header and the rows are
    written by separate functions, so a row that does not match its header is
    reported here rather than read into the wrong columns.
    """
    if path.is_dir():
        found = list(path.rglob("tdata_from_t*.dat"))
        assert len(found) == 1, found
        path = found[0]
    lines = [l for l in path.read_text().splitlines() if l.strip()]
    names = lines[0].lstrip("#").split()
    assert len(set(names)) == len(names), "%s: a column name repeats: %s" % (path.name, names)
    rows = [l.split() for l in lines[1:]]
    for row in rows:
        assert len(row) == len(names), (
            "%s: %d values under %d names: %s"
            % (path.name, len(row), len(names), names[len(row):] or row[len(names):]))
    data = np.array(rows, dtype=float).reshape(len(rows), len(names))
    return {name: data[:, j] for j, name in enumerate(names)}
