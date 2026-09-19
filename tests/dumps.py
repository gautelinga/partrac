"""Reading partrac's HDF5 dumps.

partrac sorts its particles by cell as it runs, so a particle's slot moves
between dumps and between runs; only the `id` dataset names it. Anything that
compares dumps has to go through the ids, which these helpers do.
"""

import numpy as np


def by_id(group):
    """A dump group's datasets in id order, keyed by name (without `id`).

    Per-edge datasets (edges, dl, ...) are not permuted: sorting the particles
    keeps the edge order and only renames the nodes the edges refer to.
    """
    order = np.argsort(np.array(group["id"])[:, 0], kind="stable")
    return {k: np.array(group[k])[order] if group[k].shape[0] == len(order)
            else np.array(group[k]) for k in group if k != "id"}


def dump_at(case_dir, t):
    """Every dataset dumped at time t, in id order, from whichever file holds it."""
    import h5py
    for f in sorted(case_dir.rglob("data_from_t*.h5")):
        with h5py.File(f, "r") as h:
            for key in h:
                # group names are the dump time printed as text
                if abs(float(key) - t) < 1e-9:
                    return by_id(h[key])
    raise AssertionError("no dump at t = %g under %s" % (t, case_dir))
