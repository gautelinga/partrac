"""The smoke-run table shared by the tests and the harness scripts, the
synthetic FELBM writer, and the cache that builds each generated input once
per test session.

One table says what each app takes (APPS) and one what each input kind needs
(KINDS); a smoke run is one of their product. The retired app names are
wrappers around the merged apps, read from the wrapper table in
apps/CMakeLists.txt. args_for(name, kind) is the one place the arguments are
put together.
"""

import fcntl
import os
import re
import shutil

from paths import REPO

EXAMPLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")

CORE = "Dm=0 dt=0.005 T=0.02 Nrw=100 Nrw_max=2000 dump_intv=1.0 stat_intv=1.0"

# input kind -> (the mode it is read with, a point well inside its domain)
KINDS = {
    "analytic":     ("analytic", "x0=0 y0=0 z0=0"),            # plane Poiseuille, |x| <= 1
    "triangle":     ("triangle", "x0=0.5 y0=0.5 z0=0"),
    "tet":          ("tet", "x0=0.5 y0=0.5 z0=0.5"),
    "trianglefreq": ("trianglefreq", "x0=0.5 y0=0.5 z0=0"),
    "tetfreq":      ("tetfreq", "x0=0.5 y0=0.5 z0=0.5"),
    "felbm":        ("felbm", "x0=8 y0=8 z0=8"),               # 16^3, solid at the z ends
    "xdmf":         ("xdmftriangle", "x0=0.5 y0=0.5 z0=0"),
}

# app -> (what it takes beyond CORE, a mode and a centre; whether it takes scheme=)
APPS = {
    "partrac":               ("init_mode=uniform_x int_order=1 ds_max=0.4 ds_min=0.1", True),
    "filaments":             ("init_mode=pairs_xy int_order=1 ds_max=0.4 ds_min=0.1 ds_init=0.05", True),
    "tracers":               ("init_mode=points_xy int_order=1", True),
    "tracervectors":         ("init_mode=points_xy", True),
    "tracertensors":         ("init_mode=points_xy", True),
    "static_space_stepper":  ("init_mode=uniform_x int_order=1 ds_max=0.4 ds_min=0.1 dx_max=0.1 dxn=0.05", False),
    "tracervectors_spatial": ("init_mode=points_xy dxn=0.005", False),
    "weighted_walkers":      ("init_mode=strip_y_x int_order=1 ds_max=2.0 La=0.5 Lb=0.0 Ln=1e9 Lt=0 refine_intv=0.05", False),
    # a probe of the fields, not a run: no time, no seed, no initial state
    "interpol":              ("Nrw=100 int_order=1", False),
}


def wrapper_table():
    """(retired name, target app, {ARGS, RENAME, COPY, DROP} lists) for each wrapper in apps/CMakeLists.txt."""
    # read from the one table that defines the wrappers, so the tests cannot
    # drift from what they check
    with open(os.path.join(REPO, "apps", "CMakeLists.txt")) as f:
        text = f.read()
    table = []
    for m in re.finditer(r"partrac_add_wrapper\((\w+)\s+TARGET\s+(\w+)\s+([^)]*)\)", text):
        parts = {"ARGS": [], "RENAME": [], "COPY": [], "DROP": []}
        kind = None
        for tok in m.group(3).split():
            if tok in parts:
                kind = tok
            else:
                parts[kind].append(tok)
        table.append((m.group(1), m.group(2), parts))
    return table


WRAPPERS = {name: (target, parts) for name, target, parts in wrapper_table()}


def rewrite(args, parts):
    """The caller's argument list as the wrapper hands it on to the merged app."""
    rename = dict(a.split("=") for a in parts["RENAME"])
    copy = dict(a.split("=") for a in parts["COPY"])
    out = list(parts["ARGS"])
    for a in args:
        key, val = a.split("=", 1)
        if key in copy:
            out.append(copy[key] + "=" + val)
        if key not in parts["DROP"]:
            out.append(rename.get(key, key) + "=" + val)
    return out


def wrapper_kind(name):
    """The input kind a retired name reads: the one its pinned mode names."""
    pinned = dict(a.split("=", 1) for a in WRAPPERS[name][1]["ARGS"])
    return next(kind for kind, (mode, _) in KINDS.items() if mode == pinned["mode"])


def args_for(name, kind=None):
    """The argument string of a smoke run of an app or retired name on an input kind.

    A retired name defaults to its own kind. It is given what its app is given,
    less what the wrapper pins or derives from another argument, and under the
    parameter names the retired app used.
    """
    if name in WRAPPERS:
        target, parts = WRAPPERS[name]
        kind = kind or wrapper_kind(name)
        pinned = {a.split("=")[0] for a in parts["ARGS"]}
        made = {a.split("=")[1] for a in parts["COPY"]}
        old_name = {new: old for old, new in (a.split("=") for a in parts["RENAME"])}
        out = []
        for a in args_for(target, kind).split():
            key, val = a.split("=", 1)
            if key in pinned or key in made:
                continue
            out.append(old_name.get(key, key) + "=" + val)
        return " ".join(out)
    own, _ = APPS[name]
    mode, centre = KINDS[kind or "analytic"]
    if name == "interpol":
        return "mode=%s %s" % (mode, own)
    return "%s mode=%s %s %s random=false seed=1" % (CORE, mode, centre, own)


def write_felbm(d, fields, solid, times=(0, 100), extra=""):
    """Write a FELBM case in d: one output_k.h5 per time holding fields[k], and the solid mask.

    fields[k] maps u_x, u_y, u_z, density and pressure to (nx, ny, nz) arrays,
    indexed [ix, iy, iz]; they are written transposed, as the solver writes
    them, so the loader reads x fastest. solid is written as given, first h5
    axis z. felbm_params.dat names the timestamps file and the mask, followed
    by the lines in extra.
    """
    import numpy as np
    import h5py
    with h5py.File(d / "output_is_solid.h5", "w") as f:
        f.create_dataset("is_solid", data=solid)
    for k, stamp in enumerate(fields):
        with h5py.File(d / ("output_%d.h5" % k), "w") as f:
            for name, a in stamp.items():
                f.create_dataset(name, data=np.transpose(a, (2, 1, 0)).astype(float))
    (d / "timestamps.dat").write_text("".join("%g\toutput_%d.h5\n" % (t, k)
                                              for k, t in enumerate(times)))
    (d / "felbm_params.dat").write_text(
        "timestamps=timestamps.dat\nis_solid_file=output_is_solid.h5\n" + extra)


# --- generated inputs, built once per session ----------------------------------

def session_root(tmp_path_factory):
    """The folder every worker of this pytest session shares.

    Under pytest-xdist each worker has its own basetemp inside the session's;
    without xdist the session's basetemp is the worker's. Nothing above it is
    per session, so a cache there could serve an input built by an older
    generator.
    """
    base = tmp_path_factory.getbasetemp()
    return base.parent if "PYTEST_XDIST_WORKER" in os.environ else base


def shared_dir(tmp_path_factory, name, build):
    """A folder `name` that build(folder) fills once per session, whichever worker asks first.

    The lock and the completion mark sit beside the folder, not in it, so a
    test that copies the folder copies only what build wrote. A build that
    failed leaves no mark, and the next caller builds again.
    """
    root = session_root(tmp_path_factory)
    d = root / ("shared_" + name)
    done = root / ("shared_" + name + ".done")
    with open(root / ("shared_" + name + ".lock"), "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)           # released when the file closes
        if not done.exists():
            if d.exists():
                shutil.rmtree(d)
            d.mkdir()
            build(d)
            done.touch()
    return d
