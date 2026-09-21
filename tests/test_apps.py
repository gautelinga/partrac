"""Every app starts, runs and writes something, on every input kind and scheme.

The apps share headers but each declares its own parameters, so a change to a
shared fragment can leave one unable to start while the rest are fine. The smoke
tests use the cheapest input each app accepts; they say nothing about whether
the physics is right, only that every app reaches every interpolator and every
scheme it offers. The analytic input is plane Poiseuille, |x| <= 1.

One table says what each app takes (APPS) and one what each input kind needs
(KINDS); the smoke runs are their product. The retired app names are wrappers
around the merged apps, read from the wrapper table in apps/CMakeLists.txt, and
must give the same bytes as the merged app with the pinned arguments.
args_for(name, kind) is the one place the arguments are put together, for these
tests and for the harness scripts.

Beyond the smoke tests: the initial refine and coarsen passes follow their own
flags; the statistics rows agree with the velocities dumped at the same time;
the tracer apps shuffle their particle slots; and the deterministic integrators
give bit-identical output on 1 and 4 threads.
"""

import os
import re
import shutil
import subprocess

import numpy as np
import pytest

from paths import REPO, app

EXAMPLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")

# dump_intv/dt must stay inside an int; 1e9/0.005 overflows the cast
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
    text = open(os.path.join(REPO, "apps", "CMakeLists.txt")).read()
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


def smoke_cells():
    """pytest params (app, kind, arguments) for every app x input kind x scheme it offers."""
    cells = []
    for name, (_, schemes) in APPS.items():
        for kind in KINDS:
            for scheme in (("explicit", "RK4") if schemes else (None,)):
                args = args_for(name, kind) + (" scheme=" + scheme if scheme else "")
                cells.append(pytest.param(name, kind, args,
                                          id="-".join(filter(None, (name, kind, scheme)))))
    return cells


# input files, which do not count as output of a run
KEEP = ("expr_params.dat", "dolfin_params.dat", "felbm_params.dat", "mesh.h5",
        "up_0.h5", "up_1.h5", "output_0.h5", "output_1.h5",
        "output_is_solid.h5", "timestamps.dat", "freqstamps.dat",
        "u.xdmf", "p.xdmf", "u.h5", "p.h5")


def case_dir(kind, tmp_path, mesh_dir, felbm_dir, xdmf_dir):
    """(folder, parameter file) of a private copy of the input, since output lands beside it.

    Every kind is read from the file's arrays, so none of them needs a build
    with dolfin; writing the mesh and XDMF inputs needs python dolfin, which
    their fixtures ask for.
    """
    d = tmp_path / "case"
    d.mkdir()
    if kind == "analytic":
        shutil.copy(EXAMPLE, d / "expr_params.dat")
        return d, d / "expr_params.dat"
    if kind in ("felbm", "xdmf"):
        src = felbm_dir if kind == "felbm" else xdmf_dir
        for f in os.listdir(src):
            shutil.copy(src / f, d / f)
        return d, d / ("felbm_params.dat" if kind == "felbm" else "dolfin_params.dat")
    src = mesh_dir(kind)
    for f in os.listdir(src):
        if f != "generate_up.py":
            shutil.copy(src / f, d / f)
    return d, d / "dolfin_params.dat"


@pytest.mark.parametrize("name,kind,args", smoke_cells())
def test_app_runs(name, kind, args, tmp_path, mesh_dir, felbm_dir, xdmf_dir):
    """The app exits 0 on this input and scheme, writes output, and its stats rows match the header.

    A failure means a user of that app cannot run on that kind of input, or
    gets a stats file whose columns are labelled wrongly.
    """
    binary = app(name)
    if not os.path.exists(binary):
        pytest.skip(name + " is not built")

    d, cfg = case_dir(kind, tmp_path, mesh_dir, felbm_dir, xdmf_dir)
    r = subprocess.run([binary, str(cfg)] + args.split(),
                       capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stdout + r.stderr
    assert [p for p in d.rglob("*") if p.is_file() and p.name not in KEEP], \
        "the app finished but wrote nothing"

    # the header and the rows are written by separate functions that must be
    # kept in step by hand, so every row must have one field per column name
    for stats in d.rglob("tdata_from_t*.dat"):
        rows = [l for l in stats.read_text().splitlines() if l.strip()]
        names = [h for h in rows[0].lstrip("# ").rstrip().split("\t") if h.strip()]
        for row in rows[1:]:
            fields = [v for v in row.rstrip().split("\t") if v.strip()]
            assert len(fields) == len(names), (
                "%s writes %d fields under %d names: %s"
                % (name, len(fields), len(names),
                   names[len(fields):] or fields[len(names):]))


@pytest.mark.skipif(not os.path.exists(app("partrac")), reason="apps are not built")
def test_every_built_app_is_listed():
    """Every executable in the bin directory is in APPS or the wrapper table.

    A new app must not arrive without smoke tests covering it.
    """
    # containment rather than equality: a dolfin-off build has only a subset.
    # Executables only: CMake also writes build_features.txt there.
    bindir = os.path.dirname(app("partrac"))
    built = {f for f in os.listdir(bindir)
             if os.path.isfile(os.path.join(bindir, f)) and os.access(os.path.join(bindir, f), os.X_OK)}
    listed = set(APPS) | set(WRAPPERS)
    assert built <= listed, built - listed


@pytest.mark.parametrize("name", sorted(WRAPPERS))
def test_a_retired_name_runs_the_merged_app(name, tmp_path, mesh_dir, felbm_dir, xdmf_dir):
    """A retired name gives byte-identical h5 output to the merged app with its arguments pinned.

    Scripts that still call a retired name must get the same results as the
    merged app, not a subtly different configuration.
    """
    if not os.path.exists(app(name)):
        pytest.skip(name + " is not built")
    target, parts = WRAPPERS[name]
    kind = wrapper_kind(name)
    args = args_for(name).split()
    out = {}
    for tag, binary, argv in (("wrapper", app(name), args),
                              ("merged", app(target), rewrite(args, parts))):
        parent = tmp_path / tag
        parent.mkdir()
        d, cfg = case_dir(kind, parent, mesh_dir, felbm_dir, xdmf_dir)
        r = subprocess.run([binary, str(cfg)] + argv,
                           capture_output=True, text=True, timeout=900)
        assert r.returncode == 0, r.stdout + r.stderr
        out[tag] = h5_datasets(d)
        assert out[tag], tag + " wrote no data"
    assert set(out["wrapper"]) == set(out["merged"])
    for k in out["wrapper"]:
        assert np.array_equal(out["wrapper"][k], out["merged"][k]), k


# --- and one thing that is not a smoke test ------------------------------------

SSS = ("Dm=0 dt=0.005 T=0.02 Nrw=100 Nrw_max=2000 dump_intv=1.0 stat_intv=1.0 "
       "mode=analytic init_mode=uniform_x int_order=1 dx_max=0.1 dxn=0.05 "
       "random=false seed=1 ds_max=0.4 ds_min=0.1")


@pytest.mark.skipif(not os.path.exists(app("static_space_stepper")),
                    reason="static_space_stepper is not built")
@pytest.mark.parametrize("refine,coarsen", [("true", "false"), ("false", "true")])
def test_the_initial_pass_follows_its_own_flag(tmp_path, refine, coarsen):
    """In static_space_stepper, the initial refinement runs iff refine=true and the initial coarsening iff coarsen=true.

    Each flag must control only its own pass, or a user who turns one off
    still gets its effect on the initial line.
    """
    d = tmp_path / "case"
    d.mkdir()
    shutil.copy(EXAMPLE, d / "expr_params.dat")
    r = subprocess.run([app("static_space_stepper"), str(d / "expr_params.dat")]
                       + SSS.split() + ["refine=" + refine, "coarsen=" + coarsen],
                       capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stdout + r.stderr
    assert ("Initial refinement" in r.stdout) == (refine == "true")
    assert ("Initial coarsening" in r.stdout) == (coarsen == "true")


@pytest.mark.skipif(not os.path.exists(app("static_space_stepper")),
                    reason="static_space_stepper is not built")
def test_the_stepper_statistics_read_the_fields_where_they_are_written(tmp_path):
    """Each static_space_stepper statistics row has the mean of the velocities dumped at the same step.

    Statistics and dumps share one refresh of the fields. Statistics written
    before the refresh would read unset velocities on the first row and stale
    ones on every row between dumps. ABC flow is used so the mean velocity
    changes along the march.
    """
    h5py = pytest.importorskip("h5py")
    d = tmp_path / "abc"
    d.mkdir()
    shutil.copy(os.path.join(REPO, "data_example", "abc_flow", "expr_params.dat"),
                d / "expr_params.dat")
    # stat_intv = dump_intv, so every statistics row has a dump to compare with
    args = ("mode=analytic init_mode=uniform_x Nrw=100 Nrw_max=2000 int_order=1 "
            "ds_max=0.4 ds_min=0.1 dt=0.05 dxn=0.05 dx_max=1 T=1e9 Ln=0.5 "
            "dump_intv=0.1 stat_intv=0.1 checkpoint_intv=1e9 random=false seed=1").split()
    r = subprocess.run([app("static_space_stepper"), str(d / "expr_params.dat")] + args,
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr

    lines = [l for l in next(d.rglob("tdata_from_t*.dat")).read_text().splitlines() if l.strip()]
    names = [h for h in lines[0].lstrip("# ").rstrip().split("\t") if h.strip()]
    rows = [dict(zip(names, map(float, l.split()))) for l in lines[1:]]
    dumped = {}
    for f in d.rglob("data_from_t*.h5"):
        with h5py.File(f, "r") as h:
            for g in h:
                dumped[float(g)] = np.array(h[g]["u"])
    assert len(rows) >= 4 and len(dumped) >= 4

    for row in rows:
        t = min(dumped, key=lambda s: abs(s - row["t"]))
        assert abs(t - row["t"]) < 1e-6, "no dump at the statistics row %g" % row["t"]
        mean = dumped[t].mean(axis=0)
        got = [row["ux_mean"], row["uy_mean"], row["uz_mean"]]
        # the statistics are printed to six significant digits
        assert np.allclose(got, mean, rtol=5e-6, atol=1e-9), \
            "row %g: u_mean %s, dumped %s" % (row["t"], got, mean)
    # the check can fail: the mean velocity changes along the march
    assert len({round(r["ux_mean"], 9) for r in rows}) > 1


@pytest.mark.skipif(not os.path.exists(app("filaments")),
                    reason="filaments is not built")
def test_the_filaments_statistics_read_the_fields_where_they_are_written(tmp_path):
    """Each filaments statistics row has the mean of the velocities dumped at the same time.

    filaments runs the shared loop, with one refresh of the fields for both
    the statistics and the dump; statistics written before it would read unset
    or stale velocities. Unsteady ABC flow makes the mean velocity change in time.
    """
    h5py = pytest.importorskip("h5py")
    d = tmp_path / "abc"
    d.mkdir()
    shutil.copy(os.path.join(REPO, "data_example", "abc_flow_unsteady", "expr_params.dat"),
                d / "expr_params.dat")
    # start the pairs at the centre of the example's [0, 2 pi]^3 domain
    pi = "3.14159265358979"
    args = ("mode=analytic init_mode=pairs_xyz Nrw=100 Nrw_max=2000 int_order=1 "
            "ds_max=0.4 ds_min=0.1 ds_init=0.1 x0=%s y0=%s z0=%s Dm=0 scheme=RK4 "
            "dt=0.01 T=0.5 dump_intv=0.1 stat_intv=0.1 checkpoint_intv=1e9 "
            "random=false seed=1" % (pi, pi, pi)).split()
    r = subprocess.run([app("filaments"), str(d / "expr_params.dat")] + args,
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr

    lines = [l for l in next(d.rglob("tdata_from_t*.dat")).read_text().splitlines() if l.strip()]
    names = [h for h in lines[0].lstrip("# ").rstrip().split("\t") if h.strip()]
    rows = [dict(zip(names, map(float, l.split()))) for l in lines[1:]]
    dumped = {}
    for f in d.rglob("data_from_t*.h5"):
        with h5py.File(f, "r") as h:
            for g in h:
                dumped[float(g)] = np.array(h[g]["u"])
    assert len(rows) >= 4 and len(dumped) >= 4

    for row in rows:
        t = min(dumped, key=lambda s: abs(s - row["t"]))
        assert abs(t - row["t"]) < 1e-6, "no dump at the statistics row %g" % row["t"]
        mean = dumped[t].mean(axis=0)
        got = [row["ux_mean"], row["uy_mean"], row["uz_mean"]]
        # the statistics are printed to six significant digits
        assert np.allclose(got, mean, rtol=5e-6, atol=1e-9), \
            "row %g: u_mean %s, dumped %s" % (row["t"], got, mean)
    # the check can fail: the mean velocity changes along the run
    assert len({round(r["ux_mean"], 9) for r in rows}) > 1


@pytest.mark.skipif(not os.path.exists(app("tracers")), reason="tracers is not built")
def test_the_tracer_apps_start_in_random_slot_order(tmp_path):
    """The tracer apps shuffle their particle slots once, before the first step.

    The core points initializer sorts its points by x, and under static OpenMP
    scheduling that would hand each thread one strip of the domain: on a mesh
    with uneven cost per particle, one thread gets the expensive strip and the
    others wait.
    """
    h5py = pytest.importorskip("h5py")
    d = tmp_path / "order"
    d.mkdir()
    shutil.copy(EXAMPLE, d / "expr_params.dat")
    # many particles and a single step, so the first dump shows the initial slot order
    argv = {}
    for a in (CORE + " mode=analytic init_mode=points_xy Nrw=4000 Nrw_max=4000 int_order=1 "
                     "dump_intv=0.005 T=0.005 random=false seed=1").split():
        argv[a.split("=")[0]] = a
    r = subprocess.run([app("tracers"), str(d / "expr_params.dat")] + list(argv.values()),
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr
    with h5py.File(next(d.rglob("data_from_t*.h5")), "r") as h:
        x = np.array(h[min(h.keys(), key=float)]["points"])[:, 0]
    span = x.max() - x.min()
    # each of 8 equal shares of the slots, as a static schedule would split
    # them, covers the domain rather than one strip of it
    for chunk in np.array_split(x, 8):
        assert chunk.max() - chunk.min() > 0.9 * span


# --- the parallel time loops ---------------------------------------------------

# These integrators draw no random numbers, so their result must not depend on
# how the work is divided between threads. Each loop has a parallel sibling in
# the same file that it is written to match; this checks that the match holds.
PARALLELISED = [
    # app, input kind, the parallel loop it exercises
    ("tracers", "triangle", "RK4Integrator::step<Point>"),
    ("filaments_triangleRK4", "triangle", "RK4Integrator::step<Point>"),
    ("tracertensors", "tet", "RK4Integrator::step<Tensor>"),
    ("tracervectors", "analytic", "RK4Integrator::step<Vector>"),
    ("static_space_stepper", "analytic", "SpatialIntegrator::step<Point>"),
    ("tracervectors_spatial", "triangle", "SpatialIntegrator::step<Vector>"),
    ("weighted_walkers", "analytic", "ExplicitIntegrator::step<Point>"),
]


def h5_datasets(d):
    """Every dataset in every h5 the run wrote (inputs excluded), keyed by file and path."""
    import h5py
    out = {}
    for f in sorted(d.rglob("*.h5")):
        if f.name in KEEP:
            continue
        h = h5py.File(f, "r")

        def walk(g, p=""):
            for k in g:
                if isinstance(g[k], h5py.Group):
                    walk(g[k], p + "/" + k)
                else:
                    out[f.name + p + "/" + k] = np.array(g[k])
        walk(h)
    return out


@pytest.mark.parametrize("name,kind,loop", PARALLELISED, ids=[a[0] for a in PARALLELISED])
def test_the_result_does_not_depend_on_the_thread_count(
        name, kind, loop, tmp_path, mesh_dir, felbm_dir, xdmf_dir):
    """A deterministic app writes bit-identical h5 output on 1 and on 4 threads.

    Without random numbers, any difference between thread counts is a race or
    an order dependence in the parallel loop, and results would change with
    the machine they run on.
    """
    if not os.path.exists(app(name)):
        pytest.skip(name + " is not built")
    args = args_for(name, kind)
    # enough particles and steps that every thread has work and several dumps
    # are compared; the spatial marchers stop at Ln instead of T
    work = ("Nrw=800 Nrw_max=8000 dt=0.005 dump_intv=0.05 stat_intv=1e9 "
            "checkpoint_intv=1e9 " +
            ("Ln=0.1 dxn=0.005 T=1e9" if name in ("static_space_stepper", "tracervectors_spatial")
             else "T=0.1"))

    out = {}
    for nthreads in (1, 4):
        parent = tmp_path / str(nthreads)
        parent.mkdir()
        d, cfg = case_dir(kind, parent, mesh_dir, felbm_dir, xdmf_dir)
        # later arguments override earlier ones with the same key
        argv = {}
        for a in (args + " " + work).split():
            argv[a.split("=")[0]] = a
        r = subprocess.run([app(name), str(cfg)] + list(argv.values()),
                           capture_output=True, text=True, timeout=900,
                           env=dict(os.environ, OMP_NUM_THREADS=str(nthreads)))
        assert r.returncode == 0, r.stdout + r.stderr
        out[nthreads] = h5_datasets(d)
        assert out[nthreads], name + " wrote no data to compare"

    assert set(out[1]) == set(out[4]), loop + " wrote different datasets"
    for k in out[1]:
        assert np.array_equal(out[1][k], out[4][k]), \
            "%s: %s differs between 1 and 4 threads" % (loop, k)
