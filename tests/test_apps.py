"""Every app starts, runs and steps, on every input kind and scheme.

The apps share headers but each declares its own parameters, so a change to a
shared fragment can leave one unable to start while the rest are fine. The smoke
tests use the cheapest input each app accepts; they say nothing about whether
the physics is right, only that every app reaches every interpolator and every
scheme it offers. The analytic input is plane Poiseuille, |x| <= 1.

The smoke runs are the product of the tables in cases.py: what each app takes
(APPS) and what each input kind needs (KINDS). The retired app names are
wrappers around the merged apps and must give the same bytes as the merged app
with the pinned arguments.

Beyond the smoke tests: the initial refine and coarsen passes follow their own
flags; the statistics rows agree with the velocities dumped at the same time
and do not depend on what the dump holds; the tracer apps shuffle their
particle slots; the deterministic integrators give bit-identical output on 1
and 4 threads, and partrac's statistics agree to their printed precision; a
diffusive run is reproducible for a given seed and thread count; and a bare
parameter file name reads the case in the working directory.
"""

import os
import re
import subprocess

import numpy as np
import pytest

from cases import (APPS, CORE, EXAMPLE, KINDS, WRAPPERS, args_for, rewrite,
                   wrapper_kind)
from dumps import all_dumps, read_stats
from paths import REPO, app
from runs import checkpoint_folder, copy_case, copy_example, run_app


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
    if kind == "analytic":
        return d, copy_example(EXAMPLE, d)
    copy_case({"felbm": felbm_dir, "xdmf": xdmf_dir}.get(kind) or mesh_dir(kind), d)
    return d, d / ("felbm_params.dat" if kind == "felbm" else "dolfin_params.dat")


# how far each app is taken so that a dump after its first step exists: the
# time-stepping apps dump at T, the marchers go by path length. T is also the
# marchers' cutoff in local time, which the slow felbm flow reaches within a
# step of tracervectors_spatial.
FIRST_STEP = {
    "static_space_stepper": "Ln=0.1 dump_intv=0.05",
    "tracervectors_spatial": "Ln=0.01 dump_intv=0.01 T=1e9",
}


@pytest.mark.parametrize("name,kind,args", smoke_cells())
def test_app_runs(name, kind, args, tmp_path, mesh_dir, felbm_dir, xdmf_dir):
    """The app exits 0 on this input and scheme, dumps a state past its first step, and its stats rows match the header.

    A failure means a user of that app cannot run on that kind of input, or
    gets a stats file whose columns are labelled wrongly. Exiting 0 is not
    enough: an app that stops before its first step still writes its initial
    dump and parameters.
    """
    binary = app(name)
    if not os.path.exists(binary):
        pytest.skip(name + " is not built")

    d, cfg = case_dir(kind, tmp_path, mesh_dir, felbm_dir, xdmf_dir)
    if name == "interpol":
        run_app(binary, cfg, args)
        assert list(d.rglob("interpolation.h5part")), "interpol finished but wrote nothing"
        return
    run_app(binary, cfg, args, FIRST_STEP.get(name, "dump_intv=0.02"))
    assert max(all_dumps(d, raw=True), default=0) > 0, "no dump after the first step"
    # read_stats checks every row against the header
    for stats in d.rglob("tdata_from_t*.dat"):
        read_stats(stats)


def unconditional_apps():
    """The apps apps/CMakeLists.txt adds outside any if block, and the retired
    names whose target is one of them (a wrapper is skipped without its target)."""
    with open(os.path.join(REPO, "apps", "CMakeLists.txt")) as f:
        text = f.read()
    names, depth = set(), 0
    for m in re.finditer(r"^\s*(if|endif|partrac_add_app|partrac_add_wrapper)\s*\(\s*(\w*)", text, re.M):
        if m.group(1) == "if":
            depth += 1
        elif m.group(1) == "endif":
            depth -= 1
        elif depth == 0:
            names.add(m.group(2))
    return {n for n in names if n not in WRAPPERS or WRAPPERS[n][0] in names}


@pytest.mark.skipif(not os.path.exists(app("partrac")), reason="apps are not built")
def test_every_built_app_is_listed():
    """Every executable in the bin directory is in APPS or the wrapper table,
    and every app built whatever the options are is there.

    A new app must not arrive without smoke tests covering it, and an app
    that stopped being built would turn its tests into skips.
    """
    # Executables only: CMake also writes build_features.txt there. The
    # BUILD_APP_ options can leave out the apps added inside an if block.
    bindir = os.path.dirname(app("partrac"))
    built = {f for f in os.listdir(bindir)
             if os.path.isfile(os.path.join(bindir, f)) and os.access(os.path.join(bindir, f), os.X_OK)}
    listed = set(APPS) | set(WRAPPERS)
    assert built <= listed, built - listed
    always = unconditional_apps()
    assert "partrac" in always and always <= listed, always - listed
    assert always <= built, always - built


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
        run_app(binary, cfg, argv)
        out[tag] = h5_datasets(d)
        assert out[tag], tag + " wrote no data"
    assert set(out["wrapper"]) == set(out["merged"])
    for k in out["wrapper"]:
        assert np.array_equal(out["wrapper"][k], out["merged"][k]), k


# --- and a few things that are not smoke tests ---------------------------------

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
    r = run_app(app("static_space_stepper"), copy_example(EXAMPLE, tmp_path / "case"),
                SSS, ["refine=" + refine, "coarsen=" + coarsen])
    assert ("Initial refinement" in r.stdout) == (refine == "true")
    assert ("Initial coarsening" in r.stdout) == (coarsen == "true")


def assert_stats_match_dumps(d):
    """Every statistics row under d has the mean of the velocities dumped at its time, and the mean changes."""
    stats = read_stats(d)
    dumped = {t: g["u"] for t, g in all_dumps(d, raw=True).items()}
    assert len(stats["t"]) >= 4 and len(dumped) >= 4
    for i, t_row in enumerate(stats["t"]):
        t = min(dumped, key=lambda s: abs(s - t_row))
        assert abs(t - t_row) < 1e-6, "no dump at the statistics row %g" % t_row
        mean = dumped[t].mean(axis=0)
        got = [stats[c][i] for c in ("ux_mean", "uy_mean", "uz_mean")]
        # the statistics are printed to six significant digits
        assert np.allclose(got, mean, rtol=5e-6, atol=1e-9), \
            "row %g: u_mean %s, dumped %s" % (t_row, got, mean)
    # the check can fail: the mean velocity changes along the run
    assert len({round(v, 9) for v in stats["ux_mean"]}) > 1


@pytest.mark.skipif(not os.path.exists(app("static_space_stepper")),
                    reason="static_space_stepper is not built")
def test_the_stepper_statistics_read_the_fields_where_they_are_written(tmp_path):
    """Each static_space_stepper statistics row has the mean of the velocities dumped at the same step.

    Statistics and dumps share one refresh of the fields. Statistics written
    before the refresh would read unset velocities on the first row and stale
    ones on every row between dumps. ABC flow is used so the mean velocity
    changes along the march.
    """
    pytest.importorskip("h5py")
    d = tmp_path / "abc"
    # stat_intv = dump_intv, so every statistics row has a dump to compare with
    run_app(app("static_space_stepper"),
            copy_example(os.path.join(REPO, "data_example", "abc_flow", "expr_params.dat"), d),
            "mode=analytic init_mode=uniform_x Nrw=100 Nrw_max=2000 int_order=1 "
            "ds_max=0.4 ds_min=0.1 dt=0.05 dxn=0.05 dx_max=1 T=1e9 Ln=0.5 "
            "dump_intv=0.1 stat_intv=0.1 checkpoint_intv=1e9 random=false seed=1")
    assert_stats_match_dumps(d)


@pytest.mark.skipif(not os.path.exists(app("filaments")),
                    reason="filaments is not built")
def test_the_filaments_statistics_read_the_fields_where_they_are_written(tmp_path):
    """Each filaments statistics row has the mean of the velocities dumped at the same time.

    filaments runs the shared loop, with one refresh of the fields for both
    the statistics and the dump; statistics written before it would read unset
    or stale velocities. Unsteady ABC flow makes the mean velocity change in time.
    """
    pytest.importorskip("h5py")
    d = tmp_path / "abc"
    # start the pairs at the centre of the example's [0, 2 pi]^3 domain
    pi = "3.14159265358979"
    run_app(app("filaments"),
            copy_example(os.path.join(REPO, "data_example", "abc_flow_unsteady", "expr_params.dat"), d),
            "mode=analytic init_mode=pairs_xyz Nrw=100 Nrw_max=2000 int_order=1 "
            "ds_max=0.4 ds_min=0.1 ds_init=0.1 x0=%s y0=%s z0=%s Dm=0 scheme=RK4 "
            "dt=0.01 T=0.5 dump_intv=0.1 stat_intv=0.1 checkpoint_intv=1e9 "
            "random=false seed=1" % (pi, pi, pi))
    assert_stats_match_dumps(d)


# Dm = 0 removes the Brownian term and random=false pins the seed, so two runs
# with the same arguments must agree bit for bit
PARTRAC_BASE = ("mode=analytic init_mode=uniform_x Nrw=100 Nrw_max=5000 ds_max=0.4 "
                "ds_min=0.1 Dm=0 int_order=1 dt=0.01 dump_intv=0.05 stat_intv=0.05 "
                "checkpoint_intv=0.05 random=false seed=3")


@pytest.mark.skipif(not os.path.exists(app("partrac")), reason="partrac is not built")
def test_statistics_do_not_depend_on_dump_verbosity(tmp_path):
    """minimal_output selects only what goes into the h5 dump. The statistics
    read the velocity whether or not it is dumped, so the mean velocities are
    identical either way, while the minimal dump still holds fewer datasets.
    Otherwise trimming the dump to save disk would change the statistics."""
    means, sizes = {}, {}
    for minimal in ("true", "false"):
        d = tmp_path / ("mo_" + minimal)
        run_app(app("partrac"), copy_example(EXAMPLE, d), PARTRAC_BASE,
                ["T=0.05", "minimal_output=" + minimal])
        st = read_stats(d)
        means[minimal] = [st[c][0] for c in ("ux_mean", "uy_mean", "uz_mean")]
        assert np.isfinite(means[minimal]).all(), means[minimal]
        sizes[minimal] = len(all_dumps(d, raw=True)[0.0])

    assert means["true"] == means["false"]
    assert abs(means["false"][2]) > 1e-6  # the flow is axial, so uz_mean is nonzero and the equality is not 0 == 0
    assert sizes["true"] < sizes["false"]


@pytest.mark.skipif(not os.path.exists(app("tracers")), reason="tracers is not built")
def test_the_tracer_apps_start_in_random_slot_order(tmp_path):
    """The tracer apps shuffle their particle slots once, before the first step.

    The core points initializer sorts its points by x, and under static OpenMP
    scheduling that would hand each thread one strip of the domain: on a mesh
    with uneven cost per particle, one thread gets the expensive strip and the
    others wait.
    """
    pytest.importorskip("h5py")
    d = tmp_path / "order"
    # many particles and a single step, so the first dump shows the initial slot order
    run_app(app("tracers"), copy_example(EXAMPLE, d), CORE,
            "mode=analytic init_mode=points_xy Nrw=4000 Nrw_max=4000 int_order=1 "
            "dump_intv=0.005 T=0.005 random=false seed=1")
    x = all_dumps(d, raw=True)[0.0]["points"][:, 0]
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
    # app, input kind, the parallel loop it exercises, arguments beyond the smoke run's
    ("tracers", "triangle", "RK4Integrator::step<Point>", ""),
    ("filaments", "triangle", "RK4Integrator::step<Point>",
     "scheme=RK4 resize=doublings resize_target=ds_init outside=reinject"),
    ("tracertensors", "tet", "RK4Integrator::step<Tensor>", ""),
    ("tracervectors", "analytic", "RK4Integrator::step<Vector>", ""),
    ("static_space_stepper", "analytic", "SpatialIntegrator::step<Point>", ""),
    ("tracervectors_spatial", "triangle", "SpatialIntegrator::step<Vector>", ""),
    ("weighted_walkers", "analytic", "ExplicitIntegrator::step<Point>", ""),
]


def h5_datasets(d):
    """Every dataset in every h5 the run wrote (inputs excluded), keyed by file and path."""
    import h5py
    out = {}

    def walk(g, name, p=""):
        for k in g:
            if isinstance(g[k], h5py.Group):
                walk(g[k], name, p + "/" + k)
            else:
                out[name + p + "/" + k] = np.array(g[k])

    for f in sorted(d.rglob("*.h5")):
        if f.name not in KEEP:
            with h5py.File(f, "r") as h:
                walk(h, f.name)
    return out


@pytest.mark.parametrize("name,kind,loop,extra", PARALLELISED, ids=[a[0] for a in PARALLELISED])
def test_the_result_does_not_depend_on_the_thread_count(
        name, kind, loop, extra, tmp_path, mesh_dir, felbm_dir, xdmf_dir):
    """A deterministic app writes bit-identical h5 output on 1 and on 4 threads.

    Without random numbers, any difference between thread counts is a race or
    an order dependence in the parallel loop, and results would change with
    the machine they run on.
    """
    if not os.path.exists(app(name)):
        pytest.skip(name + " is not built")
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
        run_app(app(name), cfg, args_for(name, kind), extra, work,
                env={"OMP_NUM_THREADS": str(nthreads)})
        out[nthreads] = h5_datasets(d)
        assert out[nthreads], name + " wrote no data to compare"

    assert set(out[1]) == set(out[4]), loop + " wrote different datasets"
    for k in out[1]:
        assert np.array_equal(out[1][k], out[4][k]), \
            "%s: %s differs between 1 and 4 threads" % (loop, k)


def final_positions(d):
    """The positions in the checkpoint the run under d ended on."""
    return np.loadtxt(checkpoint_folder(d) / "Checkpoints" / "positions.pos")


# uniform_x is a strip, so it exercises only the strip half of the statistics;
# the sheet half is a separate set of reductions and needs its own run
SHEET_STATS = ("mode=analytic init_mode=sheet_xy La=0.6 Lb=0.6 ds_init=0.05 "
               "x0=0 y0=0 z0=0 Nrw=100 Nrw_max=200000 refine=true ds_max=0.1 "
               "coarsen=false ds_min=1e-9 Dm=0 int_order=2 dt=0.01 "
               "dump_intv=1e9 checkpoint_intv=0.05 random=false seed=3")


@pytest.mark.skipif(not os.path.exists(app("partrac")), reason="partrac is not built")
@pytest.mark.parametrize("kind", ["strip", "sheet"])
def test_partrac_positions_and_statistics_do_not_depend_on_the_thread_count(tmp_path, kind):
    """Positions on 4 threads are bit-identical to one thread, for a strip and
    a sheet, since every parallel loop writes one entry per particle. The
    statistics are parallel reductions that add per-thread partials in
    whatever order the threads finish, so they may move in the last bits and
    are compared to their printed precision; a larger difference is a broken
    reduction."""
    if kind == "strip":
        args = [PARTRAC_BASE, "T=0.2 int_order=2 stat_intv=0.01"]
    else:
        args = [SHEET_STATS, "T=0.1 stat_intv=0.01"]
    ref, other = tmp_path / "t1", tmp_path / "t4"
    run_app(app("partrac"), copy_example(EXAMPLE, ref), *args, "num_threads=1")
    run_app(app("partrac"), copy_example(EXAMPLE, other), *args, "num_threads=4")

    assert np.array_equal(final_positions(ref), final_positions(other))
    a, b = read_stats(ref), read_stats(other)
    assert list(a) == list(b), "the two runs wrote different columns"
    assert len(a["t"]) == len(b["t"]), "the two runs wrote a different number of rows"
    # the file is written with six significant digits, so a value on a rounding
    # boundary can flip its last printed digit: rtol=1e-5. The atol covers the
    # columns that are zero by symmetry, such as x_mean, which hold only
    # cancellation noise.
    bad = [c for c in a if not np.allclose(a[c], b[c], rtol=1e-5, atol=1e-12)]
    assert not bad, "columns differ beyond a reduction's rounding: %s" % bad


@pytest.mark.skipif(not os.path.exists(app("partrac")), reason="partrac is not built")
def test_a_diffusive_run_is_reproducible(tmp_path):
    """With Brownian noise on 4 threads, two runs with the same seed are
    bit-identical. Each thread draws from its own generator, so this holds only
    while the schedule is static and particle i always lands on the same
    thread; the Dm = 0 tests cannot detect a change of schedule."""
    args = [PARTRAC_BASE, "Dm=1e-4 T=0.1 num_threads=4"]
    a, b = tmp_path / "d1", tmp_path / "d2"
    run_app(app("partrac"), copy_example(EXAMPLE, a), *args)
    run_app(app("partrac"), copy_example(EXAMPLE, b), *args)
    assert np.array_equal(final_positions(a), final_positions(b))


@pytest.mark.skipif(not os.path.exists(app("interpol")), reason="interpol is not built")
def test_a_bare_parameter_file_name_reads_the_case_in_the_working_directory(tmp_path, mesh_dir):
    """The loaders find the mesh and the stamps beside the parameter file, cut at
    its last slash; a bare name has none, so the working directory is taken."""
    d = tmp_path / "case"
    copy_case(mesh_dir("tet"), d)
    r = subprocess.run([app("interpol"), "dolfin_params.dat", "mode=tet", "Nrw=10", "int_order=1",
                        "random=false", "seed=1", "t0=0"],
                       cwd=d, capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
