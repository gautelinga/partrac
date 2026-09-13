"""Every app starts, runs and writes something.

The apps share headers but each declares its own parameters, so a change to a
shared fragment can leave one unable to start while the rest are fine. These are
smoke tests on the cheapest input each app accepts; they say nothing about
whether the physics is right, only that the app is reachable.

The parameter lists are spelled out per app rather than shared, because the
schemas genuinely differ: some take int_order, some do not.
"""

import os
import shutil
import subprocess

import numpy as np
import pytest

from paths import REPO, app

EXAMPLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")

# dump_intv/dt must stay inside an int; 1e9/0.005 overflows the cast
CORE = "Dm=0 dt=0.005 T=0.02 Nrw=100 Nrw_max=2000 dump_intv=1.0 stat_intv=1.0"

# name, input kind, arguments
APPS = [
    ("partrac", "analytic",
     CORE + " mode=analytic init_mode=uniform_x int_order=1 ds_max=0.4 ds_min=0.1 random=false seed=1"),
    ("filaments", "analytic",
     CORE + " mode=analytic init_mode=pairs_xyz int_order=1 ds_max=0.4 ds_min=0.1 ds_init=0.1 random=false seed=1"),
    ("interpol", "analytic",
     "mode=analytic Nrw=100 int_order=1"),
    ("static_space_stepper", "analytic",
     CORE + " mode=analytic init_mode=uniform_x int_order=1 ds_max=0.4 ds_min=0.1 dx_max=0.1 dxn=0.05 random=false seed=1"),
    ("tracervectors_analyticRK4", "analytic",
     CORE + " init_mode=points_xy random=false seed=1"),
    ("weighted_walkers", "analytic",
     CORE + " init_mode=strip_y_x int_order=1 ds_max=2.0 La=0.5 Lb=0.0 Ln=1e9 Lt=0 refine_intv=0.05 random=false seed=1"),
    ("omp_test", "triangle",
     CORE + " init_mode=points_xy int_order=1 seed=1"),
    ("tracers_triangleRK4", "triangle",
     CORE + " init_mode=points_xy int_order=1 random=false seed=1"),
    ("tracervectors_triangle_spatial", "triangle",
     CORE + " init_mode=points_xy ds_max=0.4 random=false seed=1"),
    # always uses RandomPairsInitializer, which retries a rejected pair
    # unchanged, so the starting point has to be inside the mesh
    ("filaments_triangleRK4", "triangle",
     CORE + " init_mode=pairs_xy int_order=1 ds_init=0.05 x0=0.5 y0=0.5 random=false seed=1"),
    ("tracertensors_triangleRK4", "tet",
     CORE + " init_mode=points_xy random=false seed=1"),
    ("tracervectors_trianglefreqRK4", "trianglefreq",
     CORE + " init_mode=points_xy random=false seed=1"),
    # always uses RandomPairsInitializer, as filaments_triangleRK4 does
    ("filaments_felbmRK4", "felbm",
     CORE + " init_mode=pairs_xy int_order=1 ds_init=0.5 x0=8 y0=8 z0=8 random=false seed=1"),
    ("tracervectors_triangleRK4", "xdmf",
     CORE + " init_mode=points_xy random=false seed=1"),
]

KEEP = ("expr_params.dat", "dolfin_params.dat", "felbm_params.dat", "mesh.h5",
        "up_0.h5", "up_1.h5", "output_0.h5", "output_1.h5",
        "output_is_solid.h5", "timestamps.dat", "freqstamps.dat",
        "u.xdmf", "p.xdmf", "u.h5", "p.h5")


def case_dir(kind, tmp_path, mesh_dir, felbm_dir, xdmf_dir):
    """A private copy of the input, since output lands beside it."""
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


@pytest.mark.parametrize("name,kind,args", APPS, ids=[a[0] for a in APPS])
def test_app_runs(name, kind, args, tmp_path, mesh_dir, felbm_dir, xdmf_dir):
    if kind.startswith("no_data:"):
        pytest.skip("needs " + kind.split(":", 1)[1] + ", which is not in the repository")
    binary = app(name)
    if not os.path.exists(binary):
        pytest.skip(name + " is not built")

    d, cfg = case_dir(kind, tmp_path, mesh_dir, felbm_dir, xdmf_dir)
    r = subprocess.run([binary, str(cfg)] + args.split(),
                       capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stdout + r.stderr
    assert [p for p in d.rglob("*") if p.is_file() and p.name not in KEEP], \
        "the app finished but wrote nothing"

    # the header and the row are written by two functions kept in step by hand,
    # in three separate copies of this code; they have drifted five times
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
    # a new app should not arrive without an entry here; a dolfin-off build
    # has only a subset, so this is containment rather than equality
    built = set(os.listdir(os.path.dirname(app("partrac"))))
    listed = {a[0] for a in APPS}
    assert built <= listed, built - listed


# --- and one thing that is not a smoke test ------------------------------------

SSS = ("Dm=0 dt=0.005 T=0.02 Nrw=100 Nrw_max=2000 dump_intv=1.0 stat_intv=1.0 "
       "mode=analytic init_mode=uniform_x int_order=1 dx_max=0.1 dxn=0.05 "
       "random=false seed=1 ds_max=0.4 ds_min=0.1")


@pytest.mark.skipif(not os.path.exists(app("static_space_stepper")),
                    reason="static_space_stepper is not built")
@pytest.mark.parametrize("refine,coarsen", [("true", "false"), ("false", "true")])
def test_the_initial_pass_follows_its_own_flag(tmp_path, refine, coarsen):
    # the two used to share one block gated on refine, so coarsen=true
    # refine=false got no initial coarsening and refine=true coarsen=false got
    # one anyway. partrac has always had them as two blocks, which is why only
    # the interval gating was ever wrong there.
    d = tmp_path / "case"
    d.mkdir()
    shutil.copy(EXAMPLE, d / "expr_params.dat")
    r = subprocess.run([app("static_space_stepper"), str(d / "expr_params.dat")]
                       + SSS.split() + ["refine=" + refine, "coarsen=" + coarsen],
                       capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stdout + r.stderr
    assert ("Initial refinement" in r.stdout) == (refine == "true")
    assert ("Initial coarsening" in r.stdout) == (coarsen == "true")


# --- the time loops that were made parallel ------------------------------------

# These integrators draw no random numbers, so their result may not depend on
# how the work is divided. Each of them ran on one thread until the loop was
# parallelised, and each has a parallel sibling in the same file it was written
# to match; this is what says the match is right.
PARALLELISED = [
    # app, the loop that was serial
    ("tracers_triangleRK4", "Integrator_RK4::step"),
    ("filaments_triangleRK4", "Integrator_RK4::step"),
    ("tracertensors_triangleRK4", "Integrator_RK4::step_tensor"),
    ("static_space_stepper", "Integrator_Spatial::step_vec"),
]


def h5_datasets(d):
    """Every dataset in every h5 the run wrote, keyed by file and path."""
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


@pytest.mark.parametrize("name,loop", PARALLELISED, ids=[a[0] for a in PARALLELISED])
def test_the_result_does_not_depend_on_the_thread_count(
        name, loop, tmp_path, mesh_dir, felbm_dir, xdmf_dir):
    if not os.path.exists(app(name)):
        pytest.skip(name + " is not built")
    kind = dict((a[0], a[1]) for a in APPS)[name]
    args = dict((a[0], a[2]) for a in APPS)[name]
    work = ("Nrw=800 Nrw_max=8000 dt=0.005 dump_intv=0.05 stat_intv=1e9 "
            "checkpoint_intv=1e9 " +
            ("Ln=0.1 dxn=0.005 T=1e9" if name == "static_space_stepper" else "T=0.1"))

    out = {}
    for nthreads in (1, 4):
        parent = tmp_path / str(nthreads)
        parent.mkdir()
        d, cfg = case_dir(kind, parent, mesh_dir, felbm_dir, xdmf_dir)
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
