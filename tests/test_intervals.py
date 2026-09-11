"""What happens when the remeshing intervals do not line up.

`refine_intv`, `coarsen_intv`, `filter_intv`, `inject_intv`, `tau_intv` and
`checkpoint_intv` are independent step counts through one time loop, and the
loop runs them in a fixed order: inject, refine, coarsen, filter, cull, write,
advect, integrate tau. Every combination of phases is reachable, and the ones
that bit before were all mismatches -- one part flagging an entity another had
assumed would still be there. The suite mostly set them equal, so this file
sets them against each other.

The flow is Hagen-Poiseuille throughout, because it is the one where mismatch
has a closed form to be wrong about: a steady parallel flow does not stretch
its streak surface, so dA/dA0 = 1 identically and the reference area is fixed
where it is created. Two exact laws follow, and both are phase-sensitive.

    The swept area tracks the last injection, not the clock. Injection lays
    down one interval's worth of area per firing, so at step `it` the sheet
    carries `floor(it/n) * n * dt` of sweeping, where `n` is the injection
    interval in steps. That is `t` only when `t` is an injection time. Nothing
    else in the loop may change it: refinement splits a face into two that
    share its dA0, coarsening hands a collapsed face's dA0 to its neighbours.

    A face lives `tau_max` and is replaced every `inject_intv`. So the sheet
    outlives the cull exactly when `tau_max >= inject_intv`, and below that it
    is culled to nothing between one injection and the next -- which must end
    the run cleanly, not in a crash or a file of NaN.

Runs here are a few hundred steps and take about a tenth of a second each.
"""

import os
import shutil
import subprocess

import h5py
import numpy as np
import pytest

from paths import REPO, app

PARTRAC = app("partrac")
HAGEN = os.path.join(REPO, "data_example", "hagen_poiseuille", "expr_params.dat")

U_INF, R, NRW, DT = 1.0, 1.0, 41, 0.005

BASE = ("mode=analytic init_mode=uniform_x x0=0 y0=0 z0=0 Nrw=%d Nrw_max=200000 "
        "inject=true inject_edges=true T_inject=1e10 "
        "ds_max=1e9 ds_min=1e-9 refine=false coarsen=false Dm=0 int_order=1 "
        "dt=%g stat_intv=1e9 dump_intv=1e9 checkpoint_intv=1e9 random=false "
        "seed=1" % (NRW, DT)).split()

# a refinement that actually splits the sheet without touching the inlet edges,
# which span the pipe in steps of 2R/(NRW-1) = 0.05
REMESH = ["ds_max=0.1", "ds_min=0.02", "refine=true", "coarsen=true"]

needs_partrac = pytest.mark.skipif(not os.path.exists(PARTRAC),
                                   reason="partrac is not built")


def steps_per(intv, dt=DT):
    """partrac's own rounding: an interval shorter than a step is every step."""
    n = intv / dt
    return 1 if not n > 1. else int(n)


def sweep_rate(n=NRW):
    """The rate the inlet sweeps area, as its straight edges measure it."""
    h = 2 * R / (n - 1)
    return 8. / 3. * U_INF * R - 2. / 3. * U_INF * h ** 2


def dumps_with_faces(inject_intv, dump_intv, T):
    """How many dumps fall on or after the first injection."""
    first = steps_per(inject_intv) * DT
    n = steps_per(dump_intv)
    return len([it for it in range(0, int(round(T / DT)) + 1, n)
                if it * DT >= first - DT / 2])


def swept_by(t, inject_intv):
    """Area laid down by time t: one interval's worth per injection, no more."""
    it = int(round(t / DT))
    n = steps_per(inject_intv)
    return sweep_rate() * (it // n) * n * DT


def run(tmp_path, extra):
    # partrac refuses a parameter given twice, and these cases are built by
    # layering one interval over a set of them, so the last word wins here
    tmp_path.mkdir(parents=True, exist_ok=True)
    shutil.copy(HAGEN, tmp_path / "expr_params.dat")
    argv = {}
    for a in BASE + extra:
        argv[a.split("=")[0]] = a
    return subprocess.run([PARTRAC, str(tmp_path / "expr_params.dat")]
                          + list(argv.values()),
                          capture_output=True, text=True, timeout=600)


def series(tmp_path, extra):
    """Every dump that has faces, as (t, dA, dA0)."""
    r = run(tmp_path, extra)
    assert r.returncode == 0, r.stdout + r.stderr
    dump = list(tmp_path.rglob("data_from_t*.h5"))
    assert len(dump) == 1
    h = h5py.File(dump[0], "r")
    out = []
    for k in sorted(h.keys(), key=float):
        if "dA0" not in h[k]:
            continue                   # before the first injection
        out.append((float(k),
                    np.array(h[k + "/dA"]).ravel(),
                    np.array(h[k + "/dA0"]).ravel()))
    return out


def no_nan(tmp_path):
    """Neither the dumps nor the statistics may carry a NaN or an infinity."""
    for f in tmp_path.rglob("data_from_t*.h5"):
        h = h5py.File(f, "r")
        for key in h:
            for name in h[key]:
                a = np.array(h[key + "/" + name])
                if a.dtype.kind == "f":
                    assert np.isfinite(a).all(), "%s at t = %s" % (name, key)
    for f in tmp_path.rglob("tdata_from_t*.dat"):
        text = f.read_text().lower()
        assert "nan" not in text and "inf" not in text, f


# --- injection against refinement and coarsening -------------------------------

# in steps of dt = 0.005 these are 4, 6, 10 and 14, so no pair of them shares a
# phase for long: 0.02 and 0.03 first coincide at 0.06, 0.05 and 0.07 at 0.35
@needs_partrac
@pytest.mark.parametrize("inject_intv,refine_intv,coarsen_intv", [
    (0.05, 0.03, 0.07),     # remesh faster and slower than the injection
    (0.05, 0.07, 0.02),
    (0.03, 0.05, 0.05),     # inject faster than either
    (0.03, 0.02, 0.07),
    (0.07, 0.03, 0.05),     # inject slower than either
    (0.02, 0.07, 0.03),
])
def test_remeshing_out_of_phase_moves_no_swept_area(tmp_path, inject_intv,
                                                    refine_intv, coarsen_intv):
    # dA0 is laid down at the inlet and touched by nothing afterwards, so the
    # total is fixed by the injections alone however the remeshing is phased
    s = series(tmp_path, REMESH + [
        "inject_intv=%g" % inject_intv, "refine_intv=%g" % refine_intv,
        "coarsen_intv=%g" % coarsen_intv, "T=0.3", "dump_intv=0.05"])
    assert len(s) == dumps_with_faces(inject_intv, 0.05, 0.3)
    assert len(s[-1][2]) > 400                      # it really remeshed
    for t, dA, dA0 in s:
        assert dA0.sum() == pytest.approx(swept_by(t, inject_intv), rel=1e-12)
        assert dA0.min() > 0                        # no face without area
        assert np.isfinite(dA).all()


@needs_partrac
@pytest.mark.parametrize("inject_intv", [0.02, 0.03, 0.07])
def test_the_swept_area_follows_the_injections_not_the_clock(tmp_path, inject_intv):
    # the dumps are on 0.05 and the injections are not, so most dumps land
    # between two generations and see the area of the earlier one
    s = series(tmp_path, ["inject_intv=%g" % inject_intv, "T=0.3",
                          "dump_intv=0.05"])
    off = [t for t, _, _ in s
           if int(round(t / DT)) % steps_per(inject_intv) != 0]
    assert off, "this case never lands off an injection"
    for t, _, dA0 in s:
        assert dA0.sum() == pytest.approx(swept_by(t, inject_intv), rel=1e-12)


# --- the exit-plane cull, which runs on filter_intv ----------------------------

@needs_partrac
@pytest.mark.parametrize("inject_intv,filter_intv", [
    (0.05, 0.02),           # cull faster than the injection
    (0.05, 0.07),           # and slower
    (0.03, 0.05),
    (0.07, 0.02),
])
def test_a_cull_out_of_phase_does_not_starve_the_inlet(tmp_path, inject_intv,
                                                       filter_intv):
    # the inlet spans the pipe and half of it is beyond the plane from the
    # start, so every cull is a chance to remove material the next injection
    # was going to stitch to. It has to survive whatever the phase.
    s = series(tmp_path, REMESH + [
        "inject_intv=%g" % inject_intv, "filter_intv=%g" % filter_intv,
        "refine_intv=0.03", "coarsen_intv=0.05", "exit_plane=x", "Ln=0.5",
        "T=0.4", "dump_intv=0.1", "stat_intv=0.05"])
    assert len(s) == 4
    total = [dA0.sum() for _, _, dA0 in s]
    assert all(b > a for a, b in zip(total, total[1:]))   # still being fed
    assert all(dA0.min() > 0 for _, _, dA0 in s)
    no_nan(tmp_path)


# --- the tau cull, whose interval is a time rather than a step count -----------

@needs_partrac
@pytest.mark.parametrize("inject_intv", [0.03, 0.05, 0.07])
def test_a_sheet_outlives_the_tau_cull_when_tau_max_reaches_the_next_generation(
        tmp_path, inject_intv):
    # in a parallel flow rho is 1, so tau is just the age of a face and the
    # oldest generation is exactly one injection interval old
    s = series(tmp_path, REMESH + [
        "inject_intv=%g" % inject_intv, "refine_intv=0.03", "coarsen_intv=0.07",
        "integrate_tau=true", "tau_intv=%g" % DT,
        "tau_max=%g" % (2 * inject_intv), "T=0.4", "dump_intv=0.1"])
    assert len(s) == 4
    for t, _, dA0 in s:
        assert dA0.min() > 0
    no_nan(tmp_path)


@needs_partrac
@pytest.mark.parametrize("inject_intv", [0.03, 0.05, 0.07])
def test_a_tau_cull_the_injection_cannot_outrun_stops_the_run(tmp_path,
                                                              inject_intv):
    # below one injection interval every face is culled before its replacement
    # arrives, and the sheet loses its dimension. That is a stop, not a crash,
    # and not a file of rows nobody can interpret
    r = run(tmp_path, REMESH + [
        "inject_intv=%g" % inject_intv, "refine_intv=0.03", "coarsen_intv=0.07",
        "integrate_tau=true", "tau_intv=%g" % DT,
        "tau_max=%g" % (0.5 * inject_intv), "T=0.4", "dump_intv=0.1"])
    assert r.returncode == 1
    assert "changed dimension" in r.stdout + r.stderr
    no_nan(tmp_path)


# --- the intervals that only write ---------------------------------------------

OUTPUT_ONLY = [
    ("dump_intv=0.1", "stat_intv=0.05", "checkpoint_intv=1e9"),
    ("dump_intv=%g" % DT, "stat_intv=0.05", "checkpoint_intv=1e9"),
    ("dump_intv=0.03", "stat_intv=0.05", "checkpoint_intv=1e9"),
    ("dump_intv=0", "stat_intv=0.05", "checkpoint_intv=1e9"),
    ("dump_intv=0.1", "stat_intv=%g" % DT, "checkpoint_intv=1e9"),
    ("dump_intv=0.1", "stat_intv=0", "checkpoint_intv=1e9"),
    ("dump_intv=0.1", "stat_intv=0.05", "checkpoint_intv=0.02"),
    ("dump_intv=0.1", "stat_intv=0.05", "checkpoint_intv=0.07"),
]


def final_checkpoint(tmp_path):
    """The checkpoint the run ends on, file by file."""
    cp = list(tmp_path.rglob("positions.pos"))
    assert len(cp) == 1, cp
    return {f.name: f.read_bytes()
            for f in sorted(cp[0].parent.iterdir()) if f.is_file()}


@needs_partrac
@pytest.mark.parametrize("io", OUTPUT_ONLY[1:], ids=lambda c: "_".join(c))
def test_writing_more_often_does_not_move_the_mesh(tmp_path, io):
    # dump_intv is read by the time loop twice: it dumps, and it is one of the
    # three conditions that call compute_interior. That second reading is what
    # makes this worth asserting -- an output interval steering a computation
    # is one curvature parameter away from steering the refinement with it.
    phys = REMESH + ["inject_intv=0.05", "refine_intv=0.03", "coarsen_intv=0.07",
                     "integrate_tau=true", "tau_intv=%g" % DT, "tau_max=0",
                     "T=0.4"]
    ref = run(tmp_path / "ref", phys + list(OUTPUT_ONLY[0]))
    assert ref.returncode == 0, ref.stdout + ref.stderr
    got = run(tmp_path / "got", phys + list(io))
    assert got.returncode == 0, got.stdout + got.stderr
    a, b = final_checkpoint(tmp_path / "ref"), final_checkpoint(tmp_path / "got")
    assert set(a) == set(b)
    for name in sorted(a):
        if name.endswith(".dat"):
            continue                   # the parameter dump holds the intervals
        assert a[name] == b[name], name


# --- resuming off every phase at once ------------------------------------------

@needs_partrac
def test_a_resume_off_the_phase_of_every_interval_is_identical(tmp_path):
    # test_restart covers a resume off the refinement phase; this is the same
    # question for the injecting path, where the checkpoint also has to carry
    # the inlet. The final checkpoint is written one step past T, so stopping
    # at 0.17 resumes at step 35 -- which is 5 past a refinement (6 steps), 5
    # past an injection (10) and 7 past a coarsening (14).
    phys = REMESH + ["inject_intv=0.05", "refine_intv=0.03", "coarsen_intv=0.07",
                     "integrate_tau=true", "tau_intv=%g" % DT, "tau_max=0",
                     "dump_intv=0.1"]
    end, stop = 0.4, 0.17
    cont, split = tmp_path / "cont", tmp_path / "split"
    for case, extra in ((cont, ["T=%g" % (end + DT / 2), "checkpoint_intv=1e9"]),
                        (split, ["T=%g" % stop, "checkpoint_intv=%g" % stop])):
        r = run(case, phys + extra)
        assert r.returncode == 0, r.stdout + r.stderr
    checkpoint = list(split.rglob("edges.edge"))
    assert len(checkpoint) == 1
    r = run(split, phys + ["T=%g" % (end + DT / 2), "checkpoint_intv=1e9",
                           "restart_folder=" + str(checkpoint[0].parent.parent)])
    assert r.returncode == 0, r.stdout + r.stderr

    def state_at(case, t):
        for f in sorted(case.rglob("data_from_t*.h5")):
            h = h5py.File(f, "r")
            for key in h:
                if abs(float(key) - t) < 1e-9:
                    return {n: np.array(h[key + "/" + n]) for n in h[key]}
        raise AssertionError("no dump at t = %g under %s" % (t, case))

    a, b = state_at(cont, end), state_at(split, end)
    assert len(a["dA0"]) > 400                  # it really remeshed
    assert set(a) == set(b)
    for name in sorted(a):
        assert np.array_equal(a[name], b[name]), name


# --- intervals at their limits --------------------------------------------------

@needs_partrac
@pytest.mark.parametrize("extra", [
    ["inject_intv=1e-9"],                                   # every step
    ["inject_intv=0.05", "refine_intv=1e-9"],
    ["inject_intv=0.05", "coarsen_intv=1e-9"],
    ["inject_intv=0.05", "tau_intv=1e-9", "integrate_tau=true", "tau_max=0"],
    ["inject_intv=0.05", "refine_intv=0"],                  # off
    ["inject_intv=0.05", "coarsen_intv=0"],
    ["inject_intv=0.05", "checkpoint_intv=0"],
], ids=lambda e: e[-1].replace("=", ""))
def test_an_interval_shorter_than_a_step_or_switched_off_still_runs(tmp_path, extra):
    # steps_per never returns zero, so a sub-step interval is every step; an
    # interval of 0 is off. Both ends have to leave the mesh intact, which is
    # what the every-step end tests -- it remeshes 80 times in 80 steps
    s = series(tmp_path, REMESH + ["refine_intv=0.03", "coarsen_intv=0.07",
                                   "T=0.4", "dump_intv=0.1"] + extra)
    assert s, "nothing was ever injected"
    for t, dA, dA0 in s:
        assert dA0.min() > 0
        assert np.isfinite(dA).all()
    no_nan(tmp_path)
