"""Run plane Poiseuille forwards, then backwards, and see what remeshing cost.

Reversing the flow at time T is the sharpest check there is on refinement and
coarsening, because both answers are closed form. A line laid along x has
rho = sqrt(1 + (a t)^2) with a = 3 u_inf x / R^2 on the way out; after the
reversal the slope runs back down as a(2T - t), so

    rho(2T) = 1                          the geometry returns exactly
    tau(2T) = 2 (T + a^2 T^3 / 3)        the compressed time does not

The strip stretches on the first leg, so refinement splits edges; it contracts
on the second, so coarsening merges them back. Any bookkeeping that does not
invert shows up as a departure from those two values, with no reference run
needed. This found a collapse that left the surviving node where it was while
handing its neighbours the reference length of a move that never happened.

The reversal is done by checkpointing, negating u_inf in the expression file,
and restarting -- so it also exercises tau surviving a checkpoint.
"""

import os
import re
import subprocess

import h5py
import numpy as np
import pytest

from paths import REPO, app

PARTRAC = app("partrac")
POISEUILLE = os.path.join(REPO, "data_example", "plane_poiseuille", "expr_params.dat")

T, U_INF, R = 2.0, 1.0, 1.0
BASE = ("mode=analytic init_mode=strip_x La=1.0 x0=0 y0=0 z0=0 Nrw=100 "
        "Nrw_max=200000 ds_max=0.02 refine=true refine_intv=0.01 Dm=0 "
        "int_order=1 dt=0.001 stat_intv=1e9 random=false seed=1 "
        "integrate_tau=true tau_intv=0.001 tau_max=0").split()

needs_partrac = pytest.mark.skipif(not os.path.exists(PARTRAC),
                                   reason="partrac is not built")


def there_and_back(tmp_path, extra):
    """Advect to T, reverse the flow, come back; return the final dump."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    cfg = tmp_path / "expr_params.dat"
    cfg.write_text(open(POISEUILLE).read())

    def call(args):
        argv = [a for a in BASE if a.split("=")[0]
                not in {b.split("=")[0] for b in extra + args}] + extra + args
        r = subprocess.run([PARTRAC, str(cfg)] + argv,
                           capture_output=True, text=True, timeout=900)
        assert r.returncode == 0, r.stdout + r.stderr

    call(["T=%g" % T, "dump_intv=1e9", "checkpoint_intv=%g" % T])
    checkpoint = list(tmp_path.rglob("edges.edge"))
    assert len(checkpoint) == 1
    cfg.write_text(re.sub(r"(?m)^u_inf=.*$", "u_inf=%g" % -U_INF, cfg.read_text()))
    call(["T=%g" % (2 * T + 0.002), "dump_intv=0.5", "checkpoint_intv=1e9",
          "restart_folder=" + str(checkpoint[0].parent.parent)])

    resumed = [f for f in sorted(tmp_path.rglob("data_from_t*.h5"))
               if "t0.000000" not in f.name]
    assert len(resumed) == 1
    return h5py.File(resumed[0], "r")


def errors(dump):
    key = sorted(dump.keys(), key=float)[-1]
    t = float(key)
    assert t > 2 * T - 0.01, "the return leg did not finish"
    dl = np.array(dump[key + "/dl"]).ravel()
    dl0 = np.array(dump[key + "/dl0"]).ravel()
    tau = np.array(dump[key + "/tau"]).ravel()
    e = np.array(dump[key + "/edges"])
    p = np.array(dump[key + "/points"])
    a = 3 * U_INF * 0.5 * (p[e[:, 0], 0] + p[e[:, 1], 0]) / R ** 2

    back = t - T
    rho_exact = np.sqrt(1 + (a * (2 * T - t)) ** 2)
    tau_exact = ((T + a ** 2 * T ** 3 / 3)
                 + (back + a ** 2 * (T ** 3 - (T - back) ** 3) / 3))
    return (np.abs(dl / dl0 / rho_exact - 1), np.abs(tau / tau_exact - 1),
            len(dl), dl0.sum())


@needs_partrac
def test_refinement_alone_is_reversible(tmp_path):
    rho_err, tau_err, n, s0 = errors(
        there_and_back(tmp_path, ["coarsen=false", "ds_min=1e-12"]))
    assert n > 100                                   # it did refine
    assert s0 == pytest.approx(1.0, rel=1e-12)       # and conserved the measure
    assert rho_err.max() < 1e-3
    assert tau_err.max() < 5e-2


@needs_partrac
def test_coarsening_is_reversible_too(tmp_path):
    # the collapse must put the surviving node where the reference lengths it
    # hands out say it went, or the strip comes back the wrong length
    rho_err, tau_err, n, s0 = errors(
        there_and_back(tmp_path, ["coarsen=true", "coarsen_intv=0.01",
                                  "ds_min=0.008"]))
    assert n < 120                                   # it did merge back down
    assert s0 == pytest.approx(1.0, rel=1e-12)
    assert np.median(rho_err) < 1e-12                # exact for untouched edges
    assert rho_err.max() < 1e-4
    assert tau_err.max() < 5e-2


@needs_partrac
def test_coarsening_costs_no_more_than_refinement_alone(tmp_path):
    # merging discards the collapsed edge's tau, which is only defensible while
    # neighbours agree on it; if that ever stops holding, this is where it shows
    plain = errors(there_and_back(tmp_path / "plain",
                                  ["coarsen=false", "ds_min=1e-12"]))
    merged = errors(there_and_back(tmp_path / "merged",
                                   ["coarsen=true", "coarsen_intv=0.01",
                                    "ds_min=0.008"]))
    assert merged[1].max() < 2 * plain[1].max()


# --- the same trick on the flow that actually folds ---------------------------
#
# The sine flow is a sequence of shears, so it reverses exactly too: negate the
# amplitude, replay the phases backwards, and swap the direction pair so the
# shear that ran last is undone first. Unlike plane Poiseuille it folds, which
# is the only way to make coarsening merge edges whose compressed times differ.
# The expression file is written here rather than taken from data_example, so
# the phase list is short enough to reverse by hand and cannot drift.

TAU, DT_S, N_HALF = 0.5, 0.0625, 2
CHI = [1.2154, 3.1199, 4.2865, 5.6534, 1.9023, 5.1624]

SINE_BASE = ("mode=analytic init_mode=strip_x La=0.5 x0=0.5 y0=0.5 z0=0.5 "
             "Nrw=200 Nrw_max=2000000 Dm=0 int_order=1 dt=%g stat_intv=1e9 "
             "random=false seed=1 integrate_tau=true tau_intv=%g tau_max=0"
             % (DT_S, DT_S)).split()


def write_sine(path, chi, u_inf, flowdir, depdir):
    keys = dict(t_min=0.0, t_max=1e7, x_min=0.0, y_min=0.0, z_min=0.0,
                x_max=1.0, y_max=1.0, z_max=1.0, Lx=1.0, Ly=1.0, Lz=1.0,
                expression="sine_flow", p_inf=1.0, rho=1.0, u_inf=u_inf,
                tau=TAU)
    lines = ["%s=%s" % (k, v) for k, v in keys.items()]
    lines += ["chi=" + ",".join("%.12g" % c for c in chi),
              "flowdir=" + ",".join(map(str, flowdir)),
              "depdir=" + ",".join(map(str, depdir))]
    path.write_text("\n".join(lines) + "\n")


def sine_there_and_back(tmp_path, n_half, extra):
    """Shear forwards through n_half half periods, then undo them."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    cfg = tmp_path / "expr_params.dat"
    chi = CHI[:n_half]
    u = 0.70710678118

    def call(args):
        argv = [a for a in SINE_BASE if a.split("=")[0]
                not in {b.split("=")[0] for b in extra + args}] + extra + args
        r = subprocess.run([PARTRAC, str(cfg)] + argv,
                           capture_output=True, text=True, timeout=900)
        assert r.returncode == 0, r.stdout + r.stderr

    # the checkpoint lands one step past T, so stop a step short of the boundary
    forward = n_half * TAU - DT_S
    write_sine(cfg, chi, u, (1, 0), (0, 1))
    call(["T=%g" % forward, "dump_intv=1e9",
          "checkpoint_intv=%g" % forward])
    checkpoint = list(tmp_path.rglob("edges.edge"))
    assert len(checkpoint) == 1

    write_sine(cfg, list(reversed(chi)), -u, (0, 1), (1, 0))
    call(["T=%g" % (2 * n_half * TAU + DT_S / 2),
          "dump_intv=%g" % (n_half * TAU),
          "checkpoint_intv=1e9",
          "restart_folder=" + str(checkpoint[0].parent.parent)])

    resumed = [f for f in sorted(tmp_path.rglob("data_from_t*.h5"))
               if "t0.000000" not in f.name]
    assert len(resumed) == 1
    d = h5py.File(resumed[0], "r")
    keys = sorted(d.keys(), key=float)
    assert float(keys[-1]) == pytest.approx(2 * n_half * TAU, abs=1e-9)
    # a sheet reports areas where a strip reports lengths
    size = "dA" if "dA" in d[keys[0]] else "dl"
    at = lambda k: {"w": np.array(d[k + "/" + size]).ravel(),
                    "w0": np.array(d[k + "/" + size + "0"]).ravel(),
                    "tau": np.array(d[k + "/tau"]).ravel()}
    return at(keys[0]), at(keys[-1])


@needs_partrac
def test_the_sine_flow_reverses_exactly(tmp_path):
    half, whole = sine_there_and_back(tmp_path, N_HALF,
                                      ["ds_max=1e9", "ds_min=1e-12",
                                       "refine=false", "coarsen=false"])
    assert half["w"].sum() / half["w0"].sum() > 2         # it really folded
    assert np.max(np.abs(whole["w"] / whole["w0"] - 1)) < 1e-12
    assert whole["tau"] == pytest.approx(2 * half["tau"], rel=1e-9)


@needs_partrac
def test_a_sheet_reverses_exactly_too(tmp_path):
    # the flow leaves z alone, so an x-z sheet is the x-strip extruded and its
    # faces integrate tau on dA/dA0 as the edges do on ds/ds0
    half, whole = sine_there_and_back(tmp_path, N_HALF,
                                      ["init_mode=sheet_xz", "Lb=0.5",
                                       "ds_init=0.05", "ds_max=1e9",
                                       "ds_min=1e-12", "refine=false",
                                       "coarsen=false"])
    assert half["w"].sum() / half["w0"].sum() > 2
    assert np.max(np.abs(whole["w"] / whole["w0"] - 1)) < 1e-12
    assert whole["tau"] == pytest.approx(2 * half["tau"], rel=1e-9)


@needs_partrac
def test_merging_across_cusps_preserves_the_aggregates(tmp_path):
    # six half periods fold the strip enough to make cusps, where neighbouring
    # edges finally disagree about tau; coarsening then merges across them. The
    # uncoarsened run is the reference: a collapse currently drops the removed
    # edge's tau outright, and this is what bounds what that costs.
    mesh = ["ds_max=0.005", "refine=true", "refine_intv=%g" % DT_S]
    half, whole = sine_there_and_back(tmp_path / "merged", 6, mesh
                                      + ["ds_min=0.002", "coarsen=true",
                                         "coarsen_intv=%g" % DT_S])
    ref_half, ref = sine_there_and_back(tmp_path / "plain", 6, mesh
                                        + ["ds_min=1e-12", "coarsen=false"])

    assert len(whole["w"]) < len(ref["w"]) / 4            # it really merged
    for s in (half, whole, ref_half, ref):
        assert s["w0"].sum() == pytest.approx(0.5, rel=1e-12)

    # tau must still double, both with and without the merging
    moment = lambda s: (s["w0"] * s["tau"]).sum()
    assert moment(ref) / moment(ref_half) == pytest.approx(2.0, rel=0.02)
    assert moment(whole) / moment(half) == pytest.approx(2.0, rel=0.02)

    # and the scalar variance the method reports must survive the merging,
    # over two decades of Peclet. The tolerance is what pins the merge rule: a
    # collapse conserves the sum of ds0/sqrt(tau), and the alternatives -- a
    # ds0-weighted mean of tau, or dropping the removed edge's tau outright --
    # reach 0.3% and 0.4% at the mixed end, so they fail here rather than
    # passing quietly.
    for k in (4e-3, 4e-2, 4e-1):
        c = lambda s: (s["w0"] / np.sqrt(1 + 4 * k * s["tau"])).sum()
        assert c(whole) == pytest.approx(c(ref), rel=0.002), k


@pytest.fixture(scope="module")
def folded_sheet(tmp_path_factory):
    """Four half periods of the sine flow on an x-z sheet, coarsened and not.

    Four is enough folding that the stars of the collapsed edges are properly
    non-planar, which is the case the dA0 share-out has to get right, and
    that faces disagree about tau by the time the way back merges them.
    """
    if not os.path.exists(PARTRAC):
        pytest.skip("partrac is not built")
    root = tmp_path_factory.mktemp("folded_sheet")
    sheet = ["init_mode=sheet_xz", "Lb=0.5", "ds_init=0.05"]
    mesh = ["ds_max=0.02", "refine=true", "refine_intv=%g" % DT_S]
    half, whole = sine_there_and_back(root / "merged", 4, sheet + mesh
                                      + ["ds_min=0.008", "coarsen=true",
                                         "coarsen_intv=%g" % DT_S])
    ref_half, ref = sine_there_and_back(root / "plain", 4, sheet + mesh
                                        + ["ds_min=1e-12", "coarsen=false"])
    assert half["w"].sum() / half["w0"].sum() > 2         # it really folded
    assert len(whole["w"]) < len(ref["w"]) / 4           # and really merged
    return half, whole, ref_half, ref


@needs_partrac
def test_coarsening_a_folded_sheet_conserves_its_reference_area(folded_sheet):
    # A collapse on a folded sheet really does destroy area: the star of the
    # edge is not planar, so the survivors do not grow by what the removed
    # faces held, and the difference runs about 0.4% per collapse. It must not
    # destroy material with it. Handing the removed dA0 out at the observed
    # area gain lost 0.2% of the sheet here and 3.6% over a longer fold, which
    # is a silent bias on every statistic, since all of them weight by dA0.
    for s in folded_sheet:
        assert s["w0"].sum() == pytest.approx(0.25, rel=1e-12)
        assert s["w0"].min() > 0                          # never through zero


@needs_partrac
def test_coarsening_a_folded_sheet_keeps_the_scalar_it_reports(folded_sheet):
    # A collapse shares the removed faces' variance content dA0/sqrt(tau) out
    # with their reference area, so the scalar the method reports survives it.
    # The tolerance is the cost of coarsening itself, which no merge rule
    # removes. What pins the rule is that the error is the same at every
    # Peclet number: dropping the removed tau instead drifts from 0.04% to
    # 0.34% across these k, a spread of 0.30%, and fails here.
    half, whole, ref_half, ref = folded_sheet
    assert np.isfinite(whole["tau"]).all() and whole["tau"].min() > 0
    c = lambda s, k: (s["w0"] / np.sqrt(1 + 4 * k * s["tau"])).sum()
    err = [c(whole, k) / c(ref, k) - 1 for k in (4e-3, 4e-2, 4e-1, 4.0, 40.0)]
    assert max(abs(e) for e in err) < 0.005
    assert max(err) - min(err) < 0.0015
