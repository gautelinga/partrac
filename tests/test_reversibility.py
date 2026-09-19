"""Refinement and coarsening under flow reversal: forward to T, back to 2T.

Reversing the flow is a sharp check on remeshing because the answers are closed
form. In plane Poiseuille a line laid along x has rho = sqrt(1 + (a t)^2) with
a = 3 u_inf x / R^2 on the way out; after the reversal the slope runs back
down as a(2T - t), so

    rho(2T) = 1                          the geometry returns exactly
    tau(2T) = 2 (T + a^2 T^3 / 3)        the compressed time does not

The strip stretches on the first leg, so refinement splits edges, and contracts
on the second, so coarsening merges them back. Any remeshing bookkeeping that
does not invert (for example a collapse that hands neighbours reference lengths
inconsistent with where the surviving node sits) shows up as a departure from
these values, with no reference run needed. The reversal is done by
checkpointing, negating u_inf in the expression file and restarting, so tau
must also survive the checkpoint.

The sine flow is reversed the same way and, unlike Poiseuille, folds the strip,
which is what makes coarsening merge edges whose compressed times differ.
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

T, U_INF, R, DT = 2.0, 1.0, 1.0, 0.001
BASE = ("mode=analytic init_mode=strip_x La=1.0 x0=0 y0=0 z0=0 Nrw=100 "
        "Nrw_max=200000 ds_max=0.02 refine=true refine_intv=0.01 Dm=0 "
        "int_order=1 dt=%g stat_intv=1e9 random=false seed=1 "
        "integrate_tau=true tau_intv=%g tau_max=0" % (DT, DT)).split()

needs_partrac = pytest.mark.skipif(not os.path.exists(PARTRAC),
                                   reason="partrac is not built")


def there_and_back(tmp_path, extra):
    """Advect to T, reverse u_inf by restart, return to 2T; return the resumed run's dump file."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    cfg = tmp_path / "expr_params.dat"
    cfg.write_text(open(POISEUILLE).read())

    def call(args):
        argv = [a for a in BASE if a.split("=")[0]
                not in {b.split("=")[0] for b in extra + args}] + extra + args
        r = subprocess.run([PARTRAC, str(cfg)] + argv,
                           capture_output=True, text=True, timeout=900)
        assert r.returncode == 0, r.stdout + r.stderr

    # the checkpoint lands one step past the stop time, so stop a step short:
    # the closed forms hold only if the flow reverses exactly at T
    call(["T=%g" % (T - DT), "dump_intv=1e9", "checkpoint_intv=%g" % T])
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
    """Per-edge relative errors in rho and tau against the closed forms, edge count, and sum of dl0."""
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
    """With refinement only, the strip returns to rho = 1 and tau matches the
    closed form, while the total reference length stays 1. Failure means edge
    splits corrupt dl0 or tau, biasing every stretching and mixing statistic of
    a refined run."""
    rho_err, tau_err, n, s0 = errors(
        there_and_back(tmp_path, ["coarsen=false", "ds_min=1e-12"]))
    assert n > 100                                   # it did refine
    assert s0 == pytest.approx(1.0, rel=1e-12)       # reference length conserved
    assert rho_err.max() < 1e-3
    assert tau_err.max() < 5e-2


@needs_partrac
def test_coarsening_is_reversible_too(tmp_path):
    """Coarsening on the return leg merges the strip back down and still returns
    rho = 1, exactly for edges never collapsed. A collapse must place the
    surviving node where the reference lengths it hands out say it went, or the
    strip comes back the wrong length."""
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
    """The tau error with coarsening is within a factor 2 of refinement alone.
    Merging discards the collapsed edge's tau, which is only acceptable while
    neighbouring edges agree on it."""
    plain = errors(there_and_back(tmp_path / "plain",
                                  ["coarsen=false", "ds_min=1e-12"]))
    merged = errors(there_and_back(tmp_path / "merged",
                                   ["coarsen=true", "coarsen_intv=0.01",
                                    "ds_min=0.008"]))
    assert merged[1].max() < 2 * plain[1].max()


# --- sine flow: reversed by negating u_inf and reversing the phases -------------
#
# The sine flow is a sequence of shears, so it reverses exactly: negate the
# amplitude, replay the phases backwards, and swap the direction pair so the
# shear that ran last is undone first. The expression file is written here, not
# taken from data_example, so the phase list is short and fixed.

# DT_S = TAU/8, so every half period ends exactly on a step
TAU, DT_S, N_HALF = 0.5, 0.0625, 2
CHI = [1.2154, 3.1199, 4.2865, 5.6534, 1.9023, 5.1624]

SINE_BASE = ("mode=analytic init_mode=strip_x La=0.5 x0=0.5 y0=0.5 z0=0.5 "
             "Nrw=200 Nrw_max=2000000 Dm=0 int_order=1 dt=%g stat_intv=1e9 "
             "random=false seed=1 integrate_tau=true tau_intv=%g tau_max=0"
             % (DT_S, DT_S)).split()


def write_sine(path, chi, u_inf, flowdir, depdir):
    """Write a sine_flow expression file with the given phases, amplitude and shear directions."""
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
    """Shear through n_half half periods and back; return size and tau at the turn and at the end."""
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

    # the checkpoint lands one step past the stop time, so stop a step short
    # of the half-period boundary
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
    """Without remeshing, a strip sheared out and back returns to its reference
    lengths to 1e-12 and ends with exactly twice the tau it had at the turn. This
    validates the reversal itself, so the remeshing tests below measure only the
    cost of remeshing."""
    half, whole = sine_there_and_back(tmp_path, N_HALF,
                                      ["ds_max=1e9", "ds_min=1e-12",
                                       "refine=false", "coarsen=false"])
    assert half["w"].sum() / half["w0"].sum() > 2         # it really stretched
    assert np.max(np.abs(whole["w"] / whole["w0"] - 1)) < 1e-12
    assert whole["tau"] == pytest.approx(2 * half["tau"], rel=1e-9)


@needs_partrac
def test_a_sheet_reverses_exactly_too(tmp_path):
    """An x-z sheet reverses exactly like the strip, with tau doubling on its
    faces. The flow leaves z alone, so the sheet is the strip extruded and its
    faces must integrate tau from dA/dA0 as edges do from dl/dl0."""
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
    """After six half periods out and back, coarsening that merges across cusps
    conserves total reference length, doubles the dl0-weighted tau, and keeps the
    reported scalar within 0.2% of an uncoarsened reference over two decades of
    Peclet number. This pins the merge rule that conserves sum dl0/sqrt(tau)."""
    # six half periods fold the strip into cusps, where neighbouring edges
    # disagree on tau; the uncoarsened run is the reference
    mesh = ["ds_max=0.005", "refine=true", "refine_intv=%g" % DT_S]
    half, whole = sine_there_and_back(tmp_path / "merged", 6, mesh
                                      + ["ds_min=0.002", "coarsen=true",
                                         "coarsen_intv=%g" % DT_S])
    ref_half, ref = sine_there_and_back(tmp_path / "plain", 6, mesh
                                        + ["ds_min=1e-12", "coarsen=false"])

    assert len(whole["w"]) < len(ref["w"]) / 4            # it really merged
    for s in (half, whole, ref_half, ref):
        assert s["w0"].sum() == pytest.approx(0.5, rel=1e-12)

    # tau must still double, with and without merging
    moment = lambda s: (s["w0"] * s["tau"]).sum()
    assert moment(ref) / moment(ref_half) == pytest.approx(2.0, rel=0.02)
    assert moment(whole) / moment(half) == pytest.approx(2.0, rel=0.02)

    # the 0.2% tolerance is tight enough that alternative merge rules (a
    # dl0-weighted mean of tau, or dropping the removed edge's tau) fail it
    # at the well-mixed end
    for k in (4e-3, 4e-2, 4e-1):
        c = lambda s: (s["w0"] / np.sqrt(1 + 4 * k * s["tau"])).sum()
        assert c(whole) == pytest.approx(c(ref), rel=0.002), k


@pytest.fixture(scope="module")
def folded_sheet(tmp_path_factory):
    """Sizes and tau of an x-z sheet at the turn and end of four sine half periods, coarsened and not.

    Four half periods fold the sheet enough that the stars of collapsed edges
    are not planar.
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
    """Total reference area stays 0.25 and every dA0 stays positive through
    coarsening of a folded sheet. A collapse on a non-planar star loses current
    area, but losing reference area would bias every area-weighted statistic."""
    for s in folded_sheet:
        assert s["w0"].sum() == pytest.approx(0.25, rel=1e-12)
        assert s["w0"].min() > 0                          # never through zero


@needs_partrac
def test_coarsening_a_folded_sheet_keeps_the_scalar_it_reports(folded_sheet):
    """On a folded sheet the reported scalar survives coarsening to within 0.5%
    of the uncoarsened run, with an error nearly independent of Peclet number
    over four decades. tau stays finite and positive on every face."""
    half, whole, ref_half, ref = folded_sheet
    assert np.isfinite(whole["tau"]).all() and whole["tau"].min() > 0
    c = lambda s, k: (s["w0"] / np.sqrt(1 + 4 * k * s["tau"])).sum()
    err = [c(whole, k) / c(ref, k) - 1 for k in (4e-3, 4e-2, 4e-1, 4.0, 40.0)]
    assert max(abs(e) for e in err) < 0.005
    assert max(err) - min(err) < 0.0015
