"""Every shipped analytic example runs (mode=analytic, AnalyticInterpol).

Each example whose expr_params.dat names an analytic expression must run a few
steps to a stats file of finite numbers with particles left in the domain.
"""

import os
import pathlib

import numpy as np
import pytest

from dumps import read_stats
from paths import REPO, app
from runs import copy_example, run_app

PARTRAC = app("partrac")


def analytic_examples():
    """The names of the shipped example folders whose expr_params.dat names an analytic expression."""
    root = pathlib.Path(REPO, "data_example")
    return [f.parent.name for f in sorted(root.glob("*/expr_params.dat"))
            if "\nexpression=" in "\n" + f.read_text()]


# brinkman_cylinder tabulates its profile on 10^6 points at startup
EXAMPLES = [pytest.param(c, marks=pytest.mark.slow) if c == "brinkman_cylinder" else c
            for c in analytic_examples()]


@pytest.mark.skipif(not os.path.exists(PARTRAC), reason="partrac is not built")
@pytest.mark.parametrize("case", EXAMPLES)
def test_every_shipped_analytic_example_runs(case, tmp_path):
    """Every analytic example runs two steps and writes finite statistics with particles left.

    A user starting from a shipped example should get a working run, not one
    whose particles all start outside the domain.
    """
    src = os.path.join(REPO, "data_example", case, "expr_params.dat")
    with open(src) as f:
        keys = dict(l.strip().split("=", 1) for l in f if "=" in l)
    # start the short line at the centre of the example's own domain
    centre = ["%s0=%g" % (k, (float(keys[k + "_min"]) + float(keys[k + "_max"])) / 2)
              for k in "xyz"]
    run_app(PARTRAC, copy_example(src, tmp_path),
            "mode=analytic init_mode=uniform_x La=0.2 Nrw=20 Nrw_max=20 ds_max=1e9 "
            "ds_min=1e-12 refine=false coarsen=false Dm=0 int_order=1 dt=0.01 T=0.02 "
            "stat_intv=0.01 dump_intv=1e9 checkpoint_intv=1e9 random=false seed=1",
            centre, timeout=300)

    # exiting 0 is not enough: with every particle outside the domain partrac
    # writes a stats file of NaN and still succeeds; read_stats checks every
    # row against the header
    st = read_stats(tmp_path)
    for name, values in st.items():
        assert np.isfinite(values).all(), name
    assert st["Nrw"][-1] > 0
