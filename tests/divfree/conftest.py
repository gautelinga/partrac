"""What the tests of python/divfree share across the workers of a session."""

import os

import pytest

from cases import shared_dir


@pytest.fixture(scope="session")
def channel_runs(tmp_path_factory):
    """channel_runs(dim) -> a folder holding divfree_cases.two_stamp_channel as
    a case (in/), cleaned on one rank in this process (one/, one.log) and by a
    job of three ranks on the route the writer takes by itself (three/,
    three.log). Made once a session for the tests that hold a run on three
    ranks against one rank's; a test copies what it writes into. Skipped where
    no job of several ranks runs, before anything is made."""
    import divfree_cases as C

    def get(dim):
        C.launcher_or_skip()

        def build(d):
            X, cells, per, fields, phi = C.two_stamp_channel(dim)
            cfg = C.write_case(d / "in", X, cells, fields, per, stamps=["0", "1"], phi=phi)
            (d / "one.log").write_text(C.one_rank_main([cfg, "--out", d / "one"]))
            env = dict(os.environ)
            env.pop("PARTRAC_HDF5_GATHER", None)
            r = C.mpi_run(3, [C.CLEANER, cfg, "--out", d / "three"], timeout=300, env=env)
            assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
            (d / "three.log").write_text(r.stdout)

        return shared_dir(tmp_path_factory, "channel_%dd" % dim, build)

    return get
