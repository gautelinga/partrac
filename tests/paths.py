"""Where the built apps are.

ctest passes PARTRAC_BIN; a bare `pytest tests/` falls back to the build tree
the README tells you to make, then to the old in-source location.
"""

import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _bin_dir():
    env = os.environ.get("PARTRAC_BIN")
    if env:
        return env
    for candidate in (os.path.join(REPO, "build", "bin"), os.path.join(REPO, "bin")):
        if os.path.isdir(candidate):
            return candidate
    return os.path.join(REPO, "build", "bin")


BIN = _bin_dir()


def app(name):
    return os.path.join(BIN, name)
