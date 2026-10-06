"""Where the built apps are.

ctest passes PARTRAC_BIN; a bare `pytest tests/` falls back to build/bin, the
build tree the README describes, and then to bin/ in the source tree.
"""

import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _bin_dir():
    """The directory holding the app binaries."""
    env = os.environ.get("PARTRAC_BIN")
    if env:
        return env
    for candidate in (os.path.join(REPO, "build", "bin"), os.path.join(REPO, "bin")):
        if os.path.isdir(candidate):
            return candidate
    return os.path.join(REPO, "build", "bin")


BIN = _bin_dir()


def app(name):
    """The path of the app binary `name` (which may not exist)."""
    return os.path.join(BIN, name)

def built_with_dolfin():
    """Whether the binaries have mode=fenics, the one mode that needs dolfin.

    Having dolfin importable in python says nothing about how partrac was
    configured; CMake writes build_features.txt next to the apps with the
    options it was configured with. A build tree without that file is taken
    as built without dolfin until CMake runs again.
    """
    features = os.path.join(BIN, "build_features.txt")
    if not os.path.exists(features):
        return False
    with open(features) as f:
        return "dolfin=on" in f.read().split()
