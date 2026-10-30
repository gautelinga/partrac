"""Where the built apps are.

ctest passes PARTRAC_BIN; a bare `pytest tests/` falls back to build/bin, the
build tree the README describes, and then to bin/ in the source tree.
"""

import os

import pytest

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


def configured_apps():
    """The apps the build was configured with, from build_features.txt; None
    without one (another tree, or one CMake has not run in since)."""
    features = os.path.join(BIN, "build_features.txt")
    if not os.path.exists(features):
        return None
    with open(features) as f:
        for line in f.read().split():
            if line.startswith("apps="):
                return set(line[len("apps="):].split(","))
    return None


def left_out(binary):
    """Whether the build has not got the app at `binary` by choice: one it was
    not configured with, or, where nothing says, one not there. An app it was
    configured with and lacks is not left out: its tests fail."""
    configured = configured_apps()
    if configured is None:
        return not os.path.exists(binary)
    return os.path.basename(binary) not in configured


def needs(*binaries):
    """A mark skipping the test where the build has left out any of the apps."""
    out = [os.path.basename(b) for b in binaries if left_out(b)]
    return pytest.mark.skipif(bool(out), reason=" and ".join(out) + " not built")


def need(*binaries):
    """Skip here where the build has left out any of the apps."""
    out = [os.path.basename(b) for b in binaries if left_out(b)]
    if out:
        pytest.skip(" and ".join(out) + " not built")
