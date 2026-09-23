"""Running an app on a private copy of an example.

Every app writes its output beside its parameter file, so a test runs on a
copy. The apps refuse a parameter given twice, and a test states a common set
of arguments and overrides a few of them, so the arguments are merged by key
here, the later value winning.
"""

import os
import shutil
import subprocess


def merged(*argsets):
    """The arguments of every set, one per key, a later set's value winning.

    A set is a string of whitespace-separated key=value arguments, or a list
    of such strings and sets. A word with no = after a path continues that
    path (a restart folder may hold a space); anywhere else it stays an
    argument of its own, which the app then refuses.
    """
    argv = {}
    for args in argsets:
        words = []
        for w in args.split() if isinstance(args, str) else merged(*args):
            if words and "/" in words[-1] and "=" not in w and not w.startswith("--"):
                words[-1] += " " + w
            else:
                words.append(w)
        argv.update((w.split("=")[0], w) for w in words)
    return list(argv.values())


def copy_case(src, d):
    """Copy the input folder src, less the mesh generator script, to the new folder d; return d."""
    shutil.copytree(src, d, ignore=shutil.ignore_patterns("generate_up.py"))
    return d


def copy_example(example, d):
    """Copy the parameter file `example` into the folder d (made if missing); return the copy."""
    d.mkdir(parents=True, exist_ok=True)
    shutil.copy(example, d / os.path.basename(example))
    return d / os.path.basename(example)


def run_app(binary, params, *argsets, env=None, check=True, timeout=900):
    """Run binary on the parameter file params with the merged arguments; return the process.

    env adds to the environment. With check, the run must exit 0.
    """
    r = subprocess.run([str(binary), str(params)] + merged(*argsets),
                       capture_output=True, text=True, timeout=timeout,
                       env=dict(os.environ, **(env or {})))
    if check:
        assert r.returncode == 0, r.stdout + r.stderr
    return r


def checkpoint_folder(d):
    """The folder to resume from: the parent of the one Checkpoints folder under d."""
    found = list(d.rglob("Checkpoints/positions.pos"))
    assert len(found) == 1, found
    return found[0].parent.parent


def continuous_and_resumed(binary, example, root, args, stop, end, env=None):
    """(cont, split) folders: a run to `end`, and one to `stop` resumed from its checkpoint to `end`.

    args, stop and end are argument sets; stop and end override args, and end
    serves both the uninterrupted run and the resumed one. The final
    checkpoint is written one step past the stop, and the resumed run writes
    into the tree of the run it resumes.
    """
    cont = copy_example(example, root / "cont")
    split = copy_example(example, root / "split")
    run_app(binary, cont, args, end, env=env)
    run_app(binary, split, args, stop, env=env)
    run_app(binary, split, args, end, ["restart_folder=%s" % checkpoint_folder(split.parent)], env=env)
    return cont.parent, split.parent
