"""The reference values of this fixture: OpenFOAM's own cellPoint U at the
points `interpol ... Nrw=N random=false seed=1 num_threads=1` probes
(its mt19937 per thread, seeded by seed_seq{seed, thread}, three draws a
point, right to left, in the box of the mesh points; z = 0 in 2D, where the
probes sit on the slab's mid-plane).

  python3 probe_points.py write N [2d]   system/probesCP at those points
  python3 probe_points.py collect        postProcessing/probesCP -> cellpoint_U.txt
"""
import glob
import gzip
import os
import re
import sys

import numpy as np

M32 = 0xFFFFFFFF


def seed_seq(values, n=624):
    """std::seed_seq(values).generate over n words."""
    b = [0x8b8b8b8b] * n
    s = len(values)
    t = 11
    p = (n - t) // 2
    q = p + t
    m = max(s + 1, n)
    T = lambda x: x ^ (x >> 27)
    for k in range(m):
        r1 = (1664525 * T(b[k % n] ^ b[(k + p) % n] ^ b[(k - 1) % n])) & M32
        r2 = (r1 + (s if k == 0 else (k % n + values[k - 1]) if k <= s else k % n)) & M32
        b[(k + p) % n] = (b[(k + p) % n] + r1) & M32
        b[(k + q) % n] = (b[(k + q) % n] + r2) & M32
        b[k % n] = r2
    for k in range(m, m + n):
        r3 = (1566083941 * T((b[k % n] + b[(k + p) % n] + b[(k - 1) % n]) & M32)) & M32
        r4 = (r3 - k % n) & M32
        b[(k + p) % n] ^= r3
        b[(k + q) % n] ^= r4
        b[k % n] = r4
    return b


def interpol_points(lo, hi, n, seed=1):
    key = np.array(seed_seq([seed & M32, 0]), dtype=np.uint32)
    g = np.random.MT19937()
    g.state = {"bit_generator": "MT19937", "state": {"key": key, "pos": 624}}
    w = g.random_raw(6 * n).astype(np.float64).reshape(3 * n, 2)
    r = (w[:, 0] + w[:, 1] * 4294967296.0) / 18446744073709551616.0
    r[r >= 1] = np.nextafter(1.0, 0.0)
    r = r.reshape(n, 3)
    X = np.empty((n, 3))
    for j, a in enumerate((2, 1, 0)):
        X[:, a] = r[:, j] * (hi[a] - lo[a]) + lo[a]
    return X


def read_points():
    path = "constant/polyMesh/points"
    text = gzip.open(path + ".gz", "rt").read() if os.path.exists(path + ".gz") else open(path).read()
    body = text[text.index("}", text.index("FoamFile")) + 1:]
    body = re.sub(r"//[^\n]*|/\*.*?\*/", " ", body, flags=re.S)
    v = np.array(body.replace("(", " ").replace(")", " ").split(), dtype=float)
    return v[1:].reshape(int(v[0]), 3)


def write(n, two_d):
    X = read_points()
    lo, hi = X.min(0), X.max(0)
    zmid = 0.5 * (lo[2] + hi[2])
    if two_d:
        lo[2] = hi[2] = 0.0
    P = interpol_points(lo, hi, n)
    np.savetxt("probe_points.txt", P, fmt="%.17g")
    if two_d:
        P[:, 2] = zmid
    pts = "\n".join("        (%.17g %.17g %.17g)" % tuple(p) for p in P)
    with open("system/probesCP", "w") as f:
        f.write("""FoamFile
{
    format      ascii;
    class       dictionary;
    location    "system";
    object      probesCP;
}

type            probes;
libs            ("libsampling.so");
writeControl    timeStep;
writeInterval   1;
fields          (U);
fixedLocations  true;
interpolationScheme cellPoint;
probeLocations
    (
%s
    );
""" % pts)


def collect():
    [out] = glob.glob("postProcessing/probesCP/*/U")
    lines = open(out).read().splitlines()
    data = [l for l in lines if l.strip() and not l.startswith("#")][-1]
    U = np.array(data.split(None, 1)[1].replace("(", " ").replace(")", " ").split(), dtype=float).reshape(-1, 3)
    P = np.loadtxt("probe_points.txt")
    np.savetxt("cellpoint_U.txt", np.c_[P, U], fmt="%.17g",
               header="x y z of interpol's points (Nrw=%d random=false seed=1 num_threads=1), "
                      "OpenFOAM's cellPoint U there; -1e300 outside the mesh" % len(P))
    os.remove("probe_points.txt")


if __name__ == "__main__":
    if sys.argv[1] == "write":
        write(int(sys.argv[2]), len(sys.argv) > 3 and sys.argv[3] == "2d")
    else:
        collect()
