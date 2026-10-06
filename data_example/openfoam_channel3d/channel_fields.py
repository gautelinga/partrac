"""The jittered channel of this fixture: [0,1]^3, 12^3 hexes, cyclic in x
and y, walls at z = 0 and 1; interior points moved by up to 0.2 h a
component (seed 1), cyclic-side points along their side only, both images
alike, wall points not at all. The field is linear in z, periodic in x and
y: u = z (2, 0.7, -1.1) + (1, -0.5, 0.3), at the cell centres OpenFOAM
computes (0/C) and on the walls, p = 0.

  python3 channel_fields.py jitter    move constant/polyMesh/points
  python3 channel_fields.py fields    write 0/U and 0/p from 0/C
"""
import re
import sys

import numpy as np

N, FRAC, SEED = 12, 0.2, 1
AZ = np.array([2.0, 0.7, -1.1])
BZ = np.array([1.0, -0.5, 0.3])
HEADER = """FoamFile
{
    format      ascii;
    class       %s;
    location    "%s";
    object      %s;
}
"""


def read_vectors(path):
    text = open(path).read()
    body = text[text.index("}", text.index("FoamFile")) + 1:]
    body = re.sub(r"//[^\n]*|/\*.*?\*/", " ", body, flags=re.S)
    m = re.search(r"(\d+)\s*\(", body)
    n = int(m.group(1))
    close = body.index(")\n)", m.end()) + 1 if ")\n)" in body else body.rindex(")")
    v = np.array(body[m.end():close].replace("(", " ").replace(")", " ").split(), dtype=float)
    return v[:3 * n].reshape(n, 3)


def vectors(X):
    return "%d\n(\n%s\n)\n" % (len(X), "\n".join("(%.17g %.17g %.17g)" % tuple(x) for x in X))


def jitter():
    path = "constant/polyMesh/points"
    X = read_vectors(path)
    ijk = np.rint(X * N).astype(int)
    assert np.abs(ijk / N - X).max() < 1e-9, "not a uniform lattice"
    rng = np.random.default_rng(SEED)
    D = rng.uniform(-FRAC / N, FRAC / N, (N, N, N + 1, 3))
    d = D[ijk[:, 0] % N, ijk[:, 1] % N, ijk[:, 2]].copy()
    d[(ijk[:, 0] == 0) | (ijk[:, 0] == N), 0] = 0.
    d[(ijk[:, 1] == 0) | (ijk[:, 1] == N), 1] = 0.
    d[(ijk[:, 2] == 0) | (ijk[:, 2] == N)] = 0.
    with open(path, "w") as f:
        f.write(HEADER % ("vectorField", "constant/polyMesh", "points") + "\n" + vectors(X + d))


def fields():
    C = read_vectors("0/C")
    u = C[:, 2][:, None] * AZ + BZ
    walls = {"bottom": BZ, "top": AZ + BZ}
    with open("0/U", "w") as f:
        f.write(HEADER % ("volVectorField", "0", "U") + "\ndimensions      [0 1 -1 0 0 0 0];\n"
                "internalField   nonuniform List<vector> " + vectors(u) + ";\nboundaryField\n{\n"
                + "".join("    %s { type cyclic; }\n" % s for s in ("left", "right", "front", "back"))
                + "".join("    %s { type fixedValue; value uniform (%.17g %.17g %.17g); }\n" % ((w,) + tuple(v))
                          for w, v in walls.items()) + "}\n")
    with open("0/p", "w") as f:
        f.write(HEADER % ("volScalarField", "0", "p") + "\ndimensions      [0 2 -2 0 0 0 0];\n"
                "internalField   uniform 0;\nboundaryField\n{\n"
                + "".join("    %s { type cyclic; }\n" % s for s in ("left", "right", "front", "back"))
                + "    bottom { type zeroGradient; }\n    top { type zeroGradient; }\n}\n")


if __name__ == "__main__":
    {"jitter": jitter, "fields": fields}[sys.argv[1]]()
