"""The cylinder of this fixture: the box [0,1]^2, cyclic in x and y, one
cell thick in z (front and back empty), less one circular cylinder (radius
R, centre C, no-slip); blockMesh's O-grid of four blocks, NT cells along
each side and NR from the side to the cylinder, graded to G at the
cylinder; ../stokes_case.py's solve.

  python3 cylinder_case.py dicts          write system/, constant/ and 0/
  python3 cylinder_case.py measure TIME   print the at-rest quantities at TIME (ascii U and Vc)
"""
import os
import sys

import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import stokes_case  # noqa: E402

NT, NR, G = 32, 16, 0.25
R, C, DZ = 0.25, (0.47, 0.53), 0.05


def block_mesh():
    """Box corners 0-3, the cylinder's points at 45 degrees 4-7, both planes."""
    cx, cy = C
    d = R / np.sqrt(2)
    xy = [(0, 0), (1, 0), (1, 1), (0, 1), (cx - d, cy - d), (cx + d, cy - d), (cx + d, cy + d), (cx - d, cy + d)]
    verts = "\n".join("    (%r %r %r)" % (x, y, z) for z in (0.0, DZ) for x, y in xy)
    mids = [(cx, cy - R), (cx + R, cy), (cx, cy + R), (cx - R, cy)]
    arcs = "\n".join("    arc %d %d (%r %r %r)" % (a + k, (a + 1) % 4 + k, m[0], m[1], z)
                     for k, z in ((4, 0.0), (12, DZ)) for a, m in zip(range(4), mids))
    blocks = "\n".join("    hex (%d %d %d %d %d %d %d %d) (%d %d 1) simpleGrading (1 %r 1)"
                       % (s, t, t + 4, s + 4, s + 8, t + 8, t + 12, s + 12, NT, NR, G)
                       for s, t in ((0, 1), (1, 2), (2, 3), (3, 0)))
    return """convertToMeters 1;
vertices
(
%s
);
blocks
(
%s
);
edges
(
%s
);
boundary
(
    left   { type cyclic; neighbourPatch right;  faces ((0 8 11 3)); }
    right  { type cyclic; neighbourPatch left;   faces ((1 2 10 9)); }
    bottom { type cyclic; neighbourPatch top;    faces ((0 1 9 8)); }
    top    { type cyclic; neighbourPatch bottom; faces ((3 11 10 2)); }
    cylinder { type wall; faces ((4 12 13 5) (5 13 14 6) (6 14 15 7) (7 15 12 4)); }
    frontAndBack
    {
        type empty;
        faces ((0 4 5 1) (1 5 6 2) (2 6 7 3) (3 7 4 0) (8 9 13 12) (9 10 14 13) (10 11 15 14) (11 8 12 15));
    }
);
""" % (verts, blocks, arcs)


def dicts():
    stokes_case.write("system/blockMeshDict", block_mesh())
    stokes_case.dicts("cylinder")


def measure(time):
    stokes_case.measure(time, 2, DZ)


if __name__ == "__main__":
    {"dicts": dicts, "measure": measure}[sys.argv[1]](*sys.argv[2:])
