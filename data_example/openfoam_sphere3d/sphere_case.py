"""The sphere of this fixture: the box [0,1]^3, cyclic in x, y and z, less
one sphere (radius R, centre C, no-slip), snappyHexMesh from an N^3
background, LEVEL levels of refinement at the sphere, no layers;
../stokes_case.py's solve.

  python3 sphere_case.py dicts          write system/, constant/ and 0/
  python3 sphere_case.py flatten        put the cyclic sides' points back on their planes
  python3 sphere_case.py measure TIME   print the at-rest quantities at TIME (ascii U and Vc)
"""
import os
import re
import sys

import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import stokes_case  # noqa: E402

N, LEVEL = 10, 2
R, C = 0.22, (0.47, 0.53, 0.44)
PLANES = {"left": (0, 0.0), "right": (0, 1.0), "front": (1, 0.0), "back": (1, 1.0),
          "bottom": (2, 0.0), "top": (2, 1.0)}

BLOCK_MESH = """convertToMeters 1;
vertices
(
    (0 0 0) (1 0 0) (1 1 0) (0 1 0)
    (0 0 1) (1 0 1) (1 1 1) (0 1 1)
);
blocks ( hex (0 1 2 3 4 5 6 7) (%(n)d %(n)d %(n)d) simpleGrading (1 1 1) );
boundary
(
    left   { type cyclic; neighbourPatch right;  faces ((0 4 7 3)); }
    right  { type cyclic; neighbourPatch left;   faces ((1 2 6 5)); }
    front  { type cyclic; neighbourPatch back;   faces ((0 1 5 4)); }
    back   { type cyclic; neighbourPatch front;  faces ((3 7 6 2)); }
    bottom { type cyclic; neighbourPatch top;    faces ((0 3 2 1)); }
    top    { type cyclic; neighbourPatch bottom; faces ((4 5 6 7)); }
);
"""

SNAPPY = """castellatedMesh true;
snap            true;
addLayers       false;

geometry
{
    sphere { type sphere; centre (%(c)s); radius %(r)r; }
}

castellatedMeshControls
{
    maxLocalCells 1000000;
    maxGlobalCells 2000000;
    minRefinementCells 0;
    maxLoadUnbalance 0.1;
    nCellsBetweenLevels 2;
    features ();
    refinementSurfaces
    {
        sphere { level (%(l)d %(l)d); patchInfo { type wall; } }
    }
    resolveFeatureAngle 30;
    refinementRegions {}
    insidePoint (0.0512 0.0573 0.0591);
    allowFreeStandingZoneFaces false;
}

snapControls
{
    nSmoothPatch 3;
    tolerance 2.0;
    nSolveIter 50;
    nRelaxIter 5;
}

addLayersControls
{
    relativeSizes true;
    layers {}
    expansionRatio 1.0;
    finalLayerThickness 0.3;
    minThickness 0.1;
}

meshQualityControls
{
    #includeEtc "caseDicts/mesh/generation/meshQualityDict.cfg"
}

writeFlags ();
mergeTolerance 1e-6;
"""

def dicts():
    stokes_case.write("system/blockMeshDict", BLOCK_MESH % dict(n=N))
    stokes_case.write("system/snappyHexMeshDict", SNAPPY % dict(c="%r %r %r" % C, r=R, l=LEVEL))
    stokes_case.dicts("sphere")


def flatten():
    """snappyHexMesh's snap moves the cyclic sides' points off their planes
    (they stay matched); every point of each side is put back on its plane."""
    P = "constant/polyMesh/"
    X = stokes_case.read_list(P + "points", 3)
    b = stokes_case.body(P + "faces")
    b = b[re.search(r"\d+\s*\(", b).end():]
    F = [np.array(m.split(), dtype=int) for m in re.findall(r"\d+\s*\(([^)]*)\)", b)]
    patches = {m.group(1): dict(re.findall(r"(\w+)\s+([^;]+);", m.group(2)))
               for m in re.finditer(r"(\w+)\s*\{([^}]*)\}", stokes_case.body(P + "boundary"))}
    for name, (axis, value) in PLANES.items():
        s, n = int(patches[name]["startFace"]), int(patches[name]["nFaces"])
        X[np.unique(np.concatenate(F[s:s + n])), axis] = value
    rows = "\n".join("(%.17g %.17g %.17g)" % tuple(p) for p in X)
    stokes_case.write(P + "points", "%d\n(\n%s\n)\n" % (len(X), rows), "vectorField")


def measure(time):
    stokes_case.measure(time, 3)


if __name__ == "__main__":
    {"dicts": dicts, "flatten": flatten, "measure": measure}[sys.argv[1]](*sys.argv[2:])
