"""OpenFOAM cases written with numpy, for the reader's tests: a lattice of
hexes, graded along each axis, its interior points jittered, cells removed
as obstacles; each pair of sides cyclic, or each side a wall, a symmetry
plane, a patch or empty (one cell thick: 2D). The fields are written at
OpenFOAM's own cell centres (openfoam_expected.centres), so a closed form
evaluated there is the field a case holding it would have. Everything ascii, as OpenFOAM
writes it.

  m = lattice_case(case, xs, ys, zs, sides, solid=None, jitter=0.)
  write_field(case, "0", "U", m, values, bcs)
  write_mesh(case, X, faces, owner, neighbour, patches)   a polyMesh given whole
  write_boundary(case, patches)                          its boundary file alone

sides gives each axis "cyclic", "empty" (both sides one patch), or a pair
of patch types for its low and high side. The sides are named left/right
(x), front/back (y), bottom/top (z), frontAndBack when empty; removed
cells' faces form the wall patch obstacles.
"""
import os

import numpy as np

from openfoam_expected import centres

NAMES = (("left", "right"), ("front", "back"), ("bottom", "top"))

HEADER = """FoamFile
{
    format      ascii;
    class       %s;
    location    "%s";
    object      %s;
}
"""

SYSTEM = {
    "controlDict": """application     foamRun;
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         0;
deltaT          1;
writeControl    timeStep;
writeInterval   1;
writeFormat     ascii;
writePrecision  17;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable false;
""",
    "fvSchemes": """ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""",
    "fvSolution": "solvers {}\n",
}


def num(x):
    return repr(float(x))


def vec(v):
    return "(%s)" % " ".join(num(x) for x in v)


def write(path, cls, location, obj, body):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(HEADER % (cls, location, obj) + "\n" + body)


def quad(a, s, p, q):
    """The lattice points of the face on plane s of axis a at (p, q) along the
    other two axes, in the order whose normal is +a."""
    b, c = (a + 1) % 3, (a + 2) % 3
    pts = []
    for dp, dq in ((0, 0), (1, 0), (1, 1), (0, 1)):
        ijk = [0, 0, 0]
        ijk[a], ijk[b], ijk[c] = s, p + dp, q + dq
        pts.append(tuple(ijk))
    return pts


def reverse(f):
    """The face the other way round, point 0 kept, as OpenFOAM pairs a cyclic's faces."""
    return [f[0], f[3], f[2], f[1]]


def lattice_case(case, xs, ys, zs, sides, solid=None, jitter=0., seed=1):
    """The case's system/ and constant/polyMesh; returns what write_field and
    the tests read: the cell centres C, the face centres fc and the patches
    (name, type, start, size, extra)."""
    lines = [np.asarray(v, dtype=float) for v in (xs, ys, zs)]
    n = [len(v) - 1 for v in lines]
    solid = np.zeros(n, dtype=bool) if solid is None else np.asarray(solid, dtype=bool)
    P = np.stack(np.meshgrid(*lines, indexing="ij"), -1)
    if jitter:
        # interior points by up to jitter of their smaller neighbouring spacing
        rng = np.random.default_rng(seed)
        inner = P[1:-1, 1:-1, 1:-1]
        for a in range(3):
            h = np.minimum(np.diff(lines[a])[:-1], np.diff(lines[a])[1:])
            shape = [1, 1, 1]
            shape[a] = len(h)
            inner[..., a] += jitter * h.reshape(shape) * rng.uniform(-1, 1, inner.shape[:3])
    ijk = np.argwhere(~solid)
    cell = -np.ones(n, dtype=np.int64)
    cell[tuple(ijk.T)] = np.arange(len(ijk))

    internal, bound = [], {}
    def side(name, kind, extra=None):
        bound.setdefault(name, {"type": kind, "faces": [], "extra": extra or {}})
        return bound[name]["faces"]
    obstacles = []
    for a in range(3):
        b, c = (a + 1) % 3, (a + 2) % 3
        kind = sides["xyz"[a]]
        lo, hi = NAMES[a]
        if kind == "empty":
            lo = hi = "frontAndBack"
        L = lines[a][-1] - lines[a][0]
        for p in range(n[b]):
            for q in range(n[c]):
                for s in range(n[a] + 1):
                    def at(i):
                        ix = [0, 0, 0]
                        ix[a], ix[b], ix[c] = i, p, q
                        return cell[tuple(ix)]
                    f = quad(a, s, p, q)
                    left = at(s - 1) if s > 0 else -1
                    right = at(s) if s < n[a] else -1
                    if 0 < s < n[a]:
                        if left >= 0 and right >= 0:
                            internal.append((left, right, f))
                        elif left >= 0:
                            obstacles.append((left, f))
                        elif right >= 0:
                            obstacles.append((right, reverse(f)))
                        continue
                    if kind == "cyclic":
                        across = at(n[a] - 1) if s == 0 else at(0)
                        mine = right if s == 0 else left
                        if mine < 0:
                            continue
                        if across < 0:
                            obstacles.append((mine, reverse(f) if s == 0 else f))
                            continue
                        sep = [0., 0., 0.]
                        sep[a] = -L if s == 0 else L
                        name, other = (lo, hi) if s == 0 else (hi, lo)
                        side(name, "cyclic", {"neighbourPatch": other, "transformType": "translational",
                                              "matchTolerance": "0.0001", "separation": vec(sep)}) \
                            .append((mine, reverse(f) if s == 0 else f))
                        continue
                    kinds = ("empty", "empty") if kind == "empty" else kind
                    if s == 0 and right >= 0:
                        side(lo, kinds[0]).append((right, reverse(f)))
                    elif s == n[a] and left >= 0:
                        side(hi, kinds[1]).append((left, f))
    if obstacles:
        bound["obstacles"] = {"type": "wall", "faces": obstacles, "extra": {}}

    internal.sort(key=lambda t: (t[0], t[1]))
    faces = [f for _, _, f in internal]
    owner = [o for o, _, _ in internal]
    neighbour = [nb for _, nb, _ in internal]
    patches = []
    for name, pt in bound.items():
        patches.append({"name": name, "type": pt["type"], "start": len(faces), "size": len(pt["faces"]),
                        "extra": pt["extra"]})
        faces += [f for _, f in pt["faces"]]
        owner += [o for o, _ in pt["faces"]]

    # the points the faces use, numbered in lattice order
    shape = np.array(n) + 1
    flat = np.array([[np.ravel_multi_index(v, shape) for v in f] for f in faces])
    used = np.unique(flat)
    renum = -np.ones(np.prod(shape), dtype=np.int64)
    renum[used] = np.arange(len(used))
    F = renum[flat]
    X = P.reshape(-1, 3)[used]
    owner, neighbour = np.array(owner), np.array(neighbour, dtype=np.int64)
    fc, _, C, _ = centres(X, F, owner, neighbour, len(ijk))
    write_mesh(case, X, F, owner, neighbour, patches)
    for name, text in SYSTEM.items():
        write(os.path.join(case, "system", name), "dictionary", "system", name, text)
    return {"C": C, "fc": fc, "patches": patches}


def write_mesh(case, X, faces, owner, neighbour, patches):
    """constant/polyMesh: the points X, the faces (lists of point labels),
    owner and neighbour, and the patches, each a dict of name, type, start,
    size and, optionally, extra entries."""
    poly = os.path.join(case, "constant", "polyMesh")
    write(os.path.join(poly, "points"), "vectorField", "constant/polyMesh", "points",
          "%d\n(\n%s\n)\n" % (len(X), "\n".join(vec(x) for x in X)))
    write(os.path.join(poly, "faces"), "faceList", "constant/polyMesh", "faces",
          "%d\n(\n%s\n)\n" % (len(faces), "\n".join("%d(%s)" % (len(f), " ".join(map(str, f))) for f in faces)))
    for name, v in (("owner", owner), ("neighbour", neighbour)):
        write(os.path.join(poly, name), "labelList", "constant/polyMesh", name,
              "%d\n(\n%s\n)\n" % (len(v), "\n".join(map(str, v))))
    write_boundary(case, patches)


def write_boundary(case, patches):
    """constant/polyMesh/boundary of the patches, as write_mesh takes them."""
    poly = os.path.join(case, "constant", "polyMesh")
    body = "%d\n(\n" % len(patches)
    for p in patches:
        body += "    %s\n    {\n        type %s;\n        nFaces %d;\n        startFace %d;\n" \
                % (p["name"], p["type"], p["size"], p["start"])
        body += "".join("        %s %s;\n" % kv for kv in p.get("extra", {}).items())
        body += "    }\n"
    write(os.path.join(poly, "boundary"), "polyBoundaryMesh", "constant/polyMesh", "boundary", body + ")\n")


def write_field(case, time, name, m, values, bcs=None):
    """A vector (n x 3) or scalar field at the cell centres. bcs by patch name:
    a condition's type, or (type, value) with value a vector, a scalar, or an
    array of one a face, or (type, value, entries) with entries a dict of the
    condition's other entries as written, value None for none; unnamed
    patches take their type's default (walls noSlip or zeroGradient, plain
    patches zeroGradient, coupled and constrained sides their own)."""
    values = np.asarray(values, dtype=float)
    vector = values.ndim == 2
    fmt = vec if vector else num
    body = "dimensions [0 1 -1 0 0 0 0];\n\n" if vector else "dimensions [0 2 -2 0 0 0 0];\n\n"
    body += "internalField nonuniform List<%s>\n%d\n(\n%s\n);\n\nboundaryField\n{\n" \
            % ("vector" if vector else "scalar", len(values), "\n".join(fmt(v) for v in values))
    for p in m["patches"]:
        bc = (bcs or {}).get(p["name"])
        if bc is None:
            bc = {"wall": "noSlip" if vector else "zeroGradient", "patch": "zeroGradient"}.get(p["type"], p["type"])
        body += "    %s\n    {\n" % p["name"]
        if isinstance(bc, str):
            body += "        type %s;\n" % bc
        else:
            kind, value, entries = bc if len(bc) == 3 else (bc[0], bc[1], {})
            body += "        type %s;\n" % kind
            body += "".join("        %s %s;\n" % kv for kv in entries.items())
            if value is not None:
                value = np.asarray(value, dtype=float)
                if value.ndim == (2 if vector else 1):
                    text = "nonuniform List<%s> %d(%s)" % ("vector" if vector else "scalar", len(value),
                                                           " ".join(fmt(v) for v in value))
                else:
                    text = "uniform " + fmt(value)
                body += "        value %s;\n" % text
        body += "    }\n"
    write(os.path.join(case, time, name), "volVectorField" if vector else "volScalarField", time, name,
          body + "}\n")
