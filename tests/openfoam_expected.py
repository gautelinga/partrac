"""The reference split of an OpenFOAM fixture, split_expected.h5, which the
unit tests hold the C++ split to simplex for simplex. A port of the rules,
written apart from the C++ and sharing none of it:

  faces   each fanned from findFaceBasePts' base point: the first whose fan
          has tetrahedron quality above 1e-30 with the cell centres on both
          sides (across a cyclic, the owner side's choice), else the best
  W12     every cell the fan from its centre over its faces' triangles, a hex
          whose fan has a tet under 1e-12 of the cell taking W6's split
  W6      a hex Dompierre's split from the lowest (cyclic master, id) corner
          its three faces' diagonals meet at, 5 tets when no far diagonal
          passes the opposite corner, where every tet is above 1e-12 of the
          cell and they add up to the fan's volume; else the fan
  2D      the front plane of the empty pair: each face fanned from its cell
          centre (W12), or a quad cut along its base point's diagonal (W6)

and the patch of every simplex facet that is a boundary face's triangle (2D:
a side face's edge), looked up among all of them. The geometry is OpenFOAM's
formulas (face::areaAndCentre, primitiveMesh::makeCellCentresAndVols) on the
case's polyMesh, ascii or binary, plain or gzipped.

  python3 tests/openfoam_expected.py CASE [OUT.h5]   (default CASE/split_expected.h5)

writes w12/ and w6/ with cell_of, node_kind, node_point, node_cell,
facet_patch (int32) and the simplices: cells in full up to 20000 entries,
else cells_head (the first 4000) and the attributes cells_size and
cells_fnv1a (64-bit FNV-1a over their int64 little-endian bytes).
"""
import gzip
import os
import re
import sys

import h5py
import numpy as np

VALID = 1e-12
MIN_QUALITY = 1e-30
GREAT, ROOT_VSMALL = 1e15, 1e-150


# --- the polyMesh --------------------------------------------------------------

def read_foam(path):
    """A FoamFile's header dict and the bytes after it."""
    raw = open(path, "rb").read() if os.path.exists(path) else gzip.open(path + ".gz").read()
    h = raw.index(b"{", raw.index(b"FoamFile")) + 1
    e = raw.index(b"}", h)
    head = dict(re.findall(r"(\w+)\s+([^;]+);", raw[h:e].decode()))
    return head, raw[e + 1:]


def strip(body):
    return re.sub(rb"//[^\n]*|/\*.*?\*/", b" ", body, flags=re.S)


def labels(body, binary, pos=0):
    """The labelList after pos: its values and the position after it."""
    m = re.compile(rb"(\d+)\s*\(").search(body, pos)
    n = int(m.group(1))
    if binary:
        a = np.frombuffer(body[m.end():m.end() + 4 * n], dtype="<i4").astype(np.int64)
        return a, body.index(b")", m.end() + 4 * n) + 1
    e = body.index(b")", m.end())
    return np.array(body[m.end():e].split(), dtype=np.int64), e + 1


def read_labels(path):
    head, body = read_foam(path)
    binary = head["format"] == "binary"
    return labels(body if binary else strip(body), binary)[0]


def read_mesh(case):
    poly = os.path.join(case, "constant", "polyMesh")
    head, body = read_foam(os.path.join(poly, "points"))
    if head["format"] == "binary":
        m = re.compile(rb"(\d+)\s*\(").search(body)
        n = int(m.group(1))
        X = np.frombuffer(body[m.end():m.end() + 24 * n], dtype="<f8").reshape(n, 3).copy()
    else:
        body = strip(body)
        s = body.index(b"(")
        X = np.array(body[s + 1:body.rindex(b")")].replace(b"(", b" ").replace(b")", b" ").split(),
                     dtype=float).reshape(-1, 3)
    head, body = read_foam(os.path.join(poly, "faces"))
    if head["class"] == "faceCompactList":
        binary = head["format"] == "binary"
        body = body if binary else strip(body)
        start, p = labels(body, binary)
        pts, _ = labels(body, binary, p)
        faces = [[int(v) for v in pts[start[i]:start[i + 1]]] for i in range(len(start) - 1)]
    else:
        body = strip(body).decode()
        faces = [list(map(int, f.split())) for f in re.findall(r"\d+\(([^)]*)\)", body[body.index("(") + 1:])]
    own = read_labels(os.path.join(poly, "owner"))
    nbr = read_labels(os.path.join(poly, "neighbour"))
    body = strip(read_foam(os.path.join(poly, "boundary"))[1]).decode()
    patches = []
    for name, d in re.findall(r"(\w+)\s*\{([^}]*)\}", body):
        e = dict(re.findall(r"(\w+)\s+([^;]+);", d))
        patches.append({"name": name, "type": e["type"], "start": int(e["startFace"]),
                        "size": int(e["nFaces"]), "nbr": e.get("neighbourPatch")})
    for i, p in enumerate(patches):
        p["nbr"] = next((j for j, q in enumerate(patches) if q["name"] == p["nbr"]), None)
        p["owner"] = p["nbr"] is not None and i < p["nbr"]
    own, nbr = [int(v) for v in own], [int(v) for v in nbr]
    ncells = max(own + nbr) + 1
    return {"X": X, "faces": faces, "owner": own, "nbr": nbr, "ni": len(nbr), "ncells": ncells,
            "patches": patches}


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def sub(a, b):
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]


def geometry(m):
    """Face centres and areas, cell centres and volumes, cyclic separations, OpenFOAM's way."""
    X = [list(x) for x in m["X"]]
    fc, fa = [], []
    for F in m["faces"]:
        k = len(F)
        if k == 3:
            a, b, d = X[F[0]], X[F[1]], X[F[2]]
            n = cross(sub(b, a), sub(d, a))
            fa.append([0.5 * n[i] for i in range(3)])
            fc.append([(a[i] + b[i] + d[i]) / 3. for i in range(3)])
            continue
        pav = [0., 0., 0.]
        for j in range(k):
            for i in range(3):
                pav[i] += X[F[j]][i]
        pav = [v / float(k) for v in pav]
        A, sa = [], [0., 0., 0.]
        for j in range(k):
            x, y = X[F[j]], X[F[(j + 1) % k]]
            A.append(cross(sub(y, x), sub(pav, x)))
            for i in range(3):
                sa[i] += A[j][i]
        mag = dot(sa, sa) ** 0.5
        hat = [v / mag for v in sa] if mag > 0 else [0., 0., 0.]
        num, san = [0., 0., 0.], 0.
        for j in range(k):
            x, y = X[F[j]], X[F[(j + 1) % k]]
            an = dot(A[j], hat)
            san += an
            for i in range(3):
                num[i] += an * (x[i] + y[i] + pav[i])
        fc.append([num[i] / (3 * san) if san > 1e-300 else pav[i] for i in range(3)])
        fa.append([0.5 * sa[i] for i in range(3)])
    nc, ni = m["ncells"], m["ni"]
    sides = [(f, int(m["owner"][f]), 1.) for f in range(len(m["faces"]))] + \
            [(f, int(m["nbr"][f]), -1.) for f in range(ni)]
    est, cnt = [[0., 0., 0.] for _ in range(nc)], [0] * nc
    for f, c, _ in sides:
        for i in range(3):
            est[c][i] += fc[f][i]
        cnt[c] += 1
    est = [[v / float(cnt[c]) for v in est[c]] for c in range(nc)]
    ctr, vol = [[0., 0., 0.] for _ in range(nc)], [0.] * nc
    for f, c, sign in sides:
        pyr3 = sign * dot(fa[f], sub(fc[f], est[c]))
        for i in range(3):
            ctr[c][i] += pyr3 * (0.75 * fc[f][i] + 0.25 * est[c][i])
        vol[c] += pyr3
    m["C"] = [[ctr[c][i] / vol[c] for i in range(3)] if abs(vol[c]) > 1e-300 else est[c] for c in range(nc)]
    m["V"] = [v / 3. for v in vol]
    m["fc"], m["fa"] = fc, fa
    for p in m["patches"]:
        if p["type"] != "cyclic":
            continue
        q = m["patches"][p["nbr"]]
        sep = [0., 0., 0.]
        for i in range(p["size"]):
            for d in range(3):
                sep[d] += fc[q["start"] + i][d] - fc[p["start"] + i][d]
        p["sep"] = [s / float(p["size"]) if p["size"] else 0. for s in sep]


# --- the rules -----------------------------------------------------------------

def vol6(a, b, c, d):
    return dot(cross(sub(b, a), sub(c, a)), sub(d, a))


def quality(a, b, c, d):
    """tetrahedron::quality."""
    A, B, C = sub(b, a), sub(c, a), sub(d, a)
    ba, ca = cross(B, A), cross(C, A)
    lam, mu, den = dot(C, C) - dot(A, C), dot(B, B) - dot(A, B), dot(C, ba)
    if abs(den) < ROOT_VSMALL:
        r = 3 ** 0.5 * GREAT
    else:
        v = [(A[i] + (lam * ba[i] - mu * ca[i]) / den) / 2 for i in range(3)]
        r = dot(v, v) ** 0.5
    rr = min(r, GREAT)
    return (vol6(a, b, c, d) / 6.) / ((8.0 / 27.0) * 3 ** 0.5 * rr * rr * rr + ROOT_VSMALL)


def face_patch(m):
    fp = [-1] * len(m["faces"])
    for i, p in enumerate(m["patches"]):
        for f in range(p["start"], p["start"] + p["size"]):
            fp[f] = i
    return fp


def bases(m, only=None):
    """findFaceBasePts' base point of every face (of those in only)."""
    X = [list(x) for x in m["X"]]
    fp = face_patch(m)
    base = [-1] * len(m["faces"])

    def worst(cc, F, owner_side, b):
        k = len(F)
        q = 1e300
        for i in range(1, k - 1):
            a0 = (i + b) % k
            a1 = (a0 + 1) % k
            pa, pb = (F[a0], F[a1]) if owner_side else (F[a1], F[a0])
            q = min(q, quality(cc, X[F[b]], X[pa], X[pb]))
        return q

    for f, F in enumerate(m["faces"]):
        if only is not None and f not in only:
            continue
        own, nbr = m["C"][m["owner"][f]], None
        if f < m["ni"]:
            nbr = m["C"][m["nbr"][f]]
        elif m["patches"][fp[f]]["type"] == "cyclic":
            p = m["patches"][fp[f]]
            if not p["owner"]:
                continue
            q = m["patches"][p["nbr"]]
            nbr = sub(m["C"][m["owner"][f - p["start"] + q["start"]]], p["sep"])
        best_q, best = -1e300, 0
        for b in range(len(F)):
            q = worst(own, F, True, b)
            if nbr is not None:
                q = min(q, worst(nbr, F, False, b))
            if q > MIN_QUALITY:
                base[f] = b
                break
            if q > best_q:
                best_q, best = q, b
        if base[f] < 0:
            base[f] = best
    for p in m["patches"]:
        if p["type"] == "cyclic" and not p["owner"]:
            q = m["patches"][p["nbr"]]
            for i in range(p["size"]):
                b = base[q["start"] + i]
                if b >= 0:
                    base[p["start"] + i] = b if b < 1 else len(m["faces"][p["start"] + i]) - b
    return base


def triangles(F, b):
    j = 0 if len(F) == 3 else b
    k = len(F)
    return [(F[j], F[(j + i) % k], F[(j + i + 1) % k]) for i in range(1, k - 1)]


def masters(m):
    """The smallest cyclic image of every point."""
    comp = list(range(len(m["X"])))

    def root(i):
        while comp[i] != i:
            i = comp[i]
        return i

    for p in m["patches"]:
        if p["type"] == "cyclic" and p["owner"]:
            q = m["patches"][p["nbr"]]
            for i in range(p["size"]):
                fa, fb = m["faces"][p["start"] + i], m["faces"][q["start"] + i]
                k = len(fa)
                for j in range(k):
                    r1, r2 = root(fa[j]), root(fb[(k - j) % k])
                    comp[max(r1, r2)] = min(r1, r2)
    return [root(i) for i in range(len(comp))]


def split3d(m, n):
    X = [list(x) for x in m["X"]]
    npts, nc = len(X), m["ncells"]
    base = bases(m)
    tris = [triangles(F, base[f]) for f, F in enumerate(m["faces"])]
    key = masters(m)
    cell_faces = [[] for _ in range(nc)]
    for f in range(len(m["faces"])):
        cell_faces[m["owner"][f]].append(f)
        if f < m["ni"]:
            cell_faces[m["nbr"][f]].append(f)
    cell_faces = [sorted(fs) for fs in cell_faces]

    def x(v):
        return X[v] if v < npts else m["C"][v - npts]

    def vol(t):
        return vol6(x(t[0]), x(t[1]), x(t[2]), x(t[3])) / 6.

    def fan(c, apex):
        out = []
        for side in (0, 1):
            for f in cell_faces[c]:
                own = m["owner"][f] == c
                if (side == 0) != own:
                    continue
                out += [(apex, t[0], t[1], t[2]) if own else (apex, t[0], t[2], t[1]) for t in tris[f]]
        return out

    def dompierre(c):
        F = cell_faces[c]
        P = sorted({p for f in F for p in m["faces"][f]})
        on_face = lambda i, p: p in m["faces"][F[i]]
        on_diag = lambda i, p: p in (tris[F[i]][0][0], tris[F[i]][0][2])
        cand = [p for p in P if all(on_diag(i, p) for i in range(6) if on_face(i, p))]
        if not cand:
            return None
        v0 = min(cand, key=lambda p: (key[p], p))
        v6 = [p for p in P if not any(on_face(i, v0) and on_face(i, p) for i in range(6))][-1]
        cone, through = [], False
        for i in range(6):
            if on_face(i, v0):
                continue
            through = through or on_diag(i, v6)
            f = F[i]
            own = m["owner"][f] == c
            cone += [(v0, t[0], t[1], t[2]) if own else (v0, t[0], t[2], t[1]) for t in tris[f][:2]]
        if through:
            return cone, False
        keep = [t for t in cone if v6 not in t[1:]]
        u = sorted({v for t in cone if v6 in t[1:] for v in t[1:] if v != v6})
        if vol6(X[v0], X[u[0]], X[u[1]], X[u[2]]) < 0:
            u[1], u[2] = u[2], u[1]
        return keep + [(v0, u[0], u[1], u[2]), (v6, u[0], u[2], u[1])], True

    hexes = {}
    for c in range(nc):
        F = cell_faces[c]
        hexa = len(F) == 6 and all(len(m["faces"][f]) == 4 for f in F) \
            and len({p for f in F for p in m["faces"][f]}) == 8
        if not hexa:
            continue
        vc = m["V"][c]
        if n == 12 and all(vol(t) > VALID * vc for t in fan(c, npts + c)):
            continue
        r = dompierre(c)
        if r is None:
            continue
        v = [vol(t) for t in r[0]]
        fv = sum(vol(t) for t in fan(c, npts + c))
        if all(a > VALID * vc for a in v) and abs(sum(v) - fv) <= 1e-9 * abs(fv):
            hexes[c] = r[0]
    fan_ids = [c for c in range(nc) if c not in hexes]
    apex = {c: npts + i for i, c in enumerate(fan_ids)}
    cells, cell_of = [], []
    for c in range(nc):
        t = hexes[c] if c in hexes else fan(c, apex[c])
        cells += t
        cell_of += [c] * len(t)
    fp = face_patch(m)
    boundary = {tuple(sorted(t)): fp[f] for f in range(m["ni"], len(m["faces"])) for t in tris[f]}
    facet = []
    for t in cells:
        for k in range(4):
            q = tuple(sorted(t[j] for j in range(4) if j != k))
            facet.append(boundary.get(q, -1) if max(q) < npts else -1)
    return {"cells": cells, "cell_of": cell_of, "node_kind": [0] * npts + [1] * len(fan_ids),
            "node_point": list(range(npts)) + [-1] * len(fan_ids),
            "node_cell": [-1] * npts + fan_ids, "facet_patch": facet}


def split2d(m, n):
    X = m["X"]
    fp = face_patch(m)
    empty = [f for f in range(m["ni"], len(m["faces"])) if m["patches"][fp[f]]["type"] == "empty"]
    a = m["fa"][empty[0]]
    axis = max(range(3), key=lambda d: (abs(a[d]), -d))
    inplane = [d for d in range(3) if d != axis]
    size = float(np.sqrt(((X.max(0) - X.min(0)) ** 2).sum()))
    lo = min(m["fc"][f][axis] for f in empty)
    front = {}
    for f in empty:
        if abs(m["fc"][f][axis] - lo) <= 1e-9 * size:
            front[m["owner"][f]] = f
    base = bases(m, set(front.values()))
    in2 = lambda p: (X[p][inplane[0]], X[p][inplane[1]])
    area2 = lambda a, b, c: (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    split, fanned = {}, []
    for c in range(m["ncells"]):
        f = front[c]
        poly = list(m["faces"][f])
        k = len(poly)
        ar = sum(in2(poly[i])[0] * in2(poly[(i + 1) % k])[1] - in2(poly[(i + 1) % k])[0] * in2(poly[i])[1]
                 for i in range(k))
        j = base[f]
        if ar < 0:
            poly.reverse()
            j = k - 1 - j
        if k == 3:
            split[c] = [tuple(poly)]
            continue
        if k == 4 and n == 6:
            t1 = (poly[j], poly[(j + 1) % 4], poly[(j + 2) % 4])
            t2 = (poly[j], poly[(j + 2) % 4], poly[(j + 3) % 4])
            if all(area2(*[in2(p) for p in t]) > 2 * VALID * (abs(ar) / 2) for t in (t1, t2)):
                split[c] = [t1, t2]
                continue
        fanned.append(c)
        split[c] = poly
    pts = sorted({p for f in front.values() for p in m["faces"][f]})
    node = {p: i for i, p in enumerate(pts)}
    apex = {c: len(pts) + i for i, c in enumerate(fanned)}
    cells, cell_of = [], []
    for c in range(m["ncells"]):
        if c in apex:
            poly = split[c]
            t = [(apex[c], node[poly[i]], node[poly[(i + 1) % len(poly)]]) for i in range(len(poly))]
        else:
            t = [tuple(node[p] for p in q) for q in split[c]]
        cells += t
        cell_of += [c] * len(t)
    node_point = pts + [-1] * len(fanned)
    boundary = {}
    for f in range(m["ni"], len(m["faces"])):
        if m["patches"][fp[f]]["type"] != "empty":
            boundary[tuple(sorted(p for p in m["faces"][f] if p in node))] = fp[f]
    facet = []
    for t in cells:
        for k in range(3):
            q = [node_point[t[j]] for j in range(3) if j != k]
            facet.append(boundary.get(tuple(sorted(q)), -1) if min(q) >= 0 else -1)
    return {"cells": cells, "cell_of": cell_of, "node_kind": [0] * len(pts) + [1] * len(fanned),
            "node_point": node_point, "node_cell": [-1] * len(pts) + fanned, "facet_patch": facet}


def fnv1a(a):
    h = 0xcbf29ce484222325
    for b in np.asarray(a, dtype="<i8").tobytes():
        h = ((h ^ b) * 0x100000001b3) & 0xFFFFFFFFFFFFFFFF
    return h


def main(case, out=None, full_max=20000, head=4000):
    m = read_mesh(case)
    geometry(m)
    two = any(p["type"] == "empty" and p["size"] > 0 for p in m["patches"])
    with h5py.File(out or os.path.join(case, "split_expected.h5"), "w") as h:
        h.attrs["source"] = "tests/openfoam_expected.py"
        for n in (12, 6):
            s = split2d(m, n) if two else split3d(m, n)
            g = h.create_group("w%d" % n)
            cells = np.asarray(s["cells"], dtype=np.int64).ravel()
            arrays = [(k, s[k]) for k in ("cell_of", "node_kind", "node_point", "node_cell", "facet_patch")]
            if len(cells) <= full_max:
                arrays.append(("cells", cells))
            else:
                arrays.append(("cells_head", cells[:head]))
                g.attrs["cells_size"] = len(cells)
                g.attrs["cells_fnv1a"] = np.uint64(fnv1a(cells))
            for name, a in arrays:
                g.create_dataset(name, data=np.asarray(a, dtype=np.int32).ravel(),
                                 compression="gzip", compression_opts=9, shuffle=True)
            print(case, "w%d:" % n, len(s["cell_of"]), "simplices,", len(s["node_kind"]), "nodes")


if __name__ == "__main__":
    main(*sys.argv[1:])
