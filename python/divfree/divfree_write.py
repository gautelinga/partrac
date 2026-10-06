"""The cleaner's output files written by every rank from its own rows.

The layout is `divfree_clean`'s serial writers'; what changes is who writes
it. Each rank passes its contiguous block of the rows in serial order -- cells
`[c0, c1)`, nodes `[n0, n1)` -- the bounds an exclusive scan of the counts over
the ranks; a rank may pass none. Every call of `BlockWriter` is collective and
returns once the file is complete.

Two routes, the same datasets from the same code: through h5py's MPI driver,
each rank writing its slab in pieces of at most PIECE_BYTES; or gathered to
the first rank, which then holds the whole field and writes it through the
serial driver. The first is taken where h5py has MPI and there is more than
one rank, unless PARTRAC_HDF5_GATHER is 1 or true. A failure on one rank
would leave the others in a collective, so every write runs under `guard`,
which aborts the job naming the rank and what it could not write; on one rank
the exception is raised.
"""

import contextlib
import os
import sys

import h5py
import numpy as np
from mpi4py import MPI

PIECE_BYTES = 1 << 30     # the largest single write or message
GATHER_ENV = "PARTRAC_HDF5_GATHER"
_TAG = 7301


def _clean():
    """divfree_clean, imported when first needed so that it can import this
    module at its top."""
    import divfree_clean
    return divfree_clean


def gather_forced():
    """Whether PARTRAC_HDF5_GATHER asks for the gathered route: 1 or true."""
    return os.environ.get(GATHER_ENV, "").strip().lower() in ("1", "true")


def parallel_hdf5(comm):
    """Whether the parallel route is available: h5py built with MPI, more than
    one rank, and the gather not forced. Collective where there is more than one
    rank, so every rank takes the same route."""
    if comm.size == 1:
        return False
    ok = bool(h5py.get_config().mpi) and not gather_forced()
    return bool(comm.allreduce(ok, op=MPI.LAND))


def _pieces(n, rowbytes):
    """[a, b) ranges of n rows, each at most PIECE_BYTES."""
    step = max(1, PIECE_BYTES // max(int(rowbytes), 1))
    return [(a, min(a + step, n)) for a in range(0, n, step)]


def _nloc(degree, dim):
    """The element's nodes in a cell."""
    return dim + 1 if degree == 1 else (dim + 1) * (dim + 2) // 2


def _create(g, name, n, dtype):
    """A contiguous one-dimensional dataset, never filled: every entry is
    written."""
    dcpl = h5py.h5p.create(h5py.h5p.DATASET_CREATE)
    dcpl.set_fill_time(h5py.h5d.FILL_TIME_NEVER)
    return g.create_dataset(name, (n,), dtype=dtype, dcpl=dcpl)


def _put(dset, start, rows):
    """rows written into dset from row start on, a piece at a time."""
    rowbytes = rows.itemsize * int(np.prod(rows.shape[1:], dtype=np.int64))
    for a, b in _pieces(len(rows), rowbytes):
        dset[start + a:start + b] = rows[a:b]


def _check_path(path, mode):
    """Open the file as the OS would before the collective open does: a rank
    that cannot fails here, alone, instead of inside MPI_File_open."""
    if mode == "a" and not os.path.exists(str(path)):
        # h5py makes it; an empty file made here would not open as HDF5
        folder = os.path.dirname(os.path.abspath(str(path)))
        if not os.access(folder, os.W_OK | os.X_OK):
            raise PermissionError("cannot create a file in %s" % folder)
        return
    fd = os.open(str(path), os.O_WRONLY | os.O_CREAT if mode == "w" else os.O_RDWR, 0o666)
    os.close(fd)


def _copy_attrs(a, b):
    """Every attribute of a onto b, its type and shape as they are."""
    for k in a.attrs:
        aid = a.attrs.get_id(k)
        b.attrs.create(k, a.attrs[k], shape=aid.shape, dtype=aid.dtype)


class BlockWriter:
    """The output files from every rank's row blocks, by parallel HDF5 or
    gathered to the first rank; see the module's docstring."""

    def __init__(self, comm, parallel=None):
        self.comm = comm
        self.parallel = parallel_hdf5(comm) if parallel is None else bool(parallel)
        if self.parallel and not h5py.get_config().mpi:
            raise ValueError("this h5py has no MPI driver, so it cannot write in parallel")
        self.route = "parallel HDF5" if self.parallel else "gathered to rank 0"

    # ---- failure
    @contextlib.contextmanager
    def guard(self, what):
        """A failure here takes the job down, since the other ranks would wait
        for this one inside a collective forever; on one rank it is raised."""
        try:
            yield
        except BaseException as exc:
            if self.comm.size == 1:
                raise
            sys.stderr.write("rank %d could not %s: %s: %s\n"
                             % (self.comm.rank, what, type(exc).__name__, exc))
            sys.stderr.flush()
            self.comm.Abort(1)

    @contextlib.contextmanager
    def _file(self, path, mode):
        """The file, opened by every rank on the parallel route and by the first
        alone otherwise (None on the others). It is closed only on success: a
        rank that fails aborts, and a close is collective."""
        if not (self.parallel or self.comm.rank == 0):
            yield None
            return
        _check_path(path, mode)
        if self.parallel:
            f = h5py.File(str(path), mode, driver="mpio", comm=self.comm)
        else:
            f = h5py.File(str(path), mode)
        try:
            yield f
        except BaseException:
            if self.comm.size == 1:
                f.close()
            raise
        f.close()

    # ---- blocks
    def _bounds(self, n):
        """(this rank's first row, the total) from every rank's row count."""
        counts = np.array(self.comm.allgather(int(n)), dtype=np.int64)
        return int(counts[:self.comm.rank].sum()), int(counts.sum())

    def _gather(self, rows):
        """Every rank's rows, in rank order, on the first rank (None on the
        others), sent a piece at a time."""
        comm = self.comm
        if comm.size == 1:
            return rows
        counts = comm.allgather(len(rows))
        rowbytes = rows.itemsize * int(np.prod(rows.shape[1:], dtype=np.int64))
        if comm.rank != 0:
            for a, b in _pieces(len(rows), rowbytes):
                comm.Send(np.ascontiguousarray(rows[a:b]), dest=0, tag=_TAG)
            return None
        out = np.empty((sum(counts),) + rows.shape[1:], rows.dtype)
        out[:counts[0]] = rows
        off = counts[0]
        for r in range(1, comm.size):
            for a, b in _pieces(counts[r], rowbytes):
                comm.Recv(out[off + a:off + b], source=r, tag=_TAG)
            off += counts[r]
        return out

    def _rows(self, rows):
        """(the rows this rank writes, their first row in the whole, the total):
        its own on the parallel route, all of them on the first rank otherwise."""
        start, total = self._bounds(len(rows))
        if self.parallel:
            return rows, start, total
        return self._gather(rows), 0, total

    # ---- the writers
    def write_checkpoint(self, path, field, cd_rows, gid_rows, values_rows, degree, ncomp,
                         dim, mode="w"):
        """A first file's field: the dof table, the cell labels and the values
        of every rank's blocks, as `divfree_clean.write_checkpoint` writes
        them."""
        with self.guard("write %s" % path):
            per = _nloc(degree, dim) * ncomp
            cd = np.ascontiguousarray(cd_rows, dtype=np.int32)
            if cd.size == 0:
                cd = cd.reshape(0, per)
            if cd.ndim != 2 or cd.shape[1] != per:
                raise ValueError("%s: '%s' dof rows of shape %s, not the (rows, %d) of a "
                                 "degree %d element of %d components in %dD"
                                 % (path, field, cd.shape, per, degree, ncomp, dim))
            gid = np.ascontiguousarray(gid_rows, dtype=np.uint64).ravel()
            if len(gid) != len(cd):
                raise ValueError("%s: '%s' has %d cell labels for %d dof rows"
                                 % (path, field, len(gid), len(cd)))
            vals = self._values(path, field, values_rows, ncomp)
            cd, c0, ncells = self._rows(cd.ravel().reshape(-1, per))
            gid = self._rows(gid)[0]
            vals, n0, nnodes = self._rows(vals)
            with self._file(path, mode) as f:
                if f is not None:
                    g = self._group(f, field, degree, ncomp, dim)
                    dcd = _create(g, "cell_dofs", ncells * per, np.int32)
                    dxc = _create(g, "x_cell_dofs", ncells + 1, np.uint64)
                    dgid = _create(g, "cells", ncells, np.uint64)
                    dv = _create(g, "vector_0", nnodes * ncomp, np.float64)
                    dv.attrs["partition"] = np.array([0], dtype=np.uint64)
                    _put(dcd, c0 * per, cd.ravel())
                    _put(dgid, c0, gid)
                    _put(dv, n0 * ncomp, vals.ravel())
                    m = len(cd)
                    last = self.comm.rank == self.comm.size - 1 or not self.parallel
                    for a, b in _pieces(m + (1 if last else 0), 8):
                        dxc[c0 + a:c0 + b] = np.arange(c0 + a, c0 + b, dtype=np.uint64) * per
            self._done()

    def write_vector(self, path, field, values_rows, degree, ncomp, dim, mode="w"):
        """A later stamp or component: `<field>/vector_0` and the signature, as
        `divfree_clean.write_vector` writes them."""
        with self.guard("write %s" % path):
            vals = self._values(path, field, values_rows, ncomp)
            vals, n0, nnodes = self._rows(vals)
            with self._file(path, mode) as f:
                if f is not None:
                    g = self._group(f, field, degree, ncomp, dim)
                    _put(_create(g, "vector_0", nnodes * ncomp, np.float64), n0 * ncomp,
                         vals.ravel())
            self._done()

    def _values(self, path, field, values_rows, ncomp):
        """A rank's values as float64 rows of ncomp."""
        v = np.ascontiguousarray(values_rows, dtype=np.float64)
        if v.size == 0:
            v = v.reshape(0, ncomp)
        if v.ndim == 1 and ncomp == 1:
            v = v[:, None]
        if v.ndim != 2 or v.shape[1] != ncomp:
            raise ValueError("%s: '%s' value rows of shape %s for an element of %d components"
                             % (path, field, v.shape, ncomp))
        return v

    def _group(self, f, field, degree, ncomp, dim):
        g = f.create_group(field)
        g.attrs["signature"] = np.bytes_(_clean()._signature(degree, ncomp, dim))
        return g

    def _done(self):
        """The gathered route returns on every rank once the first has written."""
        if not self.parallel and self.comm.size > 1:
            self.comm.Barrier()

    # ---- a field carried through
    def copy_group(self, src, dst, name):
        """`src[name]`, a group, copied into dst, datasets and attributes as
        they are. False, on every rank, where src holds no such group. On the
        parallel route every rank reads and writes its share of each dataset's
        first axis; datasets come out contiguous whatever their chunking and
        filters were. Soft and external links, and a second hard link to an
        object, are not reproduced; a dolfin field group holds none."""
        what = "copy '%s' of %s into %s" % (name, src, dst)
        if not self.parallel:
            ok = None
            if self.comm.rank == 0:
                with self.guard(what):
                    ok = _clean().copy_group(src, dst, name)
            return bool(self.comm.bcast(ok, root=0))
        with self.guard(what):
            with h5py.File(str(src), "r") as a:
                has = self.comm.allgather(name in a)
                if not any(has):
                    return False
                if name not in a:
                    raise ValueError("%s holds no '%s' where rank %d finds one"
                                     % (src, name, has.index(True)))
                items = [(name, a[name])]
                a[name].visititems(lambda n, o: items.append((name + "/" + n, o)))
                with self._file(dst, "a") as b:
                    for n, o in items:
                        if isinstance(o, h5py.Group):
                            _copy_attrs(o, b.create_group(n))
                            continue
                        if h5py.check_vlen_dtype(o.dtype) is not None or o.dtype.kind == "O":
                            raise ValueError("%s: '%s' holds variable-length data, which "
                                             "parallel HDF5 does not write; %s=1 copies it"
                                             % (src, n, GATHER_ENV))
                        dcpl = h5py.h5p.create(h5py.h5p.DATASET_CREATE)
                        dcpl.set_fill_time(h5py.h5d.FILL_TIME_NEVER)
                        d = b.create_dataset(n, o.shape, dtype=o.dtype, dcpl=dcpl)
                        _copy_attrs(o, d)
                        if o.shape == ():
                            if self.comm.rank == 0:
                                d[()] = o[()]
                            continue
                        r0, r1 = _clean()._split(o.shape[0], self.comm.size, self.comm.rank)
                        rowbytes = o.dtype.itemsize * int(np.prod(o.shape[1:], dtype=np.int64))
                        for p, q in _pieces(r1 - r0, rowbytes):
                            d[r0 + p:r0 + q] = o[r0 + p:r0 + q]
        return True

    # ---- text
    def write_text(self, path, text):
        """A text file -- the stamp list, the parameter file -- written by the
        first rank; every rank returns once it is there."""
        if self.comm.rank == 0:
            with self.guard("write %s" % path):
                with open(str(path), "w") as f:
                    f.write(text)
        if self.comm.size > 1:
            self.comm.Barrier()
