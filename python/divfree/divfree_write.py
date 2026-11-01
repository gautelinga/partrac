"""The cleaner's output files written by every rank from its own rows.

dolfin's checkpoint layout. Each rank passes its contiguous block of the rows in
serial order (cells `[c0, c1)`, nodes `[n0, n1)`), possibly none; every call of
`BlockWriter` is collective and returns once the file is complete. Written
through h5py's MPI driver where it has one and there is more than one rank,
else gathered through the first rank; PARTRAC_HDF5_GATHER=1 or true forces the
gather. A failure on one rank aborts the job (`guard`).
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


def _piece_rows(rowbytes):
    """The rows a piece of at most PIECE_BYTES holds."""
    return max(1, PIECE_BYTES // max(int(rowbytes), 1))


def _ceil_div(n, step):
    """n / step rounded up: the pieces of step rows that n rows take."""
    return -(-int(n) // step)


def _pieces(n, rowbytes):
    """[a, b) ranges of n rows, each at most PIECE_BYTES."""
    step = _piece_rows(rowbytes)
    return [(a, min(a + step, n)) for a in range(0, n, step)]


def _nloc(degree, dim):
    """The element's nodes in a cell."""
    return dim + 1 if degree == 1 else (dim + 1) * (dim + 2) // 2


def _create(g, name, n, dtype):
    """A contiguous one-dimensional dataset, never filled: every entry is
    written. None where g is: a rank that holds no file."""
    if g is None:
        return None
    dcpl = h5py.h5p.create(h5py.h5p.DATASET_CREATE)
    dcpl.set_fill_time(h5py.h5d.FILL_TIME_NEVER)
    return g.create_dataset(name, (n,), dtype=dtype, dcpl=dcpl)


def _write(dset, a, b, data):
    """data into dset's rows [a, b); with none, a write of nothing, which is a
    rank's share of a collective write it has no rows for."""
    if b > a:
        dset[a:b] = data
        return
    fs = dset.id.get_space()
    fs.select_none()
    ms = h5py.h5s.create_simple((1,))
    ms.select_none()
    dxpl = h5py.h5p.create(h5py.h5p.DATASET_XFER)
    dxpl.set_dxpl_mpio(h5py.h5fd.MPIO_COLLECTIVE)
    dset.id.write(ms, fs, np.empty(1, dset.dtype), dxpl=dxpl)


def _collective(dset, lo, hi, step, most, rows):
    """dset's rows [lo, hi) from rows(a, b), in `most` collective writes of at
    most step rows, as many as every other rank makes."""
    with dset.collective:
        for i in range(most):
            a, b = min(lo + i * step, hi), min(lo + (i + 1) * step, hi)
            _write(dset, a, b, rows(a, b))


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

    # ---- rows
    def _counts(self, *n):
        """Every rank's counts, a row a rank, from one exchange."""
        return np.array(self.comm.allgather([int(k) for k in n]), np.int64)

    def _put(self, dset, rows, counts, width=1):
        """Every rank's rows into dset, in rank order: rows this rank's, counts
        every rank's, width the dataset's entries a row. On the parallel route
        each rank writes its own; otherwise the first rank writes them all,
        receiving the others' a piece at a time."""
        comm, rank = self.comm, self.comm.rank
        rowbytes = rows.itemsize * int(np.prod(rows.shape[1:], dtype=np.int64))
        step = _piece_rows(rowbytes)
        start = np.concatenate([[0], np.cumsum(counts)])
        flat = rows.reshape(-1)
        if self.parallel:
            most = max(_ceil_div(c, step) for c in counts)
            s0 = int(start[rank])
            _collective(dset, s0 * width, (s0 + len(rows)) * width, step * width, most,
                        lambda a, b: flat[a - s0 * width:b - s0 * width])
            return
        if rank:
            for a, b in _pieces(len(rows), rowbytes):
                comm.Send(rows[a:b], dest=0, tag=_TAG)
            return
        buf = np.empty((min(step, int(counts[1:].max(initial=0))),) + rows.shape[1:],
                       rows.dtype)
        for r in range(comm.size):
            for a, b in _pieces(int(counts[r]), rowbytes):
                piece = rows[a:b] if r == 0 else buf[:b - a]
                if r:
                    comm.Recv(piece, source=r, tag=_TAG)
                dset[(start[r] + a) * width:(start[r] + b) * width] = piece.reshape(-1)

    def _offsets(self, dset, counts, per):
        """`x_cell_dofs`, each writing rank its own slab of it and the last the
        final row: row k is k*per."""
        size, total = self.comm.size, int(counts.sum())
        if not self.parallel:
            if dset is not None:
                for a, b in _pieces(total + 1, 8):
                    dset[a:b] = np.arange(a, b, dtype=np.uint64) * per
            return
        m = counts + (np.arange(size) == size - 1)
        lo = int(counts[:self.comm.rank].sum())
        step = _piece_rows(8)
        _collective(dset, lo, lo + int(m[self.comm.rank]), step,
                    max(_ceil_div(k, step) for k in m),
                    lambda a, b: np.arange(a, b, dtype=np.uint64) * per)

    # ---- the writers
    def write_checkpoint(self, path, field, cd_rows, gid_rows, values_rows, degree, ncomp,
                         dim, mode="w"):
        """A first file's field: the dof table, the cell labels and the values
        of every rank's blocks, in dolfin's checkpoint layout. The dof table is
        int64 where any rank's rows are, as `_cell_dofs` makes them past
        int32, and int32 otherwise."""
        with self.guard("write %s" % path):
            per = _nloc(degree, dim) * ncomp
            wide = np.asarray(cd_rows).dtype.itemsize > 4
            cd = np.ascontiguousarray(cd_rows, dtype=np.int64 if wide else np.int32)
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
            counts = self._counts(len(cd), len(vals))
            ncells, nnodes = (int(k) for k in counts.sum(axis=0))
            if self.comm.allreduce(wide, op=MPI.LOR):
                cd = cd.astype(np.int64, copy=False)
            with self._file(path, mode) as f:
                g = None if f is None else self._group(f, field, degree, ncomp, dim)
                dcd = _create(g, "cell_dofs", ncells * per, cd.dtype)
                dxc = _create(g, "x_cell_dofs", ncells + 1, np.uint64)
                dgid = _create(g, "cells", ncells, np.uint64)
                dv = _create(g, "vector_0", nnodes * ncomp, np.float64)
                if dv is not None:
                    dv.attrs["partition"] = np.array([0], dtype=np.uint64)
                self._put(dcd, cd, counts[:, 0], per)
                self._put(dgid, gid, counts[:, 0])
                self._put(dv, vals, counts[:, 1], ncomp)
                self._offsets(dxc, counts[:, 0], per)
            self._done()

    def write_vector(self, path, field, values_rows, degree, ncomp, dim, mode="w"):
        """A later stamp or component: `<field>/vector_0` and the signature, all
        a loader opens of it."""
        with self.guard("write %s" % path):
            vals = self._values(path, field, values_rows, ncomp)
            counts = self._counts(len(vals))
            with self._file(path, mode) as f:
                g = None if f is None else self._group(f, field, degree, ncomp, dim)
                self._put(_create(g, "vector_0", int(counts.sum()) * ncomp, np.float64),
                          vals, counts[:, 0], ncomp)
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
                        block = [_clean()._split(o.shape[0], self.comm.size, r)
                                 for r in range(self.comm.size)]
                        step = _piece_rows(o.dtype.itemsize
                                           * int(np.prod(o.shape[1:], dtype=np.int64)))
                        r0, r1 = block[self.comm.rank]
                        most = max(_ceil_div(j1 - j0, step) for j0, j1 in block)
                        _collective(d, r0, r1, step, most, lambda j0, j1: o[j0:j1])
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
