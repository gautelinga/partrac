import argparse
import numpy as np
import os
import h5py
from scipy import signal
from tqdm import tqdm
from helpers import mpi_print, MPI, rank, size, comm, select_timestamps, make_folder_safe, get_fluid_domain
import itertools


def parse_args():
    parser = argparse.ArgumentParser(description="Plot pos")
    parser.add_argument("folder", type=str, help="Folder")
    parser.add_argument("-t0", type=float, default=None, help="")
    parser.add_argument("-t1", type=float, default=None, help="")
    parser.add_argument("-axis", type=int, default=0, help="Axis")
    parser.add_argument("-fmax", type=float, default=0, help="Cut-off frequency of FFT")
    parser.add_argument("--show", action="store_true", help="Show plot")
    parser.add_argument("--phase_by_density", action="store_true", help="Show plot")
    parser.add_argument("--blocks", type=int, default=1, help="Number of blocks to process")
    parser.add_argument("-o", "--output_dir", type=str, default="", help="Output directory")
    args = parser.parse_args()
    return args


def load_mesh(nx, ny, is_fluid):
    mesh = None
    if rank == 0:
        nodes = np.array([(i, j) for j,  i in itertools.product(range(ny+1), range(nx+1))], dtype=float)
        elems = np.array([(i + j*(nx+1), i+1 + j*(nx+1), i+1 + (j+1)*(nx+1), i + (j+1)*(nx+1)) 
                          for j, i in itertools.product(range(ny), range(nx))], dtype=int)
        
        elems = elems[is_fluid, :]
        used_nodes = np.unique(elems)
        map_ids = np.zeros(used_nodes.max()+1, dtype=int)
        for i, j in zip(used_nodes, range(len(used_nodes))):
            map_ids[i] = j
        nodes = nodes[used_nodes, :]
        elems = map_ids[elems]
        mesh = (elems, nodes)
    return mesh


class XDMFTimeSeries:
    def __init__(self, xdmffilename, points, cells):
        assert(xdmffilename[-5:] == ".xdmf")
        assert(points.shape[1] == 2)
        assert(cells.shape[1] == 4)
        self.xdmffilename = xdmffilename
        self.h5filename = xdmffilename[:-5] + ".h5"
        self.h5filename_base = os.path.basename(self.h5filename)
        self.points = np.array(points)
        self.cells = np.array(cells)

        self.xdmffile = open(self.xdmffilename, "w")
        self.h5file = h5py.File(self.h5filename, "w")

        self.xdmffile.write(self._head())
        self.h5file.create_dataset("mesh/cells", data=self.cells)
        self.h5file.create_dataset("mesh/points", data=self.points)

        self._COUNT = 0

    def close(self):
        self.xdmffile.write(self._tail())
        self.xdmffile.close()
        self.h5file.close()

    def write(self, cell_datasets, t):
        self.xdmffile.write(self._body_top(t))
        for name, cell_data in cell_datasets:
            assert(cell_data.shape[0] == self.cells.shape[0])
            if len(cell_data.shape) == 1 or cell_data.shape[1] == 1:
                self.xdmffile.write(self._body_scalar(self._COUNT, name))
                self.h5file.create_dataset("{i}/{name}".format(i=self._COUNT, name=name), data=cell_data)
            elif cell_data.shape[1] == 2:
                self.xdmffile.write(self._body_vector(self._COUNT, name))
                self.h5file.create_dataset("{i}/{name}".format(i=self._COUNT, name=name), data=cell_data)
            else:
                exit("UNSUPPORTED")
        self.xdmffile.write(self._body_btm())
        self._COUNT += 1

    def _head(self):
        return """<Xdmf xmlns:ns0="http://www.w3.org/2003/XInclude" Version="3.0">
    <Domain>
        <Grid Name="TimeSeries_partrac" GridType="Collection" CollectionType="Temporal">"""

    def _body_top(self, t):
        return """
            <Grid>
                <ns0:include xpointer="xpointer(//Grid[@Name=&quot;mesh&quot;]/*[self::Topology or self::Geometry])" />
                <Time Value="{time}" />""".format(time=t)

    def _body_vector(self, i, name):
        return """
                <Attribute Name="{name}" AttributeType="Vector" Center="Cell">
                    <DataItem DataType="Float" Dimensions="{num_cells} 2" Format="HDF" Precision="8">{h5filename}:/{i}/{name}</DataItem>
                </Attribute>""".format(name=name, i=i, num_points=len(self.points), num_cells=len(self.cells), h5filename=self.h5filename_base)

    def _body_scalar(self, i, name):
        return """
                <Attribute Name="{name}" AttributeType="Scalar" Center="Cell">
                    <DataItem DataType="Float" Dimensions="{num_cells}" Format="HDF" Precision="8">{h5filename}:/{i}/{name}</DataItem>
                </Attribute>""".format(name=name, i=i, num_points=len(self.points), num_cells=len(self.cells), h5filename=self.h5filename_base)

    def _body_btm(self):
        return """
            </Grid>"""

    def _tail(self):
        return """
        </Grid>
        <Grid Name="mesh" GridType="Uniform">
            <Geometry GeometryType="XY">
                <DataItem DataType="Float" Dimensions="{num_points} 2" Format="HDF" Precision="8">{h5filename}:/mesh/points</DataItem>
            </Geometry>
            <Topology TopologyType="Quadrilateral" NumberOfElements="843710">
                <DataItem DataType="Int" Dimensions="{num_cells} 4" Format="HDF" Precision="8">{h5filename}:/mesh/cells</DataItem>
            </Topology>
        </Grid>
    </Domain>
</Xdmf>""".format(num_points=len(self.points), num_cells=len(self.cells), h5filename=self.h5filename_base)

if __name__ == "__main__":
    args = parse_args()

    felbm_folder = args.folder
    if args.output_dir != "":
        output_folder = args.output_dir
    else:
        output_folder = felbm_folder
    analysisfolder = os.path.join(output_folder, "Analysis")
    make_folder_safe(analysisfolder)

    timestamps = select_timestamps(felbm_folder, args.t0, args.t1)

    is_fluid_xy, (nx, ny, nz) = get_fluid_domain(felbm_folder)
    is_fluid = is_fluid_xy.flatten()
    mesh = load_mesh(nx, ny, is_fluid)
    
    is_domain_xyz = np.zeros((nz, ny, nx), dtype=bool)
    is_domain_xyz[nz // 2, :, :] = is_fluid_xy
    is_domain = is_domain_xyz.flatten()
    fluid2xy = np.argwhere(is_fluid).flatten()
    mpi_print(len(fluid2xy), fluid2xy[-1])

    Ndofs = int(np.sum(is_domain))

    rho = np.zeros_like(is_fluid, dtype=float)
    rho_loc = np.zeros_like(rho)

    ux = np.zeros_like(is_fluid, dtype=float)
    ux_loc = np.zeros_like(ux)

    uy = np.zeros_like(is_fluid, dtype=float)
    uy_loc = np.zeros_like(uy)

    du = np.zeros_like(rho)
    du_loc = np.zeros_like(du)

    freq = np.zeros_like(ux)
    freq_loc = np.zeros_like(freq)

    freqy = np.zeros_like(ux)
    freqy_loc = np.zeros_like(freqy)

    Nt = len(timestamps)

    mpi_print("porosity={}".format(Ndofs/(nx*ny)))
    mpi_print("Is fluid: {}".format(Ndofs))
    block_size = int(np.ceil(Ndofs / args.blocks))

    blocks = list(range(args.blocks))
    bstart = args.blocks * rank // size
    bstop = args.blocks * (rank + 1) // size
    blocks_proc = blocks[bstart:bstop]

    t_ = np.array([t for t, _ in timestamps])

    dt = t_[1]-t_[0]
    
    f = np.fft.rfftfreq(Nt, dt)
    Nf = len(f)
    assert(Nf == Nt // 2 + 1)
    df = f[1]-f[0]

    if args.fmax > 0:
        ff = f[f <= args.fmax]
        Nff = len(ff)
    else:
        ff = f
        Nff = Nf
    mpi_print("Number of frequencies to keep: {}".format(Nff))

    Fux = np.zeros((len(is_fluid), Nff), dtype=complex)
    Fuy = np.zeros((len(is_fluid), Nff), dtype=complex)
    Fux_loc = np.zeros_like(Fux)
    Fuy_loc = np.zeros_like(Fuy)
    
    P2ux = np.zeros(Nf)
    P2ux_loc = np.zeros(Nf)

    P2uy = np.zeros(Nf)
    P2uy_loc = np.zeros(Nf)

    comm.Barrier()
    if rank == 0:
        pbar = tqdm(total=len(blocks_proc)*len(timestamps))

    for iblock in blocks_proc:
        istart = iblock * block_size
        istop = min((iblock + 1) * block_size, Ndofs)
        uxt_block = np.zeros((istop-istart, Nt))
        uyt_block = np.zeros_like(uxt_block)
        rho_block = np.zeros_like(uxt_block)
        #Ux = np.zeros((Nt, ny, nx))
        #Uy = np.zeros_like(Ux)
        for it, timestamp in enumerate(timestamps):
            t = timestamp[0]
            h5filename = timestamp[1]
            with h5py.File(os.path.join(felbm_folder, h5filename), "r") as h5f:
                #p[:, :] = np.array(h5f["pressure"]).reshape((nz, ny, nx))[nz // 2, :, :]
                #rho[:, :] = np.array(h5f["density"]).reshape((nz, ny, nx))[nz // 2, :, :]
                #ux[:, :] = np.array(h5f["u_x"]).reshape((nz, ny, nx))[nz // 2, :, :]
                #uy[:, :] = np.array(h5f["u_y"]).reshape((nz, ny, nx))[nz // 2, :, :]
                rho_block[:, it] = np.array(h5f["density"]).flatten()[is_domain][istart:istop]
                uxt_block[:, it] = np.array(h5f["u_x"]).flatten()[is_domain][istart:istop]
                uyt_block[:, it] = np.array(h5f["u_y"]).flatten()[is_domain][istart:istop]
                #print(w)
                if rank == 0:
                    pbar.update(1)
        freq_block = np.zeros(uxt_block.shape[0])
        freqy_block = np.zeros(uxt_block.shape[0])
        rhoavg_block = rho_block.mean(axis=1)

        Fux_block = np.zeros((uxt_block.shape[0], Nff), dtype=complex)
        Fuy_block = np.zeros((uxt_block.shape[0], Nff), dtype=complex)

        for i in range(uxt_block.shape[0]):
            _Fux = np.fft.rfft(uxt_block[i, :])
            _Fuy = np.fft.rfft(uyt_block[i, :])
            P2ux_i = (np.conjugate(_Fux) * _Fux).real
            P2uy_i = (np.conjugate(_Fuy) * _Fuy).real

            P2ux_loc[:] += P2ux_i
            P2uy_loc[:] += P2uy_i

            Pux_i = np.sqrt(P2ux_i)
            Puy_i = np.sqrt(P2uy_i)

            sumPux = np.sum(Pux_i)
            sumPuy = np.sum(Puy_i)

            freq_block[i] = np.sum(f*Pux_i)/sumPux if sumPux > 0. else 0.
            freqy_block[i] = np.sum(f*Puy_i)/sumPuy if sumPuy > 0. else 0.

            Fux_block[i, :] = _Fux[:Nff]
            Fuy_block[i, :] = _Fuy[:Nff]

        ids = fluid2xy[istart:istop]
        ux_loc[ids] = np.sqrt(np.mean(uxt_block**2, axis=1))
        uy_loc[ids] = np.sqrt(np.mean(uyt_block**2, axis=1))
        freq_loc[ids] = freq_block[:]
        freqy_loc[ids] = freqy_block[:]

        rho_loc[ids] = rhoavg_block[:]
        du_loc[ids] = np.sqrt(np.mean((uxt_block - np.outer(uxt_block.mean(axis=1), np.ones(uxt_block.shape[1])))**2 
                                      + (uyt_block - np.outer(uyt_block.mean(axis=1), np.ones(uyt_block.shape[1])))**2, axis=1))
        
        Fux_loc[ids, :] = Fux_block[:, :]
        Fuy_loc[ids, :] = Fuy_block[:, :]

    if rank == 0:
        pbar.close()


    comm.Reduce(ux_loc, ux, op=MPI.SUM, root=0)
    comm.Reduce(uy_loc, uy, op=MPI.SUM, root=0)
    comm.Reduce(freq_loc, freq, op=MPI.SUM, root=0)
    comm.Reduce(freqy_loc, freqy, op=MPI.SUM, root=0)
    
    comm.Reduce(rho_loc, rho, op=MPI.SUM, root=0)
    comm.Reduce(du_loc, du, op=MPI.SUM, root=0)

    comm.Reduce(Fux_loc, Fux, op=MPI.SUM, root=0)
    comm.Reduce(Fuy_loc, Fuy, op=MPI.SUM, root=0)
    
    comm.Reduce(P2ux_loc, P2ux, op=MPI.SUM, root=0)
    comm.Reduce(P2uy_loc, P2uy, op=MPI.SUM, root=0)

    comm.Barrier()
    if rank == 0:
        cell_data = {
            "f": [freq[is_fluid]],
            "fy": [freqy[is_fluid]],
            "ux": [ux[is_fluid]],
            "uy": [uy[is_fluid]],
            "rho": [rho[is_fluid]],
            "du": [du[is_fluid]]
        }

        import meshio
        elems, nodes = mesh
        m = meshio.Mesh(
            nodes, [("quad", elems)], cell_data=cell_data
        )
        m.write(
            os.path.join(analysisfolder, "data_avg.xdmf"))

        pbar2 = tqdm(total=Nff)
        
        ts = XDMFTimeSeries(os.path.join(analysisfolder, "data_f.xdmf"), nodes, elems)
        for i, fi in enumerate(ff):
            pbar2.update(1)

            _Fuxi = Fux[is_fluid, i]
            _Fuyi = Fuy[is_fluid, i]

            cell_datasets = [
                ("Fu_real", np.vstack((_Fuxi.real, _Fuyi.real)).T),
                ("Fu_imag", np.vstack((_Fuxi.imag, _Fuyi.imag)).T) #,
                # ("P2", P2[is_fluid, i])
            ]
            ts.write(cell_datasets, fi)

        ts.close()
        pbar2.close()

    if rank == 0:
        freq = freq.reshape((ny, nx))
        freqy = freqy.reshape((ny, nx))
        ux = ux.reshape((ny, nx))
        uy = uy.reshape((ny, nx))
        rho = rho.reshape((ny, nx))
        du = du.reshape((ny, nx))

        np.savetxt(os.path.join(analysisfolder, "uxnorm_avg.dat"), ux)
        np.savetxt(os.path.join(analysisfolder, "uynorm_avg.dat"), uy)
        np.savetxt(os.path.join(analysisfolder, "freq_avg.dat"), freq)
        np.savetxt(os.path.join(analysisfolder, "freqy_avg.dat"), freqy)

        np.savetxt(os.path.join(analysisfolder, "rho_avg.dat"), rho)
        np.savetxt(os.path.join(analysisfolder, "is_fluid_xy.dat"), is_fluid_xy)
        np.savetxt(os.path.join(analysisfolder, "du.dat"), du)

        #P2ux /= P2ux.sum() * df
        #P2uy /= P2uy.sum() * df
        #assert(abs((P2ux*df).sum()-1) < 1e-5)
        P2f_data = np.vstack((f, P2ux, P2uy)).T
        np.savetxt(os.path.join(analysisfolder, "P2f.dat"), P2f_data)
