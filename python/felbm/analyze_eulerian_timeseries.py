import argparse
from meshtools import numpy_to_dolfin
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
    parser.add_argument("-fmax", type=float, default=0, help="Cut-off frequency")
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
    #params = Params(args.folder)
    #t0 = params.get_tmin()
    #Lx = params.get("Lx", t0)

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


    #is_solid = is_solid[nz // 2, :, :].flatten()
    #is_fluid = np.logical_not(is_solid)
    #is_fluid = np.logical_not(is_solid)
    
    is_domain_xyz = np.zeros((nz, ny, nx), dtype=bool)
    is_domain_xyz[nz // 2, :, :] = is_fluid_xy
    is_domain = is_domain_xyz.flatten()
    fluid2xy = np.argwhere(is_fluid).flatten()
    mpi_print(len(fluid2xy), fluid2xy[-1])

    Ndofs = int(np.sum(is_domain))

    #print(np.sum(is_domain_flat), np.sum(is_fluid_xy_flat))
    #p = np.zeros_like(is_solid, dtype=float)
    rho = np.zeros_like(is_fluid, dtype=float)
    rho_loc = np.zeros_like(rho)
    #c = np.zeros_like(is_solid, dtype=float)
    #ux = np.zeros_like(is_solid, dtype=float)
    #uy = np.zeros_like(is_solid, dtype=float)
    ux = np.zeros_like(is_fluid, dtype=float)
    ux_loc = np.zeros_like(ux)
    du = np.zeros_like(rho)
    du_loc = np.zeros_like(du)
    # x, y = np.meshgrid(np.arange(nx) + .5, np.arange(ny) + 0.5)
    freq = np.zeros_like(ux)
    freq_loc = np.zeros_like(freq)

    Nt = len(timestamps)

    #print(is_domain.shape)

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
    P2 = np.zeros((len(is_fluid), Nff))

    Fux_loc = np.zeros_like(Fux)
    Fuy_loc = np.zeros_like(Fuy)
    P2_loc = np.zeros_like(P2)

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
        rhoavg_block = rho_block.mean(axis=1)

        Fux_block = np.zeros((uxt_block.shape[0], Nff), dtype=complex)
        Fuy_block = np.zeros((uxt_block.shape[0], Nff), dtype=complex)
        P2_block = np.zeros((uxt_block.shape[0], Nff))

        for i in range(uxt_block.shape[0]):
            #f_i, P2_i = signal.welch(uxt_block[i, :], df, 'flattop', Nt, scaling='spectrum')
            #f_i, P2_i = signal.periodogram(uxt_block[i, :], df, nfft=Nt)
            _Fux = np.fft.rfft(uxt_block[i, :])
            _Fuy = np.fft.rfft(uyt_block[i, :])
            P2_i = (np.conjugate(_Fux) * _Fux).real
            #plt.plot(f_i, P2_i)
            #plt.plot(f_i, P2_a)
            #plt.show()
            P_i = np.sqrt(P2_i)
            sumP = np.sum(P_i)
            freq_block[i] = np.sum(f*P_i)/sumP if sumP > 0. else 0.
            #plt.plot(t_, uxt_block[0, :])
            #if (all(f_i == f)) and False:
            #    print("f_i:", f_i)
            #    print("f:  ", f)
            #    print("uxt:", uxt_block[i, :])
            #    exit()
            Fux_block[i, :] = _Fux[:Nff]
            Fuy_block[i, :] = _Fuy[:Nff]
            
            #P2_block[:] += P2
            P2_block[i, :] += P2_i[:Nff]


            #ax.plot(f_i, P2_i)
        #plt.show()

        ids = fluid2xy[istart:istop]
        ux_loc[ids] = np.sqrt(np.mean(uxt_block**2, axis=1))
        freq_loc[ids] = freq_block[:]
        rho_loc[ids] = rhoavg_block[:]
        du_loc[ids] = np.sqrt(np.mean((uxt_block - np.outer(uxt_block.mean(axis=1), np.ones(uxt_block.shape[1])))**2 
                                      + (uyt_block - np.outer(uyt_block.mean(axis=1), np.ones(uyt_block.shape[1])))**2, axis=1))
        
        Fux_loc[ids, :] = Fux_block[:, :]
        Fuy_loc[ids, :] = Fuy_block[:, :]
        P2_loc[ids, :] = P2_block[:, :]

    if rank == 0:
        pbar.close()


    comm.Reduce(ux_loc, ux, op=MPI.SUM, root=0)
    comm.Reduce(freq_loc, freq, op=MPI.SUM, root=0)
    comm.Reduce(rho_loc, rho, op=MPI.SUM, root=0)
    comm.Reduce(du_loc, du, op=MPI.SUM, root=0)

    comm.Reduce(Fux_loc, Fux, op=MPI.SUM, root=0)
    comm.Reduce(Fuy_loc, Fuy, op=MPI.SUM, root=0)
    comm.Reduce(P2_loc, P2, op=MPI.SUM, root=0)

    comm.Barrier()
    if rank == 0:
        cell_data = {
            "f": [freq[is_fluid]], 
            "ux": [ux[is_fluid]],
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
                ("Fu_imag", np.vstack((_Fuyi.imag, _Fuyi.imag)).T),
                ("P2", P2[is_fluid, i])
            ]
            ts.write(cell_datasets, fi)

        ts.close()
        pbar2.close()

    if rank == 0:
        freq = freq.reshape((ny, nx))
        ux = ux.reshape((ny, nx))
        rho = rho.reshape((ny, nx))
        du = du.reshape((ny, nx))

        np.savetxt(os.path.join(analysisfolder, "uxnorm_avg.dat"), ux)
        np.savetxt(os.path.join(analysisfolder, "freq_avg.dat"), freq)
        np.savetxt(os.path.join(analysisfolder, "rho_avg.dat"), rho)
        np.savetxt(os.path.join(analysisfolder, "is_fluid_xy.dat"), is_fluid_xy)
        np.savetxt(os.path.join(analysisfolder, "du.dat"), du)

        P2f = P2[is_fluid, :].sum(axis=0)
        P2f /= P2f.sum() * df
        P2fdata = np.vstack((ff, P2f)).T
        np.savetxt(os.path.join(analysisfolder, "P2f.dat"), P2fdata)
