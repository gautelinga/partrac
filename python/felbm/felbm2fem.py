from numpy import unique
from analyze_eulerian_timeseries import *

import matplotlib.pyplot as plt
import dolfin as df
import itertools


def parse_args():
    parser = argparse.ArgumentParser(description="FEM interpolated FELBM data")
    parser.add_argument("folder", type=str, help="Folder")
    parser.add_argument("-t0", type=float, default=None, help="")
    parser.add_argument("-t1", type=float, default=None, help="")
    parser.add_argument("-axis", type=int, default=0, help="Axis")
    parser.add_argument("--show", action="store_true", help="Show plot")
    parser.add_argument("--phase_by_density", action="store_true", help="Show plot")
    args = parser.parse_args()
    return args


class PBC(df.SubDomain):
    def __init__(self, Lx, Ly):
        self.Lx = Lx
        self.Ly = Ly
        df.SubDomain.__init__(self)

    def inside(self, x, on_boundary):
        # return True if on left or bottom boundary AND NOT on one of the two slave edges
        return bool(
            (df.near(x[0], 0.) or df.near(x[1], 0.))
            and (not (df.near(x[0], self.Lx) or df.near(x[1], self.Ly)))
            and on_boundary)

    def map(self, x, y):
        if df.near(x[0], self.Lx) and df.near(x[1], self.Ly):
            y[0] = x[0] - self.Lx
            y[1] = x[1] - self.Ly
        elif df.near(x[0], self.Lx):
            y[0] = x[0] - self.Lx
            y[1] = x[1]
        else:  # near(x[2], Lz/2.):
            y[0] = x[0]
            y[1] = x[1] - self.Ly

class Wall(df.SubDomain):
    def __init__(self, is_fluid_xy):
        self.is_fluid_xy = is_fluid_xy
        self.ny, self.nx = is_fluid_xy.shape
        super().__init__()

    def inside(self, x, on_boundary):
        if bool(x[0] > df.DOLFIN_EPS_LARGE and x[0] < self.nx - df.DOLFIN_EPS_LARGE and
                x[1] > df.DOLFIN_EPS_LARGE and x[1] < self.ny - df.DOLFIN_EPS_LARGE):
            return on_boundary
        elif on_boundary:
            return x[0].is_integer() and x[1].is_integer()
        
def make_xdmf_mesh(nx, ny, is_fluid, tmpfilename):
    if rank == 0:
        nodes = np.array([(i-0.5, j-0.5) for j, i in itertools.product(range(ny+1), range(nx+1))], dtype=float)
        elems = np.array([(i + j*(nx+1), i+1 + j*(nx+1), i+1 + (j+1)*(nx+1), i + (j+1)*(nx+1)) 
                          for j, i in itertools.product(range(ny), range(nx))], dtype=int)
        
        elems = elems[is_fluid, :]
        used_nodes = np.unique(elems)
        map_ids = np.zeros(used_nodes.max()+1, dtype=int)
        for i, j in zip(used_nodes, range(len(used_nodes))):
            map_ids[i] = j
        nodes = nodes[used_nodes, :]
        elems = map_ids[elems]

        import meshio
        m = meshio.Mesh(nodes, [("quad", elems)])
        m.write(tmpfilename)

    comm.Barrier()


def load_mesh(nx, ny, is_fluid, analysisfolder):
    tmpfilename = os.path.join(analysisfolder, "foo.xdmf")
    make_xdmf_mesh(nx, ny, is_fluid, tmpfilename)

    mesh = df.Mesh()
    with df.XDMFFile(tmpfilename) as xdmff:
        xdmff.read(mesh)

    return mesh



if __name__ == "__main__":
    args = parse_args()

    felbm_folder = args.folder
    timestamps = select_timestamps(felbm_folder, args.t0, args.t1)

    analysisfolder = os.path.join(felbm_folder, "Analysis")
    make_folder_safe(analysisfolder)

    is_fluid_xy, (nx, ny, nz) = get_fluid_domain(felbm_folder)
    is_fluid = is_fluid_xy.flatten()

    #print("nx = {}, ny = {}".format(nx, ny))

    rho = np.zeros_like(is_fluid_xy, dtype=float)
    ux = np.zeros_like(rho)
    uy = np.zeros_like(rho)
    #p = np.zeros_like(rho)
    ux_ = ux.flatten()[is_fluid]
    uy_ = np.zeros_like(ux_)

    mesh = load_mesh(nx, ny, is_fluid, analysisfolder)

    pbc = PBC(nx, ny)

    subd = df.MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
    wall = Wall(is_fluid_xy)
    wall.mark(subd, 1)
    pbc.mark(subd, 2)
    with df.XDMFFile(mesh.mpi_comm(), os.path.join(analysisfolder, "subd.xdmf")) as xdmff:
        xdmff.write(subd)

    S0 = df.FunctionSpace(mesh, "DG", 0, constrained_domain=pbc)

    dof_coord_loc = np.array(S0.tabulate_dof_coordinates()-0.0, dtype=int)
    dof_coord_glob = np.array([(i, j) for j, i in itertools.product(range(ny), range(nx))], dtype=int)[is_fluid, :]
    coord2glob = dict()
    for i, coord in enumerate(dof_coord_glob):
        coord2glob[tuple(coord)] = i
    loc2glob = np.array([coord2glob[tuple(coord)] for coord in dof_coord_loc], dtype=int)

    ux0 = df.Function(S0, name="ux0")
    uy0 = df.Function(S0, name="uy0")
    u0 = df.as_vector((ux0, uy0))

    #u_el = df.FiniteElement("RT", mesh.ufl_cell(), 1)
    u_el = df.VectorElement("CG", mesh.ufl_cell(), 1)
    V = df.FunctionSpace(mesh, u_el, constrained_domain=pbc)
    s_el = df.FiniteElement("CG", mesh.ufl_cell(), 1)
    S = df.FunctionSpace(mesh, "CG", 1, constrained_domain=pbc)
    u = df.TrialFunction(V)
    v = df.TestFunction(V)
    u_ = df.Function(V, name="u")
    phi_ = df.Function(S, name="phi")
    phi = df.TrialFunction(S)
    q = df.TestFunction(S)
    noslip = df.DirichletBC(V, df.Constant((0., 0.)), subd, 1)

    beta = df.Constant(10.)
    F_u = df.dot(u - u0, v) * df.dx + beta * df.div(u) * df.div(v) * df.dx
    a_u, L_u = df.lhs(F_u), df.rhs(F_u)
    problem_u = df.LinearVariationalProblem(a_u, L_u, u_, bcs=noslip)
    solver_u = df.LinearVariationalSolver(problem_u)
    solver_u.parameters["linear_solver"] = "gmres"
    solver_u.parameters["preconditioner"] = "hypre_amg"
    solver_u.parameters["krylov_solver"]["relative_tolerance"] = 1e-9
    solver_u.parameters["krylov_solver"]["monitor_convergence"] = True

    if False:
        u_el = df.VectorElement("CG", mesh.ufl_cell(), 1)
        phi_el = df.FiniteElement("CG", mesh.ufl_cell(), 1)
        W = df.FunctionSpace(mesh, df.MixedElement([u_el, phi_el]), constrained_domain=pbc)
        w_ = df.Function(W)
        uw, phiw = df.TrialFunctions(W)
        vw, qw = df.TestFunctions(W)
        Fw = df.dot(uw - u0, vw) * df.dx - df.div(vw) * phiw * df.dx - df.div(uw) * qw * df.dx
        noslipw = df.DirichletBC(W.sub(0), df.Constant((0., 0.)), subd, 1)

        a = df.lhs(Fw)  # df.dot(df.grad(phi), df.grad(q)) * df.dx
        L = df.rhs(Fw)  # q * df.div(u_) * df.dx
        problem = df.LinearVariationalProblem(a, L, w_, bcs=noslipw)
        solver = df.LinearVariationalSolver(problem)
        solver.parameters["linear_solver"] = "gmres"
        solver.parameters["preconditioner"] = "hypre_amg"
        solver.parameters["krylov_solver"]["relative_tolerance"] = 1e-9
        solver.parameters["krylov_solver"]["monitor_convergence"] = True

    #A = df.assemble(a)
    #solver = df.KrylovSolver(A, "gmres")
    # Create vector that spans the null space
    #null_vec = df.Vector(w_.vector())
    #S.dofmap().set(null_vec, 1.0)
    #null_vec *= 1.0/null_vec.norm("l2")

    # Create null space basis object and attach to Krylov solver
    #null_space = df.VectorSpaceBasis([null_vec])
    #df.as_backend_type(A).set_nullspace(null_space)

    #isFluid.vector()[:] = ux  # is_fluid_xy.flatten()
    #with df.XDMFFile(mesh.mpi_comm(), os.path.join(analysisfolder, "isFluid.xdmf")) as xdmff:
    #    xdmff.write(isFluid)
    xdmff = df.XDMFFile(mesh.mpi_comm(), os.path.join(analysisfolder, "data.xdmf"))
    xdmff.parameters["rewrite_function_mesh"] = False
    xdmff.parameters["flush_output"] = True
    xdmff.parameters["functions_share_mesh"] = True

    mpi_print("Test")

    t_ = np.array([t for t, _ in timestamps])

    #istart = (rank * len(timestamps)) // size
    #istop = ((rank + 1) * len(timestamps)) // size
    #timestamps_block = list(enumerate(timestamps))[istart:istop]
    timestamps_block = list(enumerate(timestamps))

    comm.Barrier()
    if rank == 0:
        pbar = tqdm(total=len(timestamps_block))

    Sc = df.FunctionSpace(mesh, s_el)

    # psi stuff
    psi = df.TrialFunction(Sc)
    q = df.TestFunction(Sc)
    n = df.FacetNormal(mesh)

    vort = df.curl(u_)

    a = df.dot(df.grad(q), df.grad(psi)) * df.dx
    L = q * vort * df.dx + q * (n[1] * u_[0] - n[0] * u_[1]) * df.ds
    psi_ = df.Function(Sc, name="psi")
    bcs = [
        df.DirichletBC(
            Sc, df.Constant(0.),
            "on_boundary && sqrt(x[0]*x[0] + x[1]*x[1]) < DOLFIN_EPS_LARGE"
            )
    ]

    for it, timestamp in timestamps_block:
        t = timestamp[0]
        h5filename = timestamp[1]
        mpi_print(t, h5filename)
        with h5py.File(os.path.join(felbm_folder, h5filename), "r") as h5f:
            #p[:, :] = np.array(h5f["pressure"]).reshape((nz, ny, nx))[nz // 2, :, :]
            #rho[:, :] = np.array(h5f["density"]).reshape((nz, ny, nx))[nz // 2, :, :]
            ux[:, :] = np.array(h5f["u_x"]).reshape((nz, ny, nx))[nz // 2, :, :]
            uy[:, :] = np.array(h5f["u_y"]).reshape((nz, ny, nx))[nz // 2, :, :]
        ux0.vector().set_local(ux.flatten()[is_fluid][loc2glob])
        uy0.vector().set_local(uy.flatten()[is_fluid][loc2glob])

        # TODO: only works when there is a mass density difference between phases!!
        #rho_mean = rho[is_fluid_xy].mean()
        #phase1 = np.logical_and(rho < rho_mean, is_fluid_xy)
        #phase2 = np.logical_and(np.logical_not(phase1), is_fluid_xy)
        # u.assign(df.project(u0, V=V, bcs=noslip, solver_type="gmres", preconditioner_type="amg"))
        solver_u.solve()

        #        A = df.assemble(a)
        #    b = df.assemble(L)
        #[bc.apply(A, b) for bc in bcs]
        #df.solve(A, psi_.vector(), b, "gmres", "hypre_amg")

        #b = df.assemble(L)
        #null_space.orthogonalize(b)
        #solver.solve(w_.vector(), b)
        #solver.solve()
        #U_, Phi_ = w_.split(deepcopy=True)
        #U_.rename("U", "U")

        #df.solve(a == L, phi_, solver_parameters={"linear_solver": "gmres", "preconditioner": "hypre_amg"})
        #solver.solve(phi_.vector(), b)

        #xdmff.write(ux0, t)
        #xdmff.write(uy0, t)
        xdmff.write(u_, t)
        #xdmff.write(U_, t)
        #xdmff.write(psi_, t)

        if rank == 0:
            pbar.update(1)
        
    #comm.Reduce(uxt_loc, uxt, op=MPI.SUM, root=0)
    #comm.Reduce(uyt_loc, uyt, op=MPI.SUM, root=0)
    #comm.Reduce(S1t_loc, S1t, op=MPI.SUM, root=0)
    #comm.Reduce(S2t_loc, S2t, op=MPI.SUM, root=0)

    if rank == 0:
        #tdata = np.vstack((t_, uxt, uyt, S1t, S2t))
        #np.savetxt(os.path.join(analysisfolder, "tdata.dat"), tdata)
        pass