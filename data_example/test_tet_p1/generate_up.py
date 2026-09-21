import dolfin as df

DIRS = (0, 1, 2)   # the periodic directions dolfin_params.dat declares


class PBC(df.SubDomain):
    """The min faces of the periodic directions are the masters; each max face
    maps onto its image, so a node and its image share one dof."""
    def inside(self, x, on):
        return bool(on and any(df.near(x[k], 0) for k in DIRS)
                    and not any(df.near(x[k], 1) for k in DIRS))

    def map(self, x, y):
        for k in range(len(x)):
            y[k] = x[k] - 1 if k in DIRS and df.near(x[k], 1) else x[k]


mesh = df.UnitCubeMesh(10, 10, 10)
pbc = PBC()
V = df.VectorFunctionSpace(mesh, "CG", 1, constrained_domain=pbc)
P = df.FunctionSpace(mesh, "CG", 1, constrained_domain=pbc)

u = df.interpolate(
    df.Expression(("sin(2*M_PI*x[1])",
                   "sin(2*M_PI*x[2])",
                   "sin(2*M_PI*x[0])"), degree=1),
    V)
p = df.interpolate(
    df.Expression("sin(2*M_PI*(x[0]+x[1]+x[2]))", degree=2), P)

with df.HDF5File(mesh.mpi_comm(), "mesh.h5", "w") as h5f_mesh:
    h5f_mesh.write(mesh, "mesh")

with df.HDF5File(mesh.mpi_comm(), "up_0.h5", "w") as h5f_up:
    h5f_up.write(u, "u")
    h5f_up.write(p, "p")
