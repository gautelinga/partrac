import dolfin as df

class PBC(df.SubDomain):
    def inside(self, x, on_boundary):
        return bool((df.near(x[0], 0) or
                     df.near(x[1], 0) or
                     df.near(x[2], 0)) and
                    (not ((df.near(x[0], 1) and df.near(x[2], 1)) or
                          (df.near(x[0], 1) and df.near(x[1], 1)) or
                          (df.near(x[1], 1) and df.near(x[2], 1))))
                    and on_boundary)
    def map(self, x, y):
        if df.near(x[0], 1) and df.near(x[1], 1) and df.near(x[2], 1):
            y[0] = x[0] - 1
            y[1] = x[1] - 1
            y[2] = x[2] - 1
        if df.near(x[0], 1) and df.near(x[2], 1):
            y[0] = x[0] - 1
            y[1] = x[1]
            y[2] = x[2] - 1
        elif df.near(x[1], 1) and df.near(x[2], 1):
            y[0] = x[0]
            y[1] = x[1] - 1
            y[2] = x[2] - 1
        elif df.near(x[0], 1) and df.near(x[1], 1):
            y[0] = x[0] - 1
            y[1] = x[1] - 1
            y[2] = x[2]
        elif df.near(x[0], 1):
            y[0] = x[0] - 1
            y[1] = x[1]
            y[2] = x[2]
        elif df.near(x[1], 1):
            y[0] = x[0]
            y[1] = x[1] - 1
            y[2] = x[2]
        elif df.near(x[2], 1):
            y[0] = x[0]
            y[1] = x[1]
            y[2] = x[2] - 1


mesh = df.UnitCubeMesh(10, 10, 10)
pbc = PBC()
V = df.VectorFunctionSpace(mesh, "CG", 2, constrained_domain=pbc)
P = df.FunctionSpace(mesh, "CG", 2, constrained_domain=pbc)

u = df.interpolate(
    df.Expression(("sin(2*M_PI*x[1])",
                   "sin(2*M_PI*x[2])",
                   "sin(2*M_PI*x[0])"), degree=2),
    V)
p = df.interpolate(
    df.Expression("sin(2*M_PI*(x[0]+x[1]+x[2]))", degree=2), P)

with df.HDF5File(mesh.mpi_comm(), "mesh.h5", "w") as h5f_mesh:
    h5f_mesh.write(mesh, "mesh")

with df.HDF5File(mesh.mpi_comm(), "up_0.h5", "w") as h5f_up:
    h5f_up.write(u, "u")
    h5f_up.write(p, "p")
