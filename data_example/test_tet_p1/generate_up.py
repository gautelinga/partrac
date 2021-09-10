import dolfin as df

mesh = df.UnitCubeMesh(10, 10, 10)
V = df.VectorFunctionSpace(mesh, "CG", 1)
P = df.FunctionSpace(mesh, "CG", 1)

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
