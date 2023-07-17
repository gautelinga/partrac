import dolfin as df

mesh = df.UnitSquareMesh(10, 10)
V = df.VectorFunctionSpace(mesh, "CG", 2)
P = df.FunctionSpace(mesh, "CG", 1)

u0_ = df.interpolate(df.Expression(("sin(2*M_PI*x[1])", "0.0"), degree=2), V)
u1_ = df.interpolate(df.Expression(("0.0", "sin(2*M_PI*x[0])"), degree=2), V)
p0_ = df.interpolate(df.Expression("0", degree=2), P)

with df.HDF5File(mesh.mpi_comm(), "mesh.h5", "w") as h5f_mesh:
    h5f_mesh.write(mesh, "mesh")

with df.HDF5File(mesh.mpi_comm(), "up_0.h5", "w") as h5f_up:
    h5f_up.write(u0_, "u")
    h5f_up.write(p0_, "p")

with df.HDF5File(mesh.mpi_comm(), "up_1.h5", "w") as h5f_up:
    h5f_up.write(u1_, "u")
    h5f_up.write(p0_, "p")