import dolfin as df
import argparse

parser = argparse.ArgumentParser(description="Make plane Poiseuille flow")
parser.add_argument("-dim", default=0, type=int, help="Dimension of flow")
args = parser.parse_args()

mesh = df.UnitSquareMesh(10, 10)
V = df.VectorFunctionSpace(mesh, "CG", 2)
P = df.FunctionSpace(mesh, "CG", 1)

tdim = 1 if args.dim == 0 else 0
u_expr = ["0.0", "0.0"]
u_expr[args.dim] = f"6.*x[{tdim}]*(1-x[{tdim}])"

u_ = df.interpolate(df.Expression(u_expr, degree=2), V)
p_ = df.interpolate(df.Expression("0", degree=2), P)

with df.HDF5File(mesh.mpi_comm(), "mesh.h5", "w") as h5f_mesh:
    h5f_mesh.write(mesh, "mesh")

with df.HDF5File(mesh.mpi_comm(), "up_0.h5", "w") as h5f_up:
    h5f_up.write(u_, "u")
    h5f_up.write(p_, "p")

dfprms = f"""velocity_space=P2
pressure_space=P1
timestamps=timestamps.dat
mesh=mesh.h5
periodic_x={"true" if args.dim == 0 else "false"}
periodic_y={"true" if args.dim == 1 else "false"}
periodic_z=false
rho=1.0"""

with open("dolfin_params.dat", "w") as ofile:
    ofile.write(dfprms)