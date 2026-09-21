"""A tetrahedral case for mode=tetfreq: a velocity given as Fourier modes in time.

The loader reads u(x, t) = sum_k a_k cos(omega0 (k t + t_k)) u_k(x), with
omega0 = 2 pi / tau, one dolfin checkpoint per mode k and one line per mode in
freqstamps.dat holding t_k, a_k and the file name, in order of k. t_k is a
phase in time units and must be 0 on the mean (k = 0). A file may instead give
each line as `omega_k phi_k a_k file`, a_k cos(omega_k t + phi_k) u_k(x), with
no tau, in any order and with a frequency on several lines: a travelling
Fourier mode A cos(omega t) + B sin(omega t) is two lines, phi = 0 for A and
phi = -pi/2 for B.

The modes here are polynomials of degree two, so P2 holds them exactly and a
test can compare the probed velocity and its gradient with the closed form to
round-off. Nothing depends on x, the one direction declared periodic; u_y and
u_z vanish on the y and z faces, so a tracer started inside stays inside.
"""
import dolfin as df

DIRS = (0,)   # the periodic directions dolfin_params.dat declares

TAU = 1.0     # base period: omega0 = 2 pi / tau

# mode k -> (time shift t_k, amplitude a_k, velocity, pressure)
MODES = [
    (0.0, 1.0, ("1 + x[1]*x[2]", "0.2*x[1]*(1-x[1])", "0.0"), "x[1] + 2*x[2]"),
    (0.125, 0.7, ("x[2]*(1-x[2])", "0.0", "0.3*x[2]*(1-x[2])"), "x[2]"),
    (-0.25, 0.4, ("x[1] - 0.5", "0.1*x[1]*(1-x[1])", "-0.2*x[2]*(1-x[2])"), "x[1]"),
]


class PBC(df.SubDomain):
    """The min faces of the periodic directions are the masters; each max face
    maps onto its image, so a node and its image share one dof."""
    def inside(self, x, on):
        return bool(on and any(df.near(x[k], 0) for k in DIRS)
                    and not any(df.near(x[k], 1) for k in DIRS))

    def map(self, x, y):
        for k in range(len(x)):
            y[k] = x[k] - 1 if k in DIRS and df.near(x[k], 1) else x[k]


mesh = df.UnitCubeMesh(6, 6, 6)
pbc = PBC()
V = df.VectorFunctionSpace(mesh, "CG", 2, constrained_domain=pbc)
P = df.FunctionSpace(mesh, "CG", 1, constrained_domain=pbc)

with df.HDF5File(mesh.mpi_comm(), "mesh.h5", "w") as h5f_mesh:
    h5f_mesh.write(mesh, "mesh")

stamps = []
for k, (t_k, a_k, u_expr, p_expr) in enumerate(MODES):
    u_ = df.interpolate(df.Expression(u_expr, degree=2), V)
    p_ = df.interpolate(df.Expression(p_expr, degree=1), P)
    fname = f"up_{k}.h5"
    with df.HDF5File(mesh.mpi_comm(), fname, "w") as h5f_up:
        h5f_up.write(u_, "u")
        h5f_up.write(p_, "p")
    stamps.append(f"{t_k} {a_k} {fname}")

with open("freqstamps.dat", "w") as ofile:
    ofile.write("\n".join(stamps) + "\n")

dfprms = f"""velocity_space=P2
pressure_space=P1
freqstamps=freqstamps.dat
mesh=mesh.h5
periodic_x={"true" if 0 in DIRS else "false"}
periodic_y={"true" if 1 in DIRS else "false"}
periodic_z={"true" if 2 in DIRS else "false"}
rho=1.0
tau={TAU}
t_max=1e8
t_min=0"""

with open("dolfin_params.dat", "w") as ofile:
    ofile.write(dfprms)
