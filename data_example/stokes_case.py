"""The solve the solved OpenFOAM fixtures share (openfoam_cylinder2d,
openfoam_sphere3d): a Stokes-like steady flow (nu = 1) under a body force of
1 at 20 degrees to x in the x-y plane, the box cyclic all round, one wall
patch no-slip, converged to residuals 1e-9. Each fixture's script writes its
mesh dictionaries and calls these:

  dicts(wall)               write system/'s solver dictionaries, constant/ and 0/
  measure(time, dim, depth) print the at-rest quantities at time (ascii U and Vc)
"""
import os
import re

import numpy as np

ANGLE = 20.0
FORCE = (np.cos(np.radians(ANGLE)), np.sin(np.radians(ANGLE)), 0.0)
HEADER = """FoamFile
{
    format      ascii;
    class       %s;
    location    "%s";
    object      %s;
}
"""


def write(rel, body, cls="dictionary"):
    os.makedirs(os.path.dirname(rel), exist_ok=True)
    loc, obj = os.path.split(rel)
    with open(rel, "w") as f:
        f.write(HEADER % (cls, loc, obj) + "\n" + body)


CONTROL = """application     foamRun;
solver          incompressibleFluid;
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         50000;
deltaT          1;
writeControl    timeStep;
writeInterval   50000;
purgeWrite      1;
writeFormat     ascii;
writePrecision  17;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable false;
"""

SCHEMES = """ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linear;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
"""

SOLUTION = """solvers
{
    p
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-14;
        relTol          0.01;
    }
    U
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-14;
        relTol          0.01;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 2;
    consistent      yes;
    pRefCell        0;
    pRefValue       0;
    residualControl
    {
        p               1e-09;
        U               1e-09;
    }
}

relaxationFactors
{
    equations
    {
        U               0.9;
        ".*"            0.9;
    }
}
"""

BODY_FORCE = """momentumSource
{
    type            semiImplicitSource;
    select          all;
    volumeMode      specific;
    sources
    {
        U
        {
            explicit    (%r %r %r);
            implicit    0;
        }
    }
}
"""

U0 = """dimensions      [0 1 -1 0 0 0 0];
internalField   uniform (0 0 0);
boundaryField
{
    %s { type noSlip; }
    #includeEtc "caseDicts/setConstraintTypes"
}
"""

P0 = """dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 0;
boundaryField
{
    %s { type zeroGradient; }
    #includeEtc "caseDicts/setConstraintTypes"
}
"""


def dicts(wall):
    write("system/controlDict", CONTROL)
    write("system/fvSchemes", SCHEMES)
    write("system/fvSolution", SOLUTION)
    write("constant/physicalProperties", "viscosityModel  constant;\nnu              [0 2 -1 0 0 0 0] 1;\n")
    write("constant/momentumTransport", "simulationType laminar;\n")
    write("constant/fvModels", BODY_FORCE % tuple(float(v) for v in FORCE))
    write("0/U", U0 % wall, "volVectorField")
    write("0/p", P0 % wall, "volScalarField")


def body(path):
    text = open(path).read()
    return re.sub(r"//[^\n]*|/\*.*?\*/", " ", text[text.index("}", text.index("FoamFile")) + 1:], flags=re.S)


def read_list(path, width):
    """The first list of an ascii mesh or field file, n x width."""
    b = body(path)
    m = re.search(r"(\d+)\s*\(", b)
    n = int(m.group(1))
    v = np.array(b[m.end():].replace("(", " ").replace(")", " ").split()[:width * n], dtype=float)
    return v.reshape(n, width)


def measure(time, dim, depth=1.0):
    """Weighted by the cells' volumes, in 2D their areas (volume over the
    depth): the mean speed along the force, the largest speed, the mean cell
    size (mean volume^(1/3), area^(1/2)), dt = one mean cell at the fastest
    cell, the steps of five throughflows, the at-rest bound."""
    U = read_list(time + "/U", 3)
    V = read_list(time + "/Vc", 1)[:, 0] / depth
    along = float((V * (U @ np.array(FORCE))).sum() / V.sum())
    umax = float(np.linalg.norm(U, axis=1).max())
    h = float(V.mean() ** (1 / dim))
    dt = h / umax
    print("time %s cells %d max|u| %.17g mean u.f %.17g h %.17g dt %.17g steps %d rest %.17g"
          % (time, len(V), umax, along, h, dt, int(np.ceil(5 / along / dt)), 1e-3 * along))
