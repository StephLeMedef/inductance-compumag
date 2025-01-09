import ngsolve as ngs
import numpy as np
from netgen.geom2d import CSG2d, Rectangle
from netgen.geom2d import EdgeInfo as EI, PointInfo as PI, Solid2d
import matplotlib.pyplot as plt
from IPython import display

R = ngs.CoefficientFunction((0, 1, -1, 0), dims=(2, 2))  # Rotation matrix of -90°
Id = ngs.CoefficientFunction((1, 0, 0, 1), dims=(2, 2))  # Identity matrix

# Magnetostatic constants
mu0 = 4e-7 * np.pi
nu0 = 1 / mu0


def rot(a):
    return R * ngs.grad(a)


def create_plots(num_plots):
    fig, axes = plt.subplots(1, num_plots, figsize=(10 * num_plots, 8))
    if num_plots == 1:
        axes = [axes]
    for i, ax in enumerate(axes):
        ax.grid(True)
    fig.tight_layout()
    return fig, axes, display.display("", display_id=True)


def update_plots(fig, axes, hdisplay, data, curve_labels, semilogy=False):
    for i, ax in enumerate(axes):
        ax.clear()  # Clear the current axis
        for j, (curve, label) in enumerate(zip(data[i], curve_labels[i])):
            ax.plot(np.arange(len(curve)), curve, label=label, marker="o", linestyle="-", markersize=6)

        # Dynamically set limits based on data
        ax.set_xlim(0, len(curve) - 1)
        min_y = min([min(c) for c in data[i]])
        max_y = max([max(c) for c in data[i]])
        ax.set_ylim(min_y, max_y)

        if semilogy:
            ax.set_yscale("log")
        else:
            ax.set_yscale("linear")

        ax.grid(True)
        ax.legend(loc="best")

    fig.tight_layout()
    hdisplay.update(fig)  # Update the display with the new plot


## GEOMETRY

from netgen.geom2d import SplineGeometry

# Geometry definition

a = 1e-2
ha = 1e-2
ba = 1e-2
e = 5e-3


def referenceVelocity(mesh):
    fes1 = ngs.VectorH1(mesh, dirichlet="airgap1")
    v11 = ngs.GridFunction(fes1)
    v11.Set(ngs.CF((0, 1 - ngs.x / (a / 2))), ngs.BND)
    v12 = ngs.GridFunction(fes1)
    v12.Set(ngs.CF((0, ngs.x / (a / 2))), ngs.BND)
    fes2 = ngs.VectorH1(mesh, dirichlet="airgap2")
    v21 = ngs.GridFunction(fes2)
    v21.Set(ngs.CF((0, 1 - (ngs.x - (a / 2 + ba)) / (a / 2))), ngs.BND)
    v22 = ngs.GridFunction(fes2)
    v22.Set(ngs.CF((0, (ngs.x - (a / 2 + ba)) / (a / 2))), ngs.BND)
    return v11, v12, v21, v22


def gen_mesh9(air_gap, maxh=2e-3):
    """Gives a triangular mesh"""
    r = 0.04

    maxhFine = maxh / 25
    maxhMed = maxh / 5
    geo = SplineGeometry()
    pnts = [
        (0, air_gap / 2),  # p1
        (a / 2, air_gap / 2),  # p2
        (a / 2, e / 2 + ha / 2),  # p3
        (a / 2 + ba, e / 2 + ha / 2),  # p4
        (a / 2 + ba, air_gap / 2),  # p5
        (a + ba, air_gap / 2),  # p6
        (a + ba, e / 2 + ha / 2 + a / 2),  # p7
        (0, e / 2 + ha / 2 + a / 2),  # p8
        (a + ba, 0),  # p001
        (a + ba, e / 2 + ha / 2),  # p002
        (0, e / 2 + ha / 2),  # p003
        (0, 0),  # p00
        (r, 0),  # p01
        (0, r),  # p02
        (a / 2, 0),  # p03
        (a / 2 + ba, 0),  # p04
        (r, r),  # p05
    ]

    pointH = [
        maxh,
        maxhFine,
        maxhFine,
        maxhFine,
        maxhFine,
        maxhFine,
        maxh,
        maxh,
        maxh,
        maxh,
        maxh,
        maxh,
        maxh,
        maxh,
        maxhFine,
        maxhFine,
        maxh,
    ]

    (p1, p2, p3, p4, p5, p6, p7, p8, p001, p002, p003, p00, p01, p02, p03, p04, p05) = [
        geo.AppendPoint(*pnts[i], pointH[i]) for i in range(len(pnts))
    ]

    # List of lines with boundary conditions and domains
    lines = [
        [["line", p1, p2], {"bc": "front", "leftdomain": 2, "rightdomain": 4, "maxh": maxhFine}],
        [["line", p2, p3], {"bc": "coilVert", "leftdomain": 2, "rightdomain": 3, "maxh": maxhMed}],
        [["line", p3, p4], {"bc": "coilHor", "leftdomain": 2, "rightdomain": 3, "maxh": maxhMed}],
        [["line", p4, p5], {"bc": "coilVert", "leftdomain": 2, "rightdomain": 3, "maxh": maxhMed}],
        [["line", p5, p6], {"bc": "front", "leftdomain": 2, "rightdomain": 1, "maxh": maxhFine}],
        [["line", p6, p002], {"bc": "front", "leftdomain": 2, "rightdomain": 1, "maxh": maxhMed}],
        [["line", p002, p7], {"bc": "front", "leftdomain": 2, "rightdomain": 1, "maxh": maxhMed}],
        [["line", p7, p8], {"bc": "front", "leftdomain": 2, "rightdomain": 1, "maxh": maxhMed}],
        [["line", p00, p03], {"bc": "segment1", "leftdomain": 4, "rightdomain": 0, "maxh": maxhFine}],
        [["line", p03, p04], {"bc": "segment1", "leftdomain": 3, "rightdomain": 0, "maxh": maxhMed}],
        [["line", p04, p001], {"bc": "segment1", "leftdomain": 1, "rightdomain": 0, "maxh": maxhFine}],
        [["line", p001, p01], {"bc": "segment1", "leftdomain": 1, "rightdomain": 0}],
        [["line", p02, p8], {"bc": "segment2", "leftdomain": 1, "rightdomain": 0}],
        [["line", p8, p003], {"bc": "segment2", "leftdomain": 2, "rightdomain": 0, "maxh": maxhMed}],
        [["line", p003, p1], {"bc": "segment2", "leftdomain": 2, "rightdomain": 0, "maxh": maxhMed}],
        [["line", p1, p00], {"bc": "segment2", "leftdomain": 4, "rightdomain": 0, "maxh": maxhMed}],
        [["line", p04, p5], {"bc": "coilVert", "leftdomain": 3, "rightdomain": 1, "maxh": maxhMed}],
        [["line", p2, p03], {"bc": "coilVert", "leftdomain": 3, "rightdomain": 4, "maxh": maxhMed}],
        [["spline3", p01, p05, p02], {"bc": "arc", "leftdomain": 1, "rightdomain": 0}],
    ]

    # Append all lines to the geometry
    for line, props in lines:
        geo.Append(line, **props)

    # Set materials and meshing parameters
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "core")
    geo.SetMaterial(3, "coil")
    geo.SetMaterial(4, "air")
    ngmesh = geo.GenerateMesh(maxh=maxh)
    return ngs.Mesh(ngmesh)


def gen_mesh2(airgap, maxh=2e-3):
    """Gives a triangular mesh"""
    e11, e12, e21, e22 = airgap[0], airgap[1], airgap[2], airgap[3]
    r = 0.04
    maxhFine = maxh / 25
    maxhMed = maxh / 5
    geo = SplineGeometry()
    pnts = [
        (0, e11 / 2),  # p1
        (a / 2, e12 / 2),  # p2
        (a / 2, e / 2 + ha / 2),  # p3
        (a / 2 + ba, e / 2 + ha / 2),  # p4
        (a / 2 + ba, e21 / 2),  # p5
        (a + ba, e22 / 2),  # p6
        (a + ba, e / 2 + ha / 2 + a / 2),  # p7
        (0, e / 2 + ha / 2 + a / 2),  # p8
        (a + ba, 0),  # p001
        (a + ba, e / 2 + ha / 2),  # p002
        (0, e / 2 + ha / 2),  # p003
        (0, 0),  # p00
        (r, 0),  # p01
        (0, r),  # p02
        (a / 2, 0),  # p03
        (a / 2 + ba, 0),  # p04
        (r, r),  # p05
    ]

    pointH = [
        maxh,
        maxhFine,
        maxhFine,
        maxhFine,
        maxhFine,
        maxhFine,
        maxh,
        maxh,
        maxh,
        maxh,
        maxh,
        maxh,
        maxh,
        maxh,
        maxhFine,
        maxhFine,
        maxh,
    ]

    (p1, p2, p3, p4, p5, p6, p7, p8, p001, p002, p003, p00, p01, p02, p03, p04, p05) = [
        geo.AppendPoint(*pnts[i], pointH[i]) for i in range(len(pnts))
    ]

    # List of lines with boundary conditions and domains
    lines = [
        [["line", p1, p2], {"bc": "airgap1", "leftdomain": 2, "rightdomain": 4, "maxh": maxhFine}],
        [["line", p2, p3], {"bc": "optimVert", "leftdomain": 2, "rightdomain": 3, "maxh": maxhMed}],
        [["line", p3, p4], {"bc": "default", "leftdomain": 2, "rightdomain": 3, "maxh": maxhMed}],
        [["line", p4, p5], {"bc": "optimVert", "leftdomain": 2, "rightdomain": 3, "maxh": maxhMed}],
        [["line", p5, p6], {"bc": "airgap2", "leftdomain": 2, "rightdomain": 1, "maxh": maxhFine}],
        [["line", p6, p002], {"bc": "front", "leftdomain": 2, "rightdomain": 1, "maxh": maxhMed}],
        [["line", p002, p7], {"bc": "front", "leftdomain": 2, "rightdomain": 1, "maxh": maxhMed}],
        [["line", p7, p8], {"bc": "front", "leftdomain": 2, "rightdomain": 1, "maxh": maxhMed}],
        [["line", p00, p03], {"bc": "domainHor", "leftdomain": 4, "rightdomain": 0, "maxh": maxhFine}],
        [["line", p03, p04], {"bc": "segment1", "leftdomain": 3, "rightdomain": 0, "maxh": maxhMed}],
        [["line", p04, p001], {"bc": "domainHor", "leftdomain": 1, "rightdomain": 0, "maxh": maxhFine}],
        [["line", p001, p01], {"bc": "segment1", "leftdomain": 1, "rightdomain": 0}],
        [["line", p02, p8], {"bc": "segment2", "leftdomain": 1, "rightdomain": 0}],
        [["line", p8, p003], {"bc": "domainVert", "leftdomain": 2, "rightdomain": 0, "maxh": maxhMed}],
        [["line", p003, p1], {"bc": "domainVert", "leftdomain": 2, "rightdomain": 0, "maxh": maxhMed}],
        [["line", p1, p00], {"bc": "domainVert", "leftdomain": 4, "rightdomain": 0, "maxh": maxhMed}],
        [["line", p04, p5], {"bc": "optimVert", "leftdomain": 3, "rightdomain": 1, "maxh": maxhMed}],
        [["line", p2, p03], {"bc": "optimVert", "leftdomain": 3, "rightdomain": 4, "maxh": maxhMed}],
        [["spline3", p01, p05, p02], {"bc": "arc", "leftdomain": 1, "rightdomain": 0}],
    ]

    # Append all lines to the geometry
    for line, props in lines:
        geo.Append(line, **props)

    # Set materials and meshing parameters
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "core")
    geo.SetMaterial(3, "coil")
    geo.SetMaterial(4, "air")
    ngmesh = geo.GenerateMesh(maxh=maxh)
    return ngs.Mesh(ngmesh)
