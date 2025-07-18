"""
Custom module for parametric and shape optimization
Authors : S. Gaydier, I. Zehavi, T. Cherrière, P. Gangl
ChatGPT helped writing the docstrings, which were edited and validated by the authors afterward.

functions list

- rot
- gen_mesh9
- gen_mesh2
- referenceVelocity
- create_plots
- update_plots
"""

import ngsolve as ngs
import numpy as np
import matplotlib.pyplot as plt
from IPython import display

## MAGNETIC CONSTANT
mu0 = 4e-7 * np.pi
nu0 = 1 / mu0

## OPERATORS
R = ngs.CoefficientFunction((0, 1, -1, 0), dims=(2, 2))  # Rotation matrix of -90°
Id = ngs.CoefficientFunction((1, 0, 0, 1), dims=(2, 2))  # Identity matrix


def rot(a):
    """
    Compute the rotational (curl) of a scalar or vector field in 2D.

    This function applies the predefined rotation matrix `R` (representing a -90° rotation)
    to the gradient of the input field `a`. The computation uses the `ngs.grad` function
    from NGSolve.

    Parameters
    ----------
    a : ngs.fem.CoefficientFunction
        Input field (scalar or vector) defined in NGSolve for which the rotational
        (curl) is to be computed.

    Returns
    -------
    ngs.fem.CoefficientFunction
        The result of applying the rotation matrix `R` to the gradient of the input field `a`.

    Notes
    -----
    - The rotation matrix `R` is defined as a global constant in the module:
      ```
      R = ngs.CoefficientFunction((0, 1, -1, 0), dims=(2, 2))  # -90° rotation matrix
      ```
    - This function assumes a 2D context, where the rotational is a perpendicular vector field.

    Examples
    --------
    >>> from ngsolve import x, y, grad, CoefficientFunction
    >>> a = x*y  # Example coefficient function
    >>> curl_a = rot(a)
    >>> print(curl_a)
    """
    return R * ngs.grad(a)


## GEOMETRY
from netgen.geom2d import SplineGeometry

a = 1e-2
ha = 1e-2
ba = 1e-2
e = 5e-3


def gen_mesh9(air_gap, maxh=2e-3):
    """
    Generate a 2D triangular mesh for a specific geometry with predefined boundary conditions and materials.

    This function creates a triangular mesh using `SplineGeometry` to define the geometry of the system.
    It includes points, lines, and material regions, and applies specific boundary conditions.
    The generated mesh is used for finite element analysis in NGSolve.

    Parameters
    ----------
    air_gap : float
        The size of the air gap in the geometry (in meters).
    maxh : float, optional
        The maximum element size for meshing (default is 2e-3 meters).

    Returns
    -------
    ngs.Mesh
        The generated 2D triangular mesh in NGSolve format.

    Geometry Details
    ----------------
    - The geometry consists of a rectangular core, coils, and surrounding air regions.
    - Points are defined to create the structure with specific mesh size controls.
    - Boundary conditions:
        * `front`: Applied to external air boundaries.
        * `coilVert` and `coilHor`: Applied to vertical and horizontal coil segments, respectively.
        * `segment1` and `segment2`: Define additional boundaries for inner regions.
        * `arc`: A curved arc boundary.

    Materials
    ---------
    - `air`: External air region.
    - `core`: Core material.
    - `coil`: Coil region.

    Notes
    -----
    - Mesh refinement is controlled by `maxh`, and finer subdivisions are applied to specific regions like the coils.
    - Splines and lines are used to accurately define the geometry.
    - This function relies on the `netgen.geom2d.SplineGeometry` and `ngsolve` libraries.

    Examples
    --------
    >>> mesh = gen_mesh9(air_gap=0.001, maxh=0.0005)
    >>> ngs.Draw(mesh)
    """
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
    """
    Generate a 2D triangular mesh for a geometry with variable air gaps and predefined boundary conditions.

    This function creates a triangular mesh using `SplineGeometry` to define a geometry with non-constant air gaps.
    Points, lines, and material regions are defined, along with specific boundary conditions.
    The generated mesh can be used for finite element simulations in NGSolve.

    Parameters
    ----------
    airgap : list of float
        A list of four values `[e11, e12, e21, e22]` representing the air gap sizes at different points of the geometry (in meters).
    maxh : float, optional
        The maximum element size for meshing (default is 2e-3 meters).

    Returns
    -------
    ngs.Mesh
        The generated 2D triangular mesh in NGSolve format.

    Geometry Details
    ----------------
    - The geometry is defined with adjustable air gap parameters: `e11`, `e12`, `e21`, and `e22`.
    - Points are defined with associated mesh refinement levels:
        * Fine mesh for coil areas.
        * Medium mesh for optimization regions.
        * Coarser mesh for outer regions.
    - Boundary conditions:
        * `airgap1`, `airgap2`: Represent air gap segments.
        * `optimVert`, `default`: Boundaries for vertical and horizontal optimization regions.
        * `domainHor`, `domainVert`: Boundaries for general domains.
        * `arc`: Represents a curved arc boundary.

    Materials
    ---------
    - `air`: External air region.
    - `core`: Core material.
    - `coil`: Coil region.

    Notes
    -----
    - This geometry is more flexible than `gen_mesh9`, allowing different air gap sizes for the left and right sections.
    - Mesh refinement is dynamically controlled by `maxh`, with finer subdivisions applied to specific areas.
    - The geometry is built using splines and lines for accurate definition.

    Examples
    --------
    >>> airgap = [0.001, 0.0015, 0.001, 0.0015]
    >>> mesh = gen_mesh2(airgap, maxh=0.0005)
    >>> ngs.Draw(mesh)
    """
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


def gen_meshN(airgap, maxh=2e-3):
    
    NControlPoints = len(airgap)
    r = 0.04
    maxhFine = maxh / 25
    maxhMed = maxh / 5
    geo = SplineGeometry()
    
    x_airgap = np.linspace(0, a/2, NControlPoints//2)
    pnt_airgap_1 = [(x_airgap[i], airgap[i]) for i in range(NControlPoints//2)]
    pnt_airgap_2 = [(x_airgap[i]+a/2+ba, airgap[i+NControlPoints//2]) for i in range(NControlPoints//2)]

    pnts = pnt_airgap_1 +[
        #(0, e11 / 2),  # p1
        #(a / 2, e12 / 2),  # p2
        (a / 2, e / 2 + ha / 2),  # p3
        (a / 2 + ba, e / 2 + ha / 2),  # p4
    ] + pnt_airgap_2 + [
        #(a / 2 + ba, e21 / 2),  # p5
        #(a + ba, e22 / 2),  # p6
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

    pointH = [maxhFine]*(NControlPoints//2) + [
        #maxh, 
        #maxhFine,
        maxhFine,
        maxhFine,] + [maxhFine]*(NControlPoints//2) + [
        #maxhFine,
        #maxhFine,
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

    pAirgap1 = [geo.AppendPoint(*pnts[i], pointH[i]) for i in range(0, NControlPoints//2)]
    (p3, p4) = [geo.AppendPoint(*pnts[i], pointH[i]) for i in range(NControlPoints//2,NControlPoints//2+2)]
    pAirgap2 = [geo.AppendPoint(*pnts[i], pointH[i]) for i in range(NControlPoints//2+2, NControlPoints+2)]

    (p7, p8, p001, p002, p003, p00, p01, p02, p03, p04, p05) = [
        geo.AppendPoint(*pnts[i], pointH[i]) for i in range( NControlPoints+2, len(pnts))
    ]

    # List of lines with boundary conditions and domains
    lines = [[["line", pAirgap1[i], pAirgap1[i+1]], {"bc": "airgap1", "leftdomain": 2, "rightdomain": 4, "maxh": maxhFine}] for i in range(NControlPoints//2-1)] +[
        
        #[["line", p1, p2], {"bc": "airgap1", "leftdomain": 2, "rightdomain": 4, "maxh": maxhFine}],
        #[["line", p2, p3], {"bc": "optimVert", "leftdomain": 2, "rightdomain": 3, "maxh": maxhMed}],
        [["line", pAirgap1[-1], p3], {"bc": "optimVert", "leftdomain": 2, "rightdomain": 3, "maxh": maxhMed}],
        [["line", p3, p4], {"bc": "default", "leftdomain": 2, "rightdomain": 3, "maxh": maxhMed}],
        [["line", p4, pAirgap2[0]], {"bc": "optimVert", "leftdomain": 2, "rightdomain": 3, "maxh": maxhMed}],
        #[["line", p4, p5], {"bc": "optimVert", "leftdomain": 2, "rightdomain": 3, "maxh": maxhMed}],
        ] + [[["line", pAirgap2[i], pAirgap2[i+1]], {"bc": "airgap2", "leftdomain": 2, "rightdomain": 1, "maxh": maxhFine}] for i in range(NControlPoints//2-1)] +[
        #[["line", p5, p6], {"bc": "airgap2", "leftdomain": 2, "rightdomain": 1, "maxh": maxhFine}],
        [["line", pAirgap2[-1], p002], {"bc": "front", "leftdomain": 2, "rightdomain": 1, "maxh": maxhMed}],
        [["line", p002, p7], {"bc": "front", "leftdomain": 2, "rightdomain": 1, "maxh": maxhMed}],
        [["line", p7, p8], {"bc": "front", "leftdomain": 2, "rightdomain": 1, "maxh": maxhMed}],
        [["line", p00, p03], {"bc": "domainHor", "leftdomain": 4, "rightdomain": 0, "maxh": maxhFine}],
        [["line", p03, p04], {"bc": "segment1", "leftdomain": 3, "rightdomain": 0, "maxh": maxhMed}],
        [["line", p04, p001], {"bc": "domainHor", "leftdomain": 1, "rightdomain": 0, "maxh": maxhFine}],
        [["line", p001, p01], {"bc": "segment1", "leftdomain": 1, "rightdomain": 0}],
        [["line", p02, p8], {"bc": "segment2", "leftdomain": 1, "rightdomain": 0}],
        [["line", p8, p003], {"bc": "domainVert", "leftdomain": 2, "rightdomain": 0, "maxh": maxhMed}],
        #[["line", p003, p1], {"bc": "domainVert", "leftdomain": 2, "rightdomain": 0, "maxh": maxhMed}],
        [["line", p003, pAirgap1[0]], {"bc": "domainVert", "leftdomain": 2, "rightdomain": 0, "maxh": maxhMed}],
        #[["line", p1, p00], {"bc": "domainVert", "leftdomain": 4, "rightdomain": 0, "maxh": maxhMed}],
        [["line", pAirgap1[0], p00], {"bc": "domainVert", "leftdomain": 4, "rightdomain": 0, "maxh": maxhMed}],
        #[["line", p04, p5], {"bc": "optimVert", "leftdomain": 3, "rightdomain": 1, "maxh": maxhMed}],
        [["line", p04, pAirgap2[0]], {"bc": "optimVert", "leftdomain": 3, "rightdomain": 1, "maxh": maxhMed}],
        #[["line", p2, p03], {"bc": "optimVert", "leftdomain": 3, "rightdomain": 4, "maxh": maxhMed}],
        [["line", pAirgap1[-1], p03], {"bc": "optimVert", "leftdomain": 3, "rightdomain": 4, "maxh": maxhMed}],
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


def referenceVelocity(mesh):
    """
    Generate reference velocity fields for the air gap regions of a mesh. This is useful to compute the parametric sensitivities
    by a scalar product with shape derivatives.

    This function defines vectorial velocity fields in the air gap regions of a mesh.
    These reference velocities are represented as `GridFunction` objects associated
    with `VectorH1` finite element spaces. The velocity fields are specifically defined
    on the boundaries of the air gap regions `airgap1` and `airgap2`.

    Parameters
    ----------
    mesh : ngsolve.Mesh
        The finite element mesh on which the velocity fields are defined.

    Returns
    -------
    tuple of ngs.GridFunction
        A tuple containing the following `GridFunction` objects:
        - `v11` : Reference velocity in the `airgap1` region, boundary defined by `(0, 1 - ngs.x / (a / 2))`.
        - `v12` : Reference velocity in the `airgap1` region, boundary defined by `(0, ngs.x / (a / 2))`.
        - `v21` : Reference velocity in the `airgap2` region, boundary defined by `(0, 1 - (ngs.x - (a / 2 + ba)) / (a / 2))`.
        - `v22` : Reference velocity in the `airgap2` region, boundary defined by `(0, (ngs.x - (a / 2 + ba)) / (a / 2))`.

    Notes
    -----
    - The function uses `VectorH1` finite element spaces to define the velocity fields.
    - Dirichlet boundary conditions are applied to ensure the velocity is defined only on the specified boundaries (`airgap1` and `airgap2`).
    - The velocity fields are computed as functions of the `x`-coordinate of the mesh nodes, scaled based on the geometric parameters `a` and `ba`.
    - `fes1` and `fes2` correspond to the `VectorH1` spaces for `airgap1` and `airgap2`, respectively.

    Examples
    --------
    >>> from ngsolve import Mesh
    >>> mesh = ngs.Mesh(...)
    >>> v11, v12, v21, v22 = referenceVelocity(mesh)
    >>> ngs.Draw(v11)
    >>> ngs.Draw(v12)
    """
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


## DISPLAY


def create_plots(num_plots):
    """
    Create a specified number of side-by-side plots with grid enabled.

    This function sets up a figure with `num_plots` subplots arranged in a single row,
    with grids enabled on each subplot. It also prepares a display handle for use in
    Jupyter notebooks.

    Parameters
    ----------
    num_plots : int
        The number of subplots to create in the figure.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created Matplotlib figure.
    axes : list of matplotlib.axes._subplots.AxesSubplot
        A list of subplot axes, one for each plot.
    display_handle : IPython.core.display.DisplayHandle
        A display handle that can be updated dynamically in Jupyter notebooks.

    Notes
    -----
    - The figure width is scaled proportionally to the number of plots (`10 * num_plots`),
      with a fixed height of 8.
    - If `num_plots` is 1, the function ensures `axes` is returned as a list for consistency.

    Examples
    --------
    >>> fig, axes, display_handle = create_plots(3)
    >>> for ax, data in zip(axes, [data1, data2, data3]):
    ...     ax.plot(data)
    ...     ax.set_title("Plot Title")
    >>> display_handle.update(fig)  # Dynamically update the display in Jupyter
    """
    cm = 1/2.54
    fig, axes = plt.subplots(1, num_plots, figsize=(8*cm * num_plots, 7*cm))
    plt.rcParams.update({'font.size': 10, 'lines.linewidth' : 0.5, 'lines.markersize': 2.0, 'grid.linewidth' : 0.3})
    if num_plots == 1:
        axes = [axes]
    for i, ax in enumerate(axes):
        ax.grid(True)
    fig.tight_layout()
    return fig, axes, display.display("", display_id=True)


def update_plots(fig, axes, hdisplay, data, curve_labels, semilogy=False):
    """
    Update the content of multiple subplots with new data and refresh the display.

    This function clears the existing content of each subplot, plots new data, and updates
    the display dynamically (e.g., in a Jupyter notebook). It supports both linear and
    logarithmic scaling for the y-axis.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure containing the subplots to update.
    axes : list of matplotlib.axes._subplots.AxesSubplot
        A list of subplot axes to be updated.
    hdisplay : IPython.core.display.DisplayHandle
        A display handle used to dynamically update the output in Jupyter notebooks.
    data : list of list of list of float
        A nested list of numerical data to plot.
        `data[i][j]` represents the j-th curve for the i-th subplot.
    curve_labels : list of list of str
        A nested list of labels for each curve in the subplots.
        `curve_labels[i][j]` is the label for the j-th curve in the i-th subplot.
    semilogy : bool, optional
        If True, use a logarithmic scale for the y-axis. Defaults to False (linear scale).

    Returns
    -------
    None
        The function modifies the plots and display in-place.

    Notes
    -----
    - The function assumes that `data` and `curve_labels` have matching shapes,
      with the same number of curves for each subplot.
    - The x-axis limits are set dynamically to match the length of the curves.
    - The y-axis limits are adjusted based on the minimum and maximum values of the data.

    Examples
    --------
    >>> fig, axes, hdisplay = create_plots(2)
    >>> data = [[[1, 2, 3], [3, 2, 1]], [[4, 5, 6], [6, 5, 4]]]  # Two subplots, each with two curves
    >>> curve_labels = [["Curve A", "Curve B"], ["Curve C", "Curve D"]]
    >>> update_plots(fig, axes, hdisplay, data, curve_labels, semilogy=False)
    >>> # Dynamically update in Jupyter after modifying data
    >>> data[0][0].append(4)  # Append new data point
    >>> update_plots(fig, axes, hdisplay, data, curve_labels)
    """
    for i, ax in enumerate(axes):
        ax.clear()  # Clear the current axis
        for j, (curve, label) in enumerate(zip(data[i], curve_labels[i])):
            ax.plot(np.arange(len(curve)), curve, label=label, marker="o", linestyle="-")

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
