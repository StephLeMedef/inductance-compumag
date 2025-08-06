"""
Custom module for parametric and shape optimization
Authors : S. Gaydier, I. Zehavi, T. Cherrière, P. Gangl
ChatGPT helped writing the docstrings, which were edited and validated by the authors afterward.

functions list

- rot
- gen_mesh9
- gen_mesh2
- gen_meshN
- referenceVelocity
- create_plots
- update_plots
"""

import ngsolve as ngs
import numpy as np
import matplotlib.pyplot as plt
from IPython import display

##################################################################################################################################
## MAGNETIC CONSTANT
mu0 = 4e-7 * np.pi
nu0 = 1 / mu0

##################################################################################################################################
## OPERATORS
R = ngs.CoefficientFunction((0, 1, -1, 0), dims=(2, 2))  # Rotation matrix of -90°
Id = ngs.CoefficientFunction((1, 0, 0, 1), dims=(2, 2))  # Identity matrix
Curl = lambda u : R * ngs.grad(u)                        # 2D vector curl operator

##################################################################################################################################
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
    pnt_airgap_1 = [(x_airgap[i], airgap[i]/2) for i in range(NControlPoints//2)]
    pnt_airgap_2 = [(x_airgap[i]+a/2+ba, airgap[i+NControlPoints//2]/2) for i in range(NControlPoints//2)]

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
        maxhFine,
        maxhFine,] + [maxhFine]*(NControlPoints//2) + [
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
        [["line", pAirgap1[-1], p3], {"bc": "optimVert", "leftdomain": 2, "rightdomain": 3, "maxh": maxhMed}],
        [["line", p3, p4], {"bc": "default", "leftdomain": 2, "rightdomain": 3, "maxh": maxhMed}],
        [["line", p4, pAirgap2[0]], {"bc": "optimVert", "leftdomain": 2, "rightdomain": 3, "maxh": maxhMed}],
        ] + [[["line", pAirgap2[i], pAirgap2[i+1]], {"bc": "airgap2", "leftdomain": 2, "rightdomain": 1, "maxh": maxhFine}] for i in range(NControlPoints//2-1)] +[
        [["line", pAirgap2[-1], p002], {"bc": "front", "leftdomain": 2, "rightdomain": 1, "maxh": maxhMed}],
        [["line", p002, p7], {"bc": "front", "leftdomain": 2, "rightdomain": 1, "maxh": maxhMed}],
        [["line", p7, p8], {"bc": "front", "leftdomain": 2, "rightdomain": 1, "maxh": maxhMed}],
        [["line", p00, p03], {"bc": "domainHor", "leftdomain": 4, "rightdomain": 0, "maxh": maxhFine}],
        [["line", p03, p04], {"bc": "segment1", "leftdomain": 3, "rightdomain": 0, "maxh": maxhMed}],
        [["line", p04, p001], {"bc": "domainHor", "leftdomain": 1, "rightdomain": 0, "maxh": maxhFine}],
        [["line", p001, p01], {"bc": "segment1", "leftdomain": 1, "rightdomain": 0}],
        [["line", p02, p8], {"bc": "domainVert", "leftdomain": 1, "rightdomain": 0}],
        [["line", p8, p003], {"bc": "domainVert", "leftdomain": 2, "rightdomain": 0, "maxh": maxhMed}],
        [["line", p003, pAirgap1[0]], {"bc": "domainVert", "leftdomain": 2, "rightdomain": 0, "maxh": maxhMed}],
        [["line", pAirgap1[0], p00], {"bc": "domainVert", "leftdomain": 4, "rightdomain": 0, "maxh": maxhMed}],
        [["line", p04, pAirgap2[0]], {"bc": "optimVert", "leftdomain": 3, "rightdomain": 1, "maxh": maxhMed}],
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

    return ngs.Mesh(ngmesh), pnt_airgap_1, pnt_airgap_2


def referenceVelocity(mesh : ngs.Mesh,
                      xPoints : list[float] | np.ndarray
                      ) -> list[ngs.GridFunction]:
    """
    Compute a list of elementary node displacements (also called velocities) 
    at the airgap boundary associated with unit vertical displacements of 
    control points located along the x-axis.

    Parameters
    ----------
    mesh : ngsolve.Mesh
        The finite element mesh of the geometry, which must include boundaries 
        named "airgap1" and "airgap2" (used as reference boundaries for the velocity field).
    
    xPoints : list or array-like of float
        The x-coordinates of the control points on the airgap boundary where 
        elementary displacements are defined. Assumes uniform spacing.

    Returns
    -------
    yDisplacementList : list of ngsolve.GridFunction
        A list of vector-valued GridFunctions representing the velocity fields 
        (or displacement) corresponding to unit vertical displacements 
        at each control point in `xPoints`. Each GridFunction has nonzero 
        displacement in the y-direction (second component) near the associated 
        control point and zero elsewhere.

    Note
    ----
    - Shape functions are 1 at the associated control point, zero at the neighboring 
    control points and linear in between
    """
    
    dx = xPoints[1]-xPoints[0] # assume regular spacing
    fesVelocity = ngs.VectorH1(mesh, dirichlet="airgap1|airgap2")
    yDisplacementList = []

    def regSign(x, p = 1e10):
        """ regularized sign function """
        return 2 * ngs.atan(p*x)/ngs.pi
    
    def xdist(x,Xref, dx = dx):
        """ normalized x-distance to a reference point  """
        return ngs.sqrt((x- Xref)**2)/dx

    def nodeShapeFunction(x, Xref, dx = dx):
        """ 1D shape function in the x-direction associated to a control point """
        sd = xdist(x, Xref, dx = dx)
        xi = (regSign(1-sd) +1 )/ 2
        return xi * (1-sd)
    
    for i in range(len(xPoints)):
        xRef = xPoints[i]
        yDisplacementList.append(ngs.GridFunction(fesVelocity))
        yDisplacementList[-1].Set(ngs.CF((0,nodeShapeFunction(ngs.x, xRef, dx = dx))), ngs.BND)

    return yDisplacementList


##################################################################################################################################
## OPTIMIZATION

def gradient_descent(fun : callable,               # function to minimize
                     dfun : callable,              # gradient
                     x : np.array,                  # initial value
                     xlb : float | list = -np.inf, # lower bound(s)
                     xub : float | list = np.inf,  # upper bound(s)
                     step : float = 1.,            # step size     
                     precond : bool = True,        # preconditioning of descent direction (recommended)
                     # convergence
                     iter_max : int = 50,          # max number of iterations
                     rstep_min : float = 1e-5,     # minimum step size relative to step
                     tol : float = 1e-6,           # absolute tolerance on sqrt(|grad * descent|)
                     # step update
                     coeff_armijo : float = 0.1,   # Armijo's coefficient
                     step_increase : float = 1.2,  # increase factor of the step when accepted
                     step_decrease : float = 0.5,   # decrease factor of the step when rejected
                     # inspect
                     verbosity : int = 1           # verbosity level (0 = silent, 3 = detailed)
                     ):
    """
    Perform a simple projected gradient descent with optional preconditioning 
    and Armijo backtracking line search.

    This routine minimizes a scalar objective function `fun` with known gradient `dfun`
    under optional box constraints (bounds on `x`). The descent direction is 
    optionally preconditioned by scaling the gradient to mitigate ill-conditioning (recommended).

    Parameters
    ----------
    fun : callable
        The objective function to minimize. Must take a single argument `x` (array-like).
    
    dfun : callable
        The gradient of the objective function. Must take a single argument `x` and 
        return a NumPy array of the same shape.

    x : np.ndarray
        Initial guess for the variables.

    xlb : float or list or np.ndarray, optional
        Lower bounds for each variable. Can be a scalar or array of same shape as `x`.

    xub : float or list or np.ndarray, optional
        Upper bounds for each variable. Can be a scalar or array of same shape as `x`.

    step : float, optional
        Initial step size used in the descent.

    precond : bool, optional
        Whether to apply simple gradient preconditioning (scale by max absolute value).

    iter_max : int, optional
        Maximum number of iterations.

    rstep_min : float, optional
        Minimum allowed ratio between current and initial step size. Used to detect 
        excessively small steps.

    tol : float, optional
        Convergence tolerance based on the square root of the dot product between 
        gradient and descent direction (`sqrt(|grad · d|)`).

    coeff_armijo : float, optional
        Armijo condition coefficient (must be in (0, 1)). Controls how much decrease 
        in function value is required to accept a step.

    step_increase : float, optional
        Factor by which the step size is increased after a successful iteration.

    step_decrease : float, optional
        Factor by which the step size is decreased after an unsuccessful iteration.

    verbosity : int, optional
        Level of console output:
        - 0: Silent
        - 1: Summary (recommended)
        - 2: Includes rejected steps
        - 3: Full debug info

    Returns
    -------
    result : dict
        A dictionary containing:
        
        - `x` : list of np.ndarray  
          Sequence of iterates.
        
        - `fun` : list of float  
          Objective function values at each accepted point.

        - `step` : list of float  
          Step size at each iteration.

        - `stop_criterion` : list of float  
          Value of `sqrt(|grad · descent|)` used to monitor convergence.

        - `n_fun` : int  
          Number of function evaluations.

        - `status` : int  
          Convergence status code:
          
          | Code | Meaning                                                               |
          |------|-----------------------------------------------------------------------|
          | 0    | ✅ Converged: `sqrt(|grad · descent|)` < `tol`                        |
          | 1    | ❌ Step size too small: `step / step0` < `rstep_min`                  |
          | 2    | ❌ Max iterations reached                                             |

    Notes
    -----
    - All updates are projected onto the bounds `[xlb, xub]`.
    - Preconditioning scales the gradient direction to mitigate ill-conditioning.
    - Armijo backtracking ensures sufficient decrease in the function.
    - This implementation supports box-constrained optimization only.

    Example
    -------
    >>> def f(x): return np.sum(x**2)
    >>> def df(x): return 2*x
    >>> x0 = np.array([1.0, -1.5])
    >>> result = gradient_descent(f, df, x0, xlb=-2.0, xub=2.0)
    >>> result["x"][-1]
    array([0., 0.])
    """
    if verbosity >=1 :
        print(f"---- Start gradient descent")
    x = np.minimum(xub, np.maximum(xlb, x))
    grad = dfun(x)
    result = {}
    result["x"] = [x]
    result["step"] = [step]
    result["fun"] = [fun(x)]
    result["stop_criterion"] = [np.sqrt(np.dot(grad, step * grad))]
    result["n_fun"] = 0 
    if verbosity >=1 : 
        print(f"it {result["n_fun"]} | fun = {result["fun"][-1]:.5e} | step = {step:.2e}")
    while result["n_fun"] < iter_max :
        result["n_fun"] += 1
        result["step"].append(step)

        # update and project
        if precond : # preconditioning
            precond_grad = grad / np.max(np.abs(grad))
            xTest = np.minimum(xub, np.maximum(xlb, x - step * precond_grad))
            # a second time to disqualify the grad associated to bound-constrained nodes
            ddir = xTest - x 
            precond_grad = precond_grad / np.max(np.sqrt(-precond_grad*ddir/step)) 
            xTest = np.minimum(xub, np.maximum(xlb, x - step * precond_grad))
        else : xTest = xTest = np.minimum(xub, np.maximum(xlb, x - step * grad))
        if verbosity >=3 :  print(f"{x = }")

        ddir = xTest - x  # descent direction
        funTest = fun(xTest)

        # check decrement:
        if funTest < result["fun"][-1] + coeff_armijo * np.dot(grad, ddir): # accept
            result["stop_criterion"].append(np.sqrt(-np.dot(grad, ddir)))
            x = xTest
            grad = dfun(x)
            result["x"] .append(x)
            result["fun"].append(funTest)
            step *= step_increase
            if verbosity >=1 : 
                print(f"it {result["n_fun"]} | fun = {funTest:.5e} | step = {step:.2e} | crit = {result["stop_criterion"][-1]:.2e} (accepted ✅)")

        else :  # reject
            step *= step_decrease
            if verbosity >=2 : 
                print(f"it {result["n_fun"]} | fun = {funTest:.5e} | step = {step:.2e} | crit = {result["stop_criterion"][-1]:.2e} (rejected ❌)")

        # stop criteria
        if result["stop_criterion"][-1] < tol:
            result["status"] = 0
            if verbosity >=1 :  print(f"---- Stopped after {result["n_fun"]} iterations because sqrt(|grad * descent|) < tol.")
            return result
        
        if step / result["step"][0] < rstep_min:
            result["status"] = 1
            if verbosity >=1 :  print(f"---- Stopped after {result["n_fun"]} iterations because step/step0 < rstep_min size to small.")
            return result
        
    result["status"] = 2
    if verbosity >=1 :  print(f"---- Stopped after {result["n_fun"]} because maximum number of iterations reached.")
    return result

##################################################################################################################################
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
