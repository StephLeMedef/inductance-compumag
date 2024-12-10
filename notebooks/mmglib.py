import ctypes
import os

# Load the shared library
if os.name == "nt":  # Windows
    try:
        mmg2d = ctypes.CDLL("lib\\mmg\\mmg2d.dll")
    except FileNotFoundError:
        mmg2d = ctypes.CDLL("..\\lib\\mmg\\mmg2d.dll")
else:
    raise OSError("This script is only configured for Windows environments.")

# Variadic argument constants
MMG5_ARG_start = 1
MMG5_ARG_ppMesh = 2
MMG5_ARG_ppMet = 4
MMG5_ARG_end = 10

# Return values
MMG5_SUCCESS = 0  # Return value for success
MMG5_LOWFAILURE = 1  # Return value if the remesh process failed but we can save a conform mesh
MMG5_STRONGFAILURE = 2  # Return value if the remesh process failed and the mesh is non-conform

# Implicit boundary in iso mode
MG_ISO = 10

# Default references in iso mode
MG_PLUS = 2  # Positive domain
MG_MINUS = 3  # Negative domain

# Variadic arguments
MMG5_ARG_start = 1  # To begin a list of variadic arguments
MMG5_ARG_ppMesh = 2  # Pointer toward a MMG5_pMesh structure
MMG5_ARG_ppLs = 3  # Pointer toward a MMG5_pSol structure storing a level-set
MMG5_ARG_ppMet = 4  # Pointer toward a MMG5_pSol structure storing a metric
MMG5_ARG_ppDisp = 5  # Pointer toward a MMG5_pSol structure storing a displacement
MMG5_ARG_ppSols = 6  # Pointer toward an array of MMG5_Sol structures
MMG5_ARG_pMesh = 7  # MMG5_pMesh structure
MMG5_ARG_pMet = 8  # MMG5_pSol structure storing a metric field
MMG5_ARG_pDisp = 9  # MMG5_pSol structure storing a displacement field
MMG5_ARG_end = 10  # To end a list of variadic arguments

# Other constants
MMG5_NSOLS_MAX = 100  # Maximal number of solutions per entity
MMG5_FILENAME_LEN_MAX = 255  # Maximal length of filenames

# Multimat mode entities
MMG5_MMAT_NoSplit = 0  # Entity that must not be split in multimat mode
MMG5_MMAT_Split = 1  # Entity that must be split in multimat mode


MMG5_Vertex = 1  # MMG5 entity type for vertex
MMG5_Scalar = 1  # Scalar solution type

MMG2D_IPARAM_verbose = 0  # [-1..10], Level of verbosity
MMG2D_IPARAM_mem = 1  # [n/-1], Max memory size in Mbytes or keep the default value
MMG2D_IPARAM_debug = 2  # [1/0], Turn on/off debug mode
MMG2D_IPARAM_angle = 3  # [1/0], Turn on/off angle detection
MMG2D_IPARAM_iso = 4  # [1/0], Enable level-set discretization
MMG2D_IPARAM_isosurf = 5  # [1/0], Enable level-set discretization on the surface part only
MMG2D_IPARAM_opnbdy = 6  # [1/0], Preserve edges at the interface of 2 domains with same reference
MMG2D_IPARAM_lag = 7  # [-1/0/1/2], Enable Lagrangian motion
MMG2D_IPARAM_3dMedit = 8  # [0/1/2], Read/write 2D mesh in 3D (Medit only)
MMG2D_IPARAM_optim = 9  # [1/0], Optimize mesh keeping its initial edge sizes
MMG2D_IPARAM_noinsert = 10  # [1/0], Avoid/allow vertex insertion
MMG2D_IPARAM_noswap = 11  # [1/0], Avoid/allow edge or face flipping
MMG2D_IPARAM_nomove = 12  # [1/0], Avoid/allow vertex relocation
MMG2D_IPARAM_nosurf = 13  # [1/0], Avoid/allow surface modifications
MMG2D_IPARAM_nreg = 14  # [0/1], Enable normal regularization
MMG2D_IPARAM_xreg = 15  # [0/1], Enable regularization by moving vertices
MMG2D_IPARAM_numsubdomain = 16  # [0/n], Save only the subdomain n (0 == all subdomains)
MMG2D_IPARAM_numberOfLocalParam = 17  # [n], Number of local parameters
MMG2D_IPARAM_numberOfLSBaseReferences = 18  # [n], Number of base references for bubble removal
MMG2D_IPARAM_numberOfMat = 19  # [n], Number of materials in level-set mode
MMG2D_IPARAM_anisosize = 20  # [1/0], Turn on/off anisotropic metric creation when no metric is provided
MMG2D_IPARAM_nosizreq = 21  # [0/1], Allow/avoid overwriting of sizes at required vertices (advanced usage)
MMG2D_DPARAM_angleDetection = 22  # [val], Threshold for angle detection
MMG2D_DPARAM_hmin = 23  # [val], Minimal edge length
MMG2D_DPARAM_hmax = 24  # [val], Maximal edge length
MMG2D_DPARAM_hsiz = 25  # [val], Constant edge length
MMG2D_DPARAM_hausd = 26  # [val], Global Hausdorff distance (on all the boundary surfaces of the mesh)
MMG2D_DPARAM_hgrad = 27  # [val], Global gradation
MMG2D_DPARAM_hgradreq = 28  # [val], Global gradation on required entities (advanced usage)
MMG2D_DPARAM_ls = 29  # [val], Function value where the level set is to be discretized
MMG2D_DPARAM_xreg = 30  # [val], Relaxation parameter for coordinate regularization (0 < val < 1)
MMG2D_DPARAM_rmc = 31  # [-1/val], Remove small disconnected components in level-set mode
MMG2D_IPARAM_nofem = 32  # [1/0], Do not attempt to make the mesh suitable for finite-element computations
MMG2D_IPARAM_isoref = 33  # [0/n], Iso-surface boundary material reference

# Define the argument and return types for MMG2D_saveSol
mmg2d.MMG2D_saveSol.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_char_p]
mmg2d.MMG2D_saveSol.restype = ctypes.c_int

# Define function prototypes
mmg2d.MMG2D_Init_mesh.argtypes = [
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_void_p),
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_void_p),
    ctypes.c_int,
]
mmg2d.MMG2D_Init_mesh.restype = ctypes.c_int

mmg2d.MMG2D_loadMesh.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
mmg2d.MMG2D_loadMesh.restype = ctypes.c_int

mmg2d.MMG2D_loadSol.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_char_p]
mmg2d.MMG2D_loadSol.restype = ctypes.c_int

mmg2d.MMG2D_mmg2dlib.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
mmg2d.MMG2D_mmg2dlib.restype = ctypes.c_int

mmg2d.MMG2D_saveMesh.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
mmg2d.MMG2D_saveMesh.restype = ctypes.c_int

mmg2d.MMG2D_saveSol.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_char_p]
mmg2d.MMG2D_saveSol.restype = ctypes.c_int

mmg2d.MMG2D_Free_all.argtypes = [
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_void_p),
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_void_p),
    ctypes.c_int,
]
mmg2d.MMG2D_Free_all.restype = None

mmg2d.MMG2D_Get_meshSize.argtypes = [
    ctypes.c_void_p,
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
]
mmg2d.MMG2D_Get_meshSize.restype = ctypes.c_int

mmg2d.MMG2D_Set_solSize.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_int]
mmg2d.MMG2D_Set_solSize.restype = ctypes.c_int

mmg2d.MMG2D_Set_scalarSol.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_int]
mmg2d.MMG2D_Set_scalarSol.restype = ctypes.c_int

mmg2d.MMG2D_Set_dparameter.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
mmg2d.MMG2D_Set_dparameter.restype = ctypes.c_int

mmg2d.MMG2D_Set_iparameter.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
mmg2d.MMG2D_Set_iparameter.restype = ctypes.c_int

mmg2d.MMG2D_Set_vertex.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_int]
mmg2d.MMG2D_Set_vertex.restype = ctypes.c_int

mmg2d.MMG2D_Set_edge.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int]
mmg2d.MMG2D_Set_edge.restype = ctypes.c_int

mmg2d.MMG2D_Set_triangle.argtypes = [
    ctypes.c_void_p,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
]
mmg2d.MMG2D_Set_triangle.restype = ctypes.c_int

mmg2d.MMG2D_Set_meshSize.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int]
mmg2d.MMG2D_Set_meshSize.restype = ctypes.c_int

mmg2d.MMG2D_Get_meshSize.argtypes = [
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_void_p,
]
mmg2d.MMG2D_Get_meshSize.restype = ctypes.c_int

mmg2d.MMG2D_Get_vertex.argtypes = [
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_void_p,
]
mmg2d.MMG2D_Get_vertex.restype = ctypes.c_int

mmg2d.MMG2D_Get_edge.argtypes = [
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_void_p,
]
mmg2d.MMG2D_Get_edge.restype = ctypes.c_int

mmg2d.MMG2D_Get_triangle.argtypes = [
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_void_p,
]
mmg2d.MMG2D_Get_triangle.restype = ctypes.c_int

mmg2d.MMG2D_mmg2dmov.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
mmg2d.MMG2D_mmg2dmov.restype = ctypes.c_int


import netgen


def init_mmgmesh_from_ngmesh(ngmesh, mmgMesh):
    nb_vertices = len(ngmesh.Points())
    nb_edges = len(ngmesh.Elements1D())
    nb_triangles = len(ngmesh.Elements2D())

    code = mmg2d.MMG2D_Set_meshSize(mmgMesh, nb_vertices, nb_triangles, 0, nb_edges)
    if code != 1:
        raise RuntimeError(f"Failed to initialize meshsize. {code=}")

    # Initialize base index and max index
    base_index = 0
    max_index = 0

    # Index map for (dimension, index)
    index_map = {}

    # Helper function for consistent indexing
    def get_index(i, dim):
        nonlocal max_index
        index = base_index + i
        max_index = max(max_index, index)
        index_map[(dim, i)] = index
        return index

    # Write Vertices
    for i, point in enumerate(ngmesh.Points()):
        x, y, z = point.p
        code = mmg2d.MMG2D_Set_vertex(mmgMesh, x, y, get_index(1, 0), i + 1)
        if code != 1:
            raise RuntimeError(f"Failed to set vertex. {code=}")

    # Update base index for edges
    base_index = max_index

    # Write Edges
    indices = set()
    for i, edge in enumerate(ngmesh.Elements1D()):
        p0, p1 = edge.points
        indices.add(edge.index)
        label = get_index(edge.index, edge.index)
        code = mmg2d.MMG2D_Set_edge(mmgMesh, p0.nr, p1.nr, label, i + 1)
        if code != 1:
            raise RuntimeError(f"Failed to set edge. {code=}")

    bc_names = dict()
    for indice in indices:
        bc_names[indice] = ngmesh.GetBCName(indice - 1)

    # Update base index for triangles
    base_index = max_index

    # Write Triangles
    for i, triangle in enumerate(ngmesh.Elements2D()):
        p0, p1, p2 = triangle.points
        label = get_index(triangle.index, 2)
        code = mmg2d.MMG2D_Set_triangle(mmgMesh, p0.nr, p1.nr, p2.nr, label, i + 1)
        if code != 1:
            raise RuntimeError(f"Failed to set triangle. {code=}")

    # Update base index for tetrahedra
    base_index = max_index

    region_names = ngmesh.GetRegionNames(2)

    return mmgMesh, region_names, bc_names


def copy_ngmesh(ngmesh):
    new_ngmesh = netgen.meshing.Mesh(dim=2)

    # Write Vertices
    pnums = []
    for i, point in enumerate(ngmesh.Points()):
        x, y, z = point.p
        pnums.append(new_ngmesh.Add(netgen.meshing.MeshPoint(netgen.meshing.Pnt(x, y, 0))))

    # Write Edges
    indices = set()
    for i, edge in enumerate(ngmesh.Elements1D()):
        p0, p1 = edge.points
        indices.add(edge.index)
        new_ngmesh.Add(netgen.meshing.Element1D([pnums[p0.nr - 1], pnums[p1.nr - 1]], index=edge.index))

    for indice in indices:
        new_ngmesh.SetBCName(indice - 1, ngmesh.GetBCName(indice - 1))

    region_names = ngmesh.GetRegionNames(2)
    for i in range(len(region_names)):
        new_ngmesh.AddRegion(region_names[i], dim=2)

    # Write Triangles
    for i, triangle in enumerate(ngmesh.Elements2D()):
        p0, p1, p2 = triangle.points
        vertices = [
            pnums[p0.nr - 1],
            pnums[p1.nr - 1],
            pnums[p2.nr - 1],
        ]
        new_ngmesh.Add(netgen.meshing.Element2D(triangle.index, vertices))

    return new_ngmesh


def init_ngmesh_from_mmgmesh(mmgMesh, region_names, bc_names):
    ngmesh = netgen.meshing.Mesh(dim=2)
    nb_vertices = ctypes.c_int(0)
    nb_edges = ctypes.c_int(0)
    nb_triangles = ctypes.c_int(0)
    nb_quads = ctypes.c_int(0)

    ctypes.byref(mmgMesh)
    code = mmg2d.MMG2D_Get_meshSize(
        mmgMesh, ctypes.byref(nb_vertices), ctypes.byref(nb_triangles), ctypes.byref(nb_quads), ctypes.byref(nb_edges)
    )
    if code != 1:
        raise RuntimeError(f"Failed to get meshsize. {code=}")

    # Add points
    pnums = []
    corner_points = 0
    required_points = 0
    for i in range(nb_vertices.value):
        x = ctypes.c_double(0)
        y = ctypes.c_double(0)
        index = ctypes.c_int(0)
        is_corder = ctypes.c_int(0)
        is_required = ctypes.c_int(0)
        code = mmg2d.MMG2D_Get_vertex(
            mmgMesh,
            ctypes.byref(x),
            ctypes.byref(y),
            ctypes.byref(index),
            ctypes.byref(is_corder),
            ctypes.byref(is_required),
        )
        if code != 1:
            raise RuntimeError(f"Failed to get vertex. {code=}")
        pnums.append(ngmesh.Add(netgen.meshing.MeshPoint(netgen.meshing.Pnt(x.value, y.value, 0))))

        if is_corder.value == 1:
            corner_points += 1
        if is_required.value == 1:
            required_points += 1
    # print(f"found {corner_points=} corner points")
    # print(f"found {required_points=} required points")

    for indice, name in bc_names.items():
        ngmesh.SetBCName(indice, name)

    # Add edges
    ridge_edges = 0
    required_edges = 0
    indices = set()
    for i in range(nb_edges.value):
        p0nr = ctypes.c_int(0)
        p1nr = ctypes.c_int(0)
        index = ctypes.c_int(0)
        is_ridge = ctypes.c_int(0)
        is_required = ctypes.c_int(0)
        code = mmg2d.MMG2D_Get_edge(
            mmgMesh,
            ctypes.byref(p0nr),
            ctypes.byref(p1nr),
            ctypes.byref(index),
            ctypes.byref(is_ridge),
            ctypes.byref(is_required),
        )
        if code != 1:
            raise RuntimeError(f"Failed to get edge. {code=}")
        indices.add(index.value)
        ngmesh.Add(netgen.meshing.Element1D([pnums[p0nr.value - 1], pnums[p1nr.value - 1]], index=index.value))
        if is_ridge.value == 1:
            ridge_edges += 1
        if is_required.value == 1:
            required_edges += 1

    # print(f"found {ridge_edges=} ridge edges")
    # print(f"found {required_edges=} required edges")

    # Recovering region names and indices
    set_index = set()
    required_triangles = 0
    for i in range(nb_triangles.value):
        p0nr = ctypes.c_int(0)
        p1nr = ctypes.c_int(0)
        p2nr = ctypes.c_int(0)
        index = ctypes.c_int(0)
        is_required = ctypes.c_int(0)
        code = mmg2d.MMG2D_Get_triangle(
            mmgMesh,
            ctypes.byref(p0nr),
            ctypes.byref(p1nr),
            ctypes.byref(p2nr),
            ctypes.byref(index),
            ctypes.byref(is_required),
        )
        if code != 1:
            raise RuntimeError(f"Failed to get triangle. {code=}")
        set_index.add(index.value)
        if is_required.value == 1:
            required_triangles += 1
    # print(f"found {required_triangles=} required triangles")

    netgen_region_ids = []
    for i in range(len(region_names)):
        netgen_region_ids.append(ngmesh.AddRegion(region_names[i], dim=2))

    mmg_region_indices = list(sorted(set_index))
    dict_region_names = dict(map(lambda i, j: (i, j), mmg_region_indices, netgen_region_ids))

    # Add triangles
    for i in range(nb_triangles.value):
        p0nr = ctypes.c_int(0)
        p1nr = ctypes.c_int(0)
        p2nr = ctypes.c_int(0)
        index = ctypes.c_int(0)
        is_required = ctypes.c_int(0)
        code = mmg2d.MMG2D_Get_triangle(
            mmgMesh,
            ctypes.byref(p0nr),
            ctypes.byref(p1nr),
            ctypes.byref(p2nr),
            ctypes.byref(index),
            ctypes.byref(is_required),
        )
        if code != 1:
            raise RuntimeError(f"Failed to get triangle. {code=}")

        vertices = [
            pnums[p0nr.value - 1],
            pnums[p1nr.value - 1],
            pnums[p2nr.value - 1],
        ]
        ngmesh.Add(netgen.meshing.Element2D(dict_region_names[index.value], vertices))

    return ngmesh


def run_adapt(ngmesh, **kwargs):
    """
    Adapt the mesh using MMG2D with flexible keyword arguments.

    Parameters:
        ngmesh: Input mesh (Netgen mesh object)
        kwargs: Key-value pairs for MMG2D parameters (e.g., hmin=0.003, hgrad=1.3)

    Returns:
        new_ngmesh: Adapted Netgen mesh
        return_code: Status code (1 for success, 0 for failure)
    """
    # Define MMG parameter mappings
    mmg_iparams = {
        "verbose": MMG2D_IPARAM_verbose,
        "mem": MMG2D_IPARAM_mem,
        "debug": MMG2D_IPARAM_debug,
        "angle": MMG2D_IPARAM_angle,
        "iso": MMG2D_IPARAM_iso,
        "isosurf": MMG2D_IPARAM_isosurf,
        "opnbdy": MMG2D_IPARAM_opnbdy,
        "lag": MMG2D_IPARAM_lag,
        "3dMedit": MMG2D_IPARAM_3dMedit,
        "optim": MMG2D_IPARAM_optim,
        "noinsert": MMG2D_IPARAM_noinsert,
        "noswap": MMG2D_IPARAM_noswap,
        "nomove": MMG2D_IPARAM_nomove,
        "nosurf": MMG2D_IPARAM_nosurf,
        "nreg": MMG2D_IPARAM_nreg,
        "xreg": MMG2D_IPARAM_xreg,
        "numsubdomain": MMG2D_IPARAM_numsubdomain,
        "numberOfLocalParam": MMG2D_IPARAM_numberOfLocalParam,
        "numberOfLSBaseReferences": MMG2D_IPARAM_numberOfLSBaseReferences,
        "numberOfMat": MMG2D_IPARAM_numberOfMat,
        "anisosize": MMG2D_IPARAM_anisosize,
        "nosizreq": MMG2D_IPARAM_nosizreq,
        "nofem": MMG2D_IPARAM_nofem,
        "isoref": MMG2D_IPARAM_isoref,
    }

    mmg_dparams = {
        "angleDetection": MMG2D_DPARAM_angleDetection,
        "hmin": MMG2D_DPARAM_hmin,
        "hmax": MMG2D_DPARAM_hmax,
        "hsiz": MMG2D_DPARAM_hsiz,
        "hausd": MMG2D_DPARAM_hausd,
        "hgrad": MMG2D_DPARAM_hgrad,
        "hgradreq": MMG2D_DPARAM_hgradreq,
        "ls": MMG2D_DPARAM_ls,
        "xreg": MMG2D_DPARAM_xreg,
        "rmc": MMG2D_DPARAM_rmc,
    }

    # Initialize MMG2D structures
    mmgMesh = ctypes.c_void_p()
    mmgSol = ctypes.c_void_p()

    code = mmg2d.MMG2D_Init_mesh(
        MMG5_ARG_start,
        MMG5_ARG_ppMesh,
        ctypes.byref(mmgMesh),
        MMG5_ARG_ppMet,
        ctypes.byref(mmgSol),
        MMG5_ARG_end,
    )
    if code != 1:
        raise RuntimeError(f"Failed to initialize MMG2D structures. {code=}")

    # Initialize a MMG mesh structure from the Netgen mesh structure
    mmgMesh, region_names, bc_names = init_mmgmesh_from_ngmesh(ngmesh, mmgMesh)

    # Set integer parameters (IPARAM)
    for param_name, param_value in kwargs.items():
        if param_name in mmg_iparams:
            param_id = mmg_iparams[param_name]
            if mmg2d.MMG2D_Set_iparameter(mmgMesh, mmgSol, param_id, param_value) != 1:
                raise RuntimeError(f"Failed to set IPARAM {param_name} (ID {param_id})")

    # Set double parameters (DPARAM)
    for param_name, param_value in kwargs.items():
        if param_name in mmg_dparams:
            param_id = mmg_dparams[param_name]
            if mmg2d.MMG2D_Set_dparameter(mmgMesh, mmgSol, param_id, param_value) != 1:
                raise RuntimeError(f"Failed to set DPARAM {param_name} (ID {param_id})")

    # Perform the mesh adaptation
    return_code = 1
    code = mmg2d.MMG2D_mmg2dlib(mmgMesh, mmgSol)
    if code == MMG5_STRONGFAILURE:
        raise RuntimeError("BAD ENDING OF MMG2DLIB: UNABLE TO SAVE MESH")
    elif code == MMG5_LOWFAILURE:
        print("BAD ENDING OF MMG2DLIB")
        return_code = 0
        new_ngmesh = ngmesh
    else:
        # Initialize a MMG Netgen structure from the MMG mesh structure
        return_code = 1
        new_ngmesh = init_ngmesh_from_mmgmesh(mmgMesh, region_names, bc_names)

    # Free MMG2D structures
    mmg2d.MMG2D_Free_all(
        MMG5_ARG_start,
        MMG5_ARG_ppMesh,
        ctypes.byref(mmgMesh),
        MMG5_ARG_ppMet,
        ctypes.byref(mmgSol),
        MMG5_ARG_end,
    )

    return new_ngmesh, return_code


def run_lagrangian_motion(ngmesh):
    # Initialize MMG2D structures
    mmgMesh = ctypes.c_void_p()
    empty_metric = ctypes.c_void_p()
    my_displacement = ctypes.c_void_p()

    # Initialize mmg mesh structure
    code = mmg2d.MMG2D_Init_mesh(
        MMG5_ARG_start,  # Start of arguments
        MMG5_ARG_ppMesh,
        ctypes.byref(mmgMesh),  # Mesh structure
        MMG5_ARG_ppMet,
        ctypes.byref(empty_metric),  # Metric structure
        MMG5_ARG_ppDisp,
        ctypes.byref(my_displacement),  # Metric structure
        MMG5_ARG_end,  # End of arguments
    )
    if code != 1:
        raise RuntimeError(f"Failed to initialize MMG2D structures. {code=}")

    # Initialize a mmg mesh structure from the netgen mesh structure, save region names information
    mmgMesh, region_names = init_mmgmesh_from_ngmesh(ngmesh, mmgMesh)

    # Set some parameters
    code = mmg2d.MMG2D_Set_iparameter(mmgMesh, my_displacement, MMG2D_IPARAM_lag, 1)
    if code != 1:
        raise RuntimeError(f"Failed to set lag {code=}")

    # Perform the mesh adaptation
    code = mmg2d.MMG2D_mmg2dmov(mmgMesh, empty_metric, my_displacement)
    if code == MMG5_STRONGFAILURE:
        raise RuntimeError("BAD ENDING OF MMG2DMOV: UNABLE TO SAVE MESH")
    elif code == MMG5_LOWFAILURE:
        print("BAD ENDING OF MMG2DLIB")

    # Initialize a mmg netgen structure from the mmg mesh structure, recover region names information
    new_ngmesh = init_ngmesh_from_mmgmesh(mmgMesh, empty_metric, region_names)

    # Free MMG2D structures
    mmg2d.MMG2D_Free_all(
        MMG5_ARG_start,  # Start of arguments
        MMG5_ARG_ppMesh,
        ctypes.byref(mmgMesh),  # Mesh structure
        MMG5_ARG_ppMet,
        ctypes.byref(empty_metric),  # Metric structure
        MMG5_ARG_ppDisp,
        ctypes.byref(my_displacement),  # Metric structure
        MMG5_ARG_end,  # End of arguments
    )

    return new_ngmesh


# Define the run_mmg2d function
def run_mmg2d(file_in, file_out):
    # Initialize MMG2D structures
    mmgMesh = ctypes.c_void_p()
    mmgSol = ctypes.c_void_p()

    # Initialize mmg mesh structure
    code = mmg2d.MMG2D_Init_mesh(
        MMG5_ARG_start,  # Start of arguments
        MMG5_ARG_ppMesh,
        ctypes.byref(mmgMesh),  # Mesh structure
        MMG5_ARG_ppMet,
        ctypes.byref(mmgSol),  # Metric structure
        MMG5_ARG_end,  # End of arguments
    )
    if code != 1:
        raise RuntimeError(f"Failed to initialize MMG2D structures. {code=}")

    print(f"MMG2D initialized successfully: mmgMesh={mmgMesh}, mmgSol={mmgSol}")

    # Set some parameters
    if mmg2d.MMG2D_Set_dparameter(mmgMesh, mmgSol, MMG2D_DPARAM_hmin, 0.003) != 1:
        raise RuntimeError("Failed to set hmin")
    if mmg2d.MMG2D_Set_dparameter(mmgMesh, mmgSol, MMG2D_DPARAM_hgrad, 1.3) != 1:
        raise RuntimeError("Failed to set hgrad")
    if mmg2d.MMG2D_Set_dparameter(mmgMesh, mmgSol, MMG2D_DPARAM_hmax, 0.3) != 1:
        raise RuntimeError("Failed to set hmax")
    if mmg2d.MMG2D_Set_iparameter(mmgMesh, mmgSol, MMG2D_IPARAM_optim, 1) != 1:
        raise RuntimeError("Failed to set optim")

    # Load initial mesh file
    code = mmg2d.MMG2D_loadMesh(mmgMesh, file_in.encode("utf-8"))
    if code != 1:
        raise RuntimeError(f"Failed to load mesh file. {code=}")

    # # Save the initial mesh and solution
    # if mmg2d.MMG2D_saveMesh(mmgMesh, file_out.encode('utf-8')) != 1:
    #     raise RuntimeError("Failed to save mesh.")

    # Save the solution
    code = mmg2d.MMG2D_saveSol(mmgMesh, mmgSol, file_out.encode("utf-8"))
    if code != 1:
        raise RuntimeError(f"Failed to save solution to file: {file_out}")
    else:
        print(f"Solution saved successfully to {file_out}")

    # Perform the mesh adaptation
    code = mmg2d.MMG2D_mmg2dlib(mmgMesh, mmgSol)
    if code == MMG5_STRONGFAILURE:
        raise RuntimeError("BAD ENDING OF MMG2DLIB: UNABLE TO SAVE MESH")
    elif code == MMG5_LOWFAILURE:
        print("BAD ENDING OF MMG2DLIB")
    else:
        print("MMG2D adaptation completed successfully.")

    # Save the adapted mesh
    result = mmg2d.MMG2D_saveMesh(mmgMesh, file_out.encode("utf-8"))
    if result != 1:
        raise RuntimeError(f"Failed to save mesh to file: {file_out}")
    else:
        print(f"Mesh saved successfully to {file_out}")

    # Free MMG2D structures
    mmg2d.MMG2D_Free_all(
        MMG5_ARG_start,
        MMG5_ARG_ppMesh,
        ctypes.byref(mmgMesh),
        MMG5_ARG_ppMet,
        ctypes.byref(mmgSol),
        MMG5_ARG_end,
    )
    print("MMG2D run completed successfully.")


# Example usage
if __name__ == "__main__":
    import ngsolve as ngs

    air_gap = 5e-3
    mesh = ngs.Mesh(ngs.unit_square.GenerateMesh(maxh=0.2))
    ngmesh = run_adapt(mesh.ngmesh)
    mesh = ngs.Mesh(ngmesh)
