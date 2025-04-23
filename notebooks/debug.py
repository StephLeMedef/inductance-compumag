from netgen.geom2d import SplineGeometry
import matplotlib.pyplot as plt
import numpy as np


def debug_geo(geo, ax=None):
    """Plot the geometry with matplotlib."""
    # Plot the points
    points = []
    nb_points = geo.GetNPoints()
    for i in range(nb_points):
        point = geo.GetPoint(i)
        points.append((point[0], point[1]))
    points = np.array(points)

    # Create a new figure if no ax is provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))

    ax.scatter(points[:, 0], points[:, 1], color="blue")

    # Plot the lines and splines
    nb_splines = geo.GetNSplines()
    for i in range(nb_splines):
        spline = geo.GetSpline(i)
        p_start, p_end = spline.StartPoint(), spline.EndPoint()
        ax.plot([p_start[0], p_end[0]], [p_start[1], p_end[1]], color="green")

        normal_start = spline.GetNormal(0)
        tangent_start = (-normal_start[1], normal_start[0])
        normal_end = spline.GetNormal(1)
        tangent_end = (-normal_end[1], normal_end[0])
        determinant = tangent_start[0] * tangent_end[1] - tangent_end[0] * tangent_start[1]
        if np.abs(determinant) > 1e-16:
            print("spline3 not yet supported, displaying as a line")

    ax.grid(True)
    ax.axis("equal")
    ax.set_title("Geometry")

    # Show the plot if no ax is provided
    if ax is None:
        plt.show()


def debug_points(geo: SplineGeometry, font_size=10, point_names=None, ax=None):
    """Plot the geometry with matplotlib. Displays for each point of index i its name point_names[i]."""
    # Plot the points
    points = []
    nb_points = geo.GetNPoints()
    for i in range(nb_points):
        point = geo.GetPoint(i)
        points.append((point[0], point[1]))
    points = np.array(points)

    # Create a new figure if no ax is provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))

    ax.scatter(points[:, 0], points[:, 1], color="blue")

    if point_names and len(point_names) == len(points):
        for i in range(len(points)):
            ax.text(points[i, 0], points[i, 1], point_names[i], color="blue", fontsize=font_size)
    elif point_names and len(point_names) != len(points):
        print("Warning : len(point_names) != len(points). points_names is ignored")
        for i in range(len(points)):
            ax.text(points[i, 0], points[i, 1], f"p[{i}]", color="blue", fontsize=font_size)
    else:
        for i in range(len(points)):
            ax.text(points[i, 0], points[i, 1], f"p[{i}]", color="blue", fontsize=font_size)

    # Plot the lines and splines
    nb_splines = geo.GetNSplines()
    for i in range(nb_splines):
        spline = geo.GetSpline(i)
        p_start, p_end = spline.StartPoint(), spline.EndPoint()
        ax.plot([p_start[0], p_end[0]], [p_start[1], p_end[1]], color="green")

        normal_start = spline.GetNormal(0)
        tangent_start = (-normal_start[1], normal_start[0])
        normal_end = spline.GetNormal(1)
        tangent_end = (-normal_end[1], normal_end[0])
        determinant = tangent_start[0] * tangent_end[1] - tangent_end[0] * tangent_start[1]
        if np.abs(determinant) > 1e-16:
            print("spline3 not yet supported, displaying as a line")

    ax.grid(True)
    ax.axis("equal")
    ax.set_title("Points names")

    # Show the plot if no ax is provided
    if ax is None:
        plt.show()


def debug_bc_names(geo: SplineGeometry, font_size=10, ax=None):
    """Plot the geometry with matplotlib. Displays for each boundary its boundary condition name."""
    # Plot the points
    points = []
    nb_points = geo.GetNPoints()
    for i in range(nb_points):
        point = geo.GetPoint(i)
        points.append((point[0], point[1]))

    points = np.array(points)

    # Create a new figure if no ax is provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))

    ax.scatter(points[:, 0], points[:, 1], color="blue")

    # Plot the lines and splines and annotate with region labels
    nb_splines = geo.GetNSplines()
    for i in range(nb_splines):
        spline = geo.GetSpline(i)
        p_start, p_end = spline.StartPoint(), spline.EndPoint()
        ax.plot([p_start[0], p_end[0]], [p_start[1], p_end[1]], color="green")

        # Midpoint for labels
        midpoint = ((p_start[0] + p_end[0]) / 2, (p_start[1] + p_end[1]) / 2)
        bc_name = geo.GetBCName(spline.bc)
        ax.text(midpoint[0], midpoint[1], bc_name, color="black", fontsize=font_size, ha="center")

        normal_start = spline.GetNormal(0)
        tangent_start = (-normal_start[1], normal_start[0])
        normal_end = spline.GetNormal(1)
        tangent_end = (-normal_end[1], normal_end[0])
        determinant = tangent_start[0] * tangent_end[1] - tangent_end[0] * tangent_start[1]
        if np.abs(determinant) > 1e-16:
            print("spline3 not yet supported, displaying as a line")

    ax.grid(True)
    ax.axis("equal")
    ax.set_title("Boundary Condition Names")

    # Show the plot if no ax is provided
    if ax is None:
        plt.show()

def debug_region_labels(geo: SplineGeometry, offset=1e-3, font_size=10, ax=None):
    """Plot the geometry with matplotlib. Displays for each boundary it's left and right region labels. 
    Use the offset parameter if the labels are not readable."""
    # Plot the points
    points = []
    nb_points = geo.GetNPoints()
    for i in range(nb_points):
        point = geo.GetPoint(i)
        points.append((point[0], point[1]))

    points = np.array(points)


    # Create a new figure if no ax is provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))

    ax.scatter(points[:, 0], points[:, 1], color="blue")

    # Plot the lines and splines and annotate with region labels
    nb_splines = geo.GetNSplines()
    for i in range(nb_splines):
        spline = geo.GetSpline(i)
        p_start, p_end = spline.StartPoint(), spline.EndPoint()
        ax.plot([p_start[0], p_end[0]], [p_start[1], p_end[1]], color="green")

        # Midpoint for labels
        midpoint = ((p_start[0] + p_end[0]) / 2, (p_start[1] + p_end[1]) / 2)
        dx, dy = p_end[0] - p_start[0], p_end[1] - p_start[1]
        normal = np.array([-dy, dx])
        normal_length = np.linalg.norm(normal)
        if normal_length != 0:
            normal = normal / normal_length
        left_domain_position = (midpoint[0] + normal[0] * offset, midpoint[1] + normal[1] * offset)
        right_domain_position = (midpoint[0] - normal[0] * offset, midpoint[1] - normal[1] * offset)
        left_domain_label = spline.leftdom
        right_domain_label = spline.rightdom
        ax.text(
            left_domain_position[0],
            left_domain_position[1],
            left_domain_label,
            color="red",
            fontsize=font_size,
            ha="center",
        )
        ax.text(
            right_domain_position[0],
            right_domain_position[1],
            right_domain_label,
            color="red",
            fontsize=font_size,
            ha="center",
        )

        normal_start = spline.GetNormal(0)
        tangent_start = (-normal_start[1], normal_start[0])
        normal_end = spline.GetNormal(1)
        tangent_end = (-normal_end[1], normal_end[0])
        determinant = tangent_start[0] * tangent_end[1] - tangent_end[0] * tangent_start[1]
        if np.abs(determinant) > 1e-16:
            print("spline3 not yet supported, displaying as a line")

    ax.grid(True)
    ax.axis("equal")
    ax.set_title("Region labels")

    # Show the plot if no ax is provided
    if ax is None:
        plt.show()
