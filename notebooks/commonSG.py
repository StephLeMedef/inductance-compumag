import ngsolve as ngs
import numpy as np
import matplotlib.pyplot as plt
from IPython import display

R = ngs.CoefficientFunction((0, 1, -1, 0), dims=(2, 2))  # Rotation matrix of -90°
Id = ngs.CoefficientFunction((1, 0, 0, 1), dims=(2, 2))  # Identity matrix

# Magnetostatic constants
mu0 = 4e-7 * np.pi
nu0 = 1 / mu0
mur = 1000
nuIron = 1 / (mur * mu0)


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
