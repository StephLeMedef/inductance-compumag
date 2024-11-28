# test
import numpy as np
import matplotlib.pyplot as plt

X = np.linspace(0, 6, 200)
Y = np.sin(X**2) ** 2

plt.plot(X, Y, "--")
plt.show()
