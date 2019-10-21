# This Python file uses the following encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt

def d(x):
    return np.exp(-4*x)

x = np.linspace(0,5, 1000)
plt.title("$y= e^{-4x}$")
plt.plot(x, d(x))
plt.xlabel("x")
plt.ylabel("f(x)")
plt.grid(b=True)
plt.axes([0, 4, 0.6, 0.01])

plt.legend()
plt.show()
