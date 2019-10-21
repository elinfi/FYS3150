import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return np.exp(-4*x)


N = 100
r = np.linspace(0, 5, N)
GL = f(r)

plt.plot(abs(r), GL)
plt.xlabel("$|\mathbf{r}|$")
plt.ylabel("$e^{-4r}$")
plt.grid()
plt.show()
