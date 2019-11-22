import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("5a")
print(data.shape)

max = np.max(data)
min = np.min(data)

bins = int((max - min)*0.05)
print(bins)
plt.hist(data.flatten(), bins, (min, max), density=True)

m = np.linspace(0, max, 10000)
beta = 1/500
plt.plot(m, beta*np.exp(-beta*m))
plt.show()


