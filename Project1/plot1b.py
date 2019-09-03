import numpy as np
import matplotlib.pyplot as plt


def plot(filename, N):
    file = open(filename, "r")

    h = 1.0/(N + 1)
    x = h*np.linspace(0, 1, N)
    print(x)

    line = file.readline().split(",")

    for i in range(0, len(line)):
        line[i] = float(line[i])


    plt.plot(x, line, label="N = "+str(N))
    plt.xlabel("x")
    plt.ylabel("v(x)")
    plt.legend()
    plt.show()


for i in (10, 100, 1000):
    plot("b"+str(i)+".txt", i)



