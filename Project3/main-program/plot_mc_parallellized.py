import numpy as np
import matplotlib.pyplot as plt

def read_file(filename):
    infile = open(filename, "r")

    N = []
    time = []

    for line in infile:
        n, t = line.split()

        N.append(float(n))
        time.append(float(t))

    return np.array(N), np.array(time)

N0, time0 = read_file("parallellized_improved_mc_O0")
N1, time1 = read_file("parallellized_improved_mc_O1")
N2, time2 = read_file("parallellized_improved_mc_O2")
N3, time3 = read_file("parallellized_improved_mc_O3")
N4, time4 = read_file("parallellized_improved_mc_Ofast")

plt.loglog(N0, time0, label="-O0")
plt.loglog(N1, time1, label="-O1")
plt.loglog(N2, time2, label="-O2")
plt.loglog(N3, time3, label="-O3")
plt.loglog(N4, time4, label="-Ofast")

plt.xlabel("N", size=14)
plt.ylabel("Time [s]", size=14)
plt.title("Time for Improved Monte Carlo with Different Compiler Flags", size=18)
plt.legend()
plt.show()
