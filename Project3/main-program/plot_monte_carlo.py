import numpy as np
import matplotlib.pyplot as plt

def read_file(filename):
    infile = open(filename, "r")

    rel_err = []
    time = []

    for line in infile:
        err, t = line.split()

        rel_err.append((float(err)))
        time.append((float(t)))

    return np.array(rel_err), np.array(time)

#rel_err1, time1 = read_file("gauss_quadrature_more")
#rel_err2, time2 = read_file("improved_gauss_quadrature_more")

rel_err3, time3 = read_file("monte_carlo_time")
rel_err4, time4 = read_file("improved_monte_carlo_time")

#rel_err5, time5 = read_file("parallellized_monte_carlo")
#rel_err6, time6 = read_file("parallellized_improved_monte_carlo")

#plt.plot(time1, rel_err1, label="Gauss Legendre")
#plt.plot(time2, rel_err2, label="Gauss Laguerre")

plt.loglog(time3, rel_err3, label="Brute Force Monte Carlo")
plt.loglog(time4, rel_err4, label="Improved Monte Carlo")

#plt.plot(time5[:-1], rel_err5[:-1], label="Brute Force Monte Carlo")
#plt.plot(time6[:-1], rel_err6[:-1], label="Improved Monte Carlo")

plt.xlabel("Time", size=14)
plt.ylabel("Relative Error", size=14)
plt.title("Relative Error for Monte Carlo", size=18)
plt.legend()
plt.show()

#plt.loglog(time3, rel_err3, label="Brute Force Monte Carlo")
#plt.loglog(time4, rel_err4, label="Improved Monte Carlo")

#plt.xlabel("Time", size=14)
#plt.ylabel("Relative Error", size=14)
#plt.title("Relative Error for Monte Carlo", size=18)
#plt.legend()
#plt.show()
