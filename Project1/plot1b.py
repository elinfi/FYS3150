import numpy as np
import matplotlib.pyplot as plt


def plot(filename, N, name):
    file = open(filename, "r")

    h = 1.0/(N + 1)
    x = h*np.linspace(0, N+1, N+2)

    line = file.readline().split(",")

    for i in range(0, len(line)):
        line[i] = float(line[i])


    plt.plot(x, line, label=name)


for i in (10, 100, 1000):
    plot("1b-"+str(i)+".txt", i, "N = "+str(i))

plot("analytic_solution-" + str(1000) + ".txt", i, "Analytic Solution")
plt.title("General Thomas Algorithm")
plt.xlabel("$x$")
plt.ylabel("$v(x)$")
plt.legend()
plt.show()


log_10_h = np.array([-1.04139, -2.00432, -3.00043, -4.00004, -5, -6])
log_10_eps = np.array([-1.1797, -3.08804, -5.08005, -7.07936, -9.0049, -6.77137])

plt.plot(log_10_h, log_10_eps)
plt.plot(log_10_h, log_10_eps, "bo")
plt.title("Maximum Relative Error")
plt.xlabel("$log_{10}(h)$")
plt.ylabel("max($\epsilon$)")
plt.show()

