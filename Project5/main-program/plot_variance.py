import numpy as np
import matplotlib.pyplot as plt

def read_file(filename):
    infile = open(filename, "r")
    T = []
    V = []
    for line in infile:
        t, var = line.split()
        T.append(float(t))
        V.append(float(var))

    plt.plot(T, V)

def nl(l):
    return 1 + (3*l)/(1 - l)


# Variance for simple model and saving model
plt.subplot(2, 2, 1)
read_file("5a_N500_m500_l0_a0_g0_3e4.txt")
plt.title("Variance for the Simple Model", size=18)
plt.xlabel("Number of Transactions", size=14)
plt.ylabel("Variance", size=14)
plt.legend(["$\lambda = 0.0$"])

plt.subplot(2, 2, 2)
read_file("5c_N500_m500_l025_a0_g0_3e4.txt")
plt.title("Variance for the Saving Model", size=18)
plt.xlabel("Number of Transactions", size=14)
plt.ylabel("Variance", size=14)
plt.legend(["$\lambda = 0.25$"])

plt.subplot(2, 2, 3)
read_file("5c_N500_m500_l05_a0_g0_3e4.txt")
plt.title("Variance for the Saving Model", size=18)
plt.xlabel("Number of Transactions", size=14)
plt.ylabel("Variance", size=14)
plt.legend(["$\lambda = 0.5$"])

plt.subplot(2, 2, 4)
read_file("5c_N500_m500_l09_a0_g0_3e4.txt")
plt.title("Variance for the Saving Model", size=18)
plt.xlabel("Number of Transactions", size=14)
plt.ylabel("Variance", size=14)
plt.legend(["$\lambda = 0.9$"])
plt.show()


# lambda = 0, varying alpha
read_file("5d_N500_m500_l0_a05_g0_middle4.txt")
read_file("5d_N500_m500_l0_a1_g0_middle4.txt")
read_file("5d_N500_m500_l0_a15_g0_middle4.txt")
read_file("5d_N500_m500_l0_a2_g0_middle4.txt")

plt.title("Variance for varying $\\alpha$, $\lambda = 0$, $N = 500$", size=18)
plt.xlabel("Time Steps", size=14)
plt.ylabel("Variance", size=14)
plt.legend(["$\\alpha = 0.5$", "$\\alpha = 1.0$", "$\\alpha = 1.5$", "$\\alpha = 2.0$"])
plt.show()


# lambda = 0.9, varying alpha, N = 500
read_file("5d_N500_m500_l09_a05_g0_middle4.txt")
read_file("5d_N500_m500_l09_a1_g0_middle4.txt")
read_file("5d_N500_m500_l09_a15_g0_middle4.txt")
read_file("5d_N500_m500_l09_a2_g0_middle4.txt")

plt.title("Variance with varying $\\alpha$, $\lambda = 0.9, N = 500$", size=18)
plt.xlabel("Time Steps", size=14)
plt.ylabel("Variance", size=14)
plt.legend(["$\\alpha = 0.5$", "$\\alpha = 1.0$", "$\\alpha = 1.5$", "$\\alpha = 2.0$"])
plt.show()


# Variance for saving model, comparing N = 500, N = 1000
read_file("5d_N500_m500_l09_a1_g0_middle4.txt")
read_file("5d_N1000_m500_l09_a1_g0_2.txt")

plt.title("Variance with $\lambda = 0.9, \\alpha = 1, \gamma = 0$", size=18)
plt.xlabel("Time Steps", size=14)
plt.ylabel("Variance", size=14)
plt.legend(["$N = 500$", "$N = 1000$"])
plt.show()


# lambda = 0, alpha = 1, N = 1000, varying gamma
read_file("5e_N1000_m500_l0_a1_g0_2.txt")
read_file("5e_N1000_m500_l0_a1_g1_2.txt")
read_file("5e_N1000_m500_l0_a1_g2_2.txt")
read_file("5e_N1000_m500_l0_a1_g3_2.txt")
read_file("5e_N1000_m500_l0_a1_g4_2.txt")

plt.title("Variance with varying $\gamma$ and $\lambda = 0, \\alpha = 1$", size=18)
plt.xlabel("Time Steps", size=14)
plt.ylabel("Variance", size=14)
plt.legend(["$\gamma = 0.0$", "$\gamma = 1.0$", "$\gamma = 2.0$", "$\gamma = 3.0$", "$\gamma = 4.0$", "$\\alpha = 2$"])
plt.show()


# lambda = 0.9, alpha = 2, N = 1000, varying gamma
read_file("5e_N1000_m500_l09_a2_g0_2.txt")
read_file("5e_N1000_m500_l09_a2_g1_2.txt")
read_file("5e_N1000_m500_l09_a2_g2_2.txt")
read_file("5e_N1000_m500_l09_a2_g3_2.txt")
read_file("5e_N1000_m500_l09_a2_g4_2.txt")

plt.title("Variance with varying $\gamma$ and $\lambda = 0.9, \\alpha = 2$", size=18)
plt.xlabel("Time Steps", size=14)
plt.ylabel("Variance", size=14)
plt.legend(["$\gamma = 0.0$", "$\gamma = 1.0$", "$\gamma = 2.0$", "$\gamma = 3.0$", "$\gamma = 4.0$"])
plt.show()
