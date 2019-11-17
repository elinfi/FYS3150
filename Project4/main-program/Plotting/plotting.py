import matplotlib.pyplot as plt
import numpy as np

def read_file(filename):
    infile = open(filename, "r")

    infile.readline()

    N = []
    E = []
    M = []

    for line in infile:
        n, e, m, mabs = line.split()
        N.append(float(n))
        E.append(float(e))
        M.append(float(mabs))

    return N, E, M

def read_file2(filename):
    infile = open(filename, "r")

    infile.readline()

    N = []
    E = []
    M = []

    for line in infile:
        n, e, m = line.split()
        N.append(float(n))
        E.append(float(e))
        M.append(float(m))

    return N, E, M

No, Eo, Mo = read_file2("4c_T1_10_3e5_100_false.txt")
Nr, Er, Mr = read_file2("4c_T1_10_3e5_100_true.txt")
#Ns, Es, Ms = read_file("4c_T1_1_3e5_1_true_L20.txt")

plt.plot(No, Eo, label="Ordered spin")
plt.plot(Nr, Er, label="Random spin")
#plt.plot(Ns, Es, ".", label="Ordered spin")
plt.xlabel("Number of Monte Carlo Cycles", size=14)
plt.ylabel("Expected Energy", size=14)
plt.title("Mean energy per spin", size=18)
plt.legend(loc="best")
plt.show()

plt.plot(No, Mo, label="Ordered spin")
#plt.plot(Ns, Ms, ".", label="Ordered spin")
plt.plot(Nr, Mr, label="Random spin")
plt.xlabel("Number of Monte Carlo Cycles", size=14)
plt.ylabel("Magnetization", size=14)
plt.title("Mean magnetization per spin", size=18)
plt.legend(loc="center right")
plt.show()
