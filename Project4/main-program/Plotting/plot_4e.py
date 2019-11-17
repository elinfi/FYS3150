import numpy as np
import matplotlib.pyplot as plt

def read_file(filename):
    infile = open(filename, "r")

    n = int(infile.readline())

    temp = np.empty(n)
    E = np.empty(n)
    M_abs = np.empty(n)
    M = np.empty(n)
    Cv = np.empty(n)
    Chi = np.empty(n)
    Chi_abs = np.empty(n)

    i = 0
    for line in infile:
        t, e, m_abs, m, cv, chi, chi_abs = line.split()
        temp[i] = float(t)
        E[i] = float(e)
        M_abs[i] = float(m_abs)
        M[i] = float(m)
        Cv[i] = float(cv)*(float(t)**2)
        Chi[i] = float(chi)*float(t)**3
        Chi_abs[i] = float(chi_abs)*float(t)**3

        i += 1

    return temp, E, M_abs, M, Cv, Chi, Chi_abs


def read_file2(filename):
    infile = open(filename, "r")

    n = int(infile.readline())

    temp = np.empty(n)
    E = np.empty(n)
    M_abs = np.empty(n)
    Cv = np.empty(n)
    Chi = np.empty(n)

    i = 0
    for line in infile:
        t, e, m_abs, cv, chi = line.split()
        temp[i] = float(t)
        E[i] = float(e)
        M_abs[i] = float(m_abs)
        Cv[i] = float(cv)*(float(t)**2)
        Chi[i] = float(chi)*(float(t)**3)

        i += 1

    return temp, E, M_abs, Cv, Chi

def read_file3(filename):
    infile = open(filename, "r")

    n = int(infile.readline())

    temp = np.empty(n)
    E = np.empty(n)
    M_abs = np.empty(n)
    M = np.empty(n)
    Cv = np.empty(n)
    Chi = np.empty(n)
    Chi_abs = np.empty(n)

    i = 0
    for line in infile:
        t, e, m_abs, m, cv, chi, chi_abs = line.split()
        temp[i] = float(t)
        E[i] = float(e)
        M_abs[i] = float(m_abs)
        M[i] = float(m)
        Cv[i] = float(cv)*(float(t)**2)
        Chi[i] = float(chi)*float(t)**2
        Chi_abs[i] = float(chi_abs)*float(t)**2

        i += 1

    return temp, E, M_abs, M, Cv, Chi, Chi_abs


T_40, E_40, M_abs_40, M_40, Cv_40, Chi_40, Chi_abs_40 = read_file("4e_L40_N5e5_21_24_005.txt")
T_60, E_60, M_abs_60, M_60, Cv_60, Chi_60, Chi_abs_60 = read_file("4e_L60_N5e5_21_24_005.txt")
T_80, E_80, M_abs_80, M_80, Cv_80, Chi_80, Chi_abs_80 = read_file("4e_L80_N5e5_21_24_005.txt")
T_100, E_100, M_abs_100, M_100, Cv_100, Chi_100, Chi_abs_100 = read_file3("4e_L100_N5e5_21_24_005.txt")

T40, E40, Mabs_40, Cv40, Chi40 = read_file2("4e_L40_N5e5_21_24_001.txt")

index_40, = np.where(Cv_40 == np.max(Cv_40))[0]
index40, = np.where(Cv40 == np.max(Cv40))[0]
index_60, = np.where(Cv_60 == np.max(Cv_60))[0]
index_80, = np.where(Cv_80 == np.max(Cv_80))[0]
index_100, = np.where(Cv_100 == np.max(Cv_100))[0]

Tc_40 = T_40[index_40]
Tc40 = T40[index40]
Tc_60 = T_60[index_60]
Tc_80 = T_80[index_80]
Tc_100 = T_100[index_100]

plt.plot(T_40, E_40, ".", label="L = 40")
plt.plot(T_60, E_60, ".", label="L = 60")
plt.plot(T_80, E_80, ".", label="L = 80")
plt.plot(T_100, E_100, ".", label="L = 100")
plt.xlabel("Temperature", size=14)
plt.ylabel("Expected Energy", size=14)
plt.title("Expected Energy per Spin", size=18)
plt.legend()
plt.show()

plt.plot(T_40, M_abs_40, ".", label = "L = 40")
plt.plot(T_60, M_abs_60, ".", label = "L = 60")
plt.plot(T_80, M_abs_80, ".", label = "L = 80")
plt.plot(T_100, M_abs_100, ".", label = "L = 100")
plt.xlabel("Temperature", size=14)
plt.ylabel("Absolute Magnetic Moment", size=14)
plt.title("Absolute Magnetic Moment per Spin", size=18)
plt.legend()
plt.show()

plt.plot(T_40, M_40, ".", label = "L = 40")
plt.plot(T_60, M_60, ".", label = "L = 60")
plt.plot(T_80, M_80, ".", label = "L = 80")
plt.plot(T_100, M_100, ".", label = "L = 100")
plt.xlabel("Temperature", size=14)
plt.ylabel("Magnetization", size=14)
plt.title("Magnetization per Spin", size=18)
plt.legend()
plt.show()

plt.plot(T_40, Cv_40, ".", label = "L = 40")
plt.plot(T_60, Cv_60, ".", label = "L = 60")
plt.plot(T_80, Cv_80, ".", label = "L = 80")
plt.plot(T_100, Cv_100, ".", label = "L = 100")
plt.plot(T_40[index_40], Cv_40[index_40], "x", color="b", label = "Tc_40 = " + str(Tc_40))
plt.plot(T_60[index_60], Cv_60[index_60], "x", color="orange", label = "Tc_60 = " + str(Tc_60))
plt.plot(T_80[index_80], Cv_80[index_80], "x", color="g", label = "Tc_80 = " + str(Tc_80))
plt.plot(T_100[index_100], Cv_100[index_100], "x", color="r", label = "Tc_100 = " + str(Tc_100))
plt.xlabel("Temperature", size=14)
plt.ylabel("Heat Capacity", size=14)
plt.title("Heat Capacity per Spin", size=18)
plt.legend()
plt.show()

plt.plot(T_40, Chi_40, ".", label = "L = 40")
plt.plot(T_60, Chi_60, ".", label = "L = 60")
plt.plot(T_80, Chi_80, ".", label = "L = 80")
plt.plot(T_100, Chi_100, ".", label = "L = 100")
plt.xlabel("Temperature", size=14)
plt.ylabel("Susceptibility", size=14)
plt.title("Susceptibility per Spin", size=18)
plt.legend()
#plt.show()

plt.plot(T_40, Chi_abs_40, ".", label = "L = 40")
plt.plot(T_60, Chi_abs_60, ".", label = "L = 60")
plt.plot(T_80, Chi_abs_80, ".", label = "L = 80")
plt.plot(T_100, Chi_abs_100, ".", label = "L = 100")
plt.xlabel("Temperature", size=14)
plt.ylabel("Susceptibility", size=14)
plt.title("Susceptibility per Spin", size=18)
plt.legend()
#plt.show()


a1 = (Tc_60 - Tc_40)/((1/60) - (1/40))
a2 = (Tc_80 - Tc_60)/((1/80) - (1/60))
a3 = (Tc_100 - Tc_80)/((1/100) - (1/80))
Tc1 = Tc_60 - a1*(1/60)
Tc2 = Tc_80 - a2*(1/80)
Tc3 = Tc_100 - a3*(1/100)
exact = 2.269
print(Tc1)
print(abs(Tc1 - exact)/exact)
print(Tc2)
print(abs(Tc2 - exact)/exact)
print(Tc3)
print(abs(Tc3 - exact)/exact)



"""
Muligens ved en temperature slutter den å flippe.
Bedre numerisk å bruke absolutt magnetisering, fordi den må kjøre veldig lenge for å nå 0. Magnetiseringer er bare i toppen.
Det tar veldig lang tid for magnetiseringen konverger, så det gir bedre resultater å kjøre med den absolutte magnetiseringen.
"""
