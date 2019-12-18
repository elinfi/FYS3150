import matplotlib.pyplot as plt
import numpy as np
import scipy.special

def nl(l):
    return 1 + (3*l)/(1 - l)

def Pn(x, n):
    g = scipy.special.gamma(n)
    a = (n**n)/g

    P = a*x**(n-1)*np.exp(-n*x)
    return P


def parametrization(l):
    m0 = 500
    X = np.linspace(0, 2.5, 10000)

    P = np.zeros(len(X))
    n = nl(l)

    i = 0
    for x in X:
        P[i] = Pn(x, n)
        i += 1

    plt.plot(500*X, 28900*P)
    return X, P

def dist_curve(filename):

    data = np.loadtxt(filename)

    max = np.max(data)
    step = 30
    bins = int(max/step)
    m_val = np.linspace(0, max, bins)
    histo = np.histogram(data, bins=m_val)

    plt.plot(m_val[:-1],(histo[0]))
    plt.xlabel('money', size=14)
    plt.ylabel('frequency', size=14)


# Plot paramterization of distribution curves
plt.subplot(2, 1, 1)
for i in [0, 0.25, 0.5, 0.9]:
    parametrization(i)
plt.title("Parametrization of distribution curves", size=18)
plt.xlabel("$\\frac{m}{\langle m \\rangle}$", size=14)
plt.ylabel("$P(x)$", size=14)
plt.legend(["$\lambda = 0.0$", "$\lambda = 0.25$", "$\lambda = 0.5$", "$\lambda = 0.9$"])
plt.show()

# Plot distribution curves
plt.subplot(2, 1, 2)
dist_curve("5c_hist_1e7_1e3_N500_m500_l0_a0_g0.txt")
dist_curve("5c_hist_1e7_1e3_N500_m500_l025_a0_g0.txt")
dist_curve("5c_hist_1e7_1e3_N500_m500_l05_a0_g0.txt")
dist_curve("5c_hist_1e7_1e3_N500_m500_l09_a0_g0.txt")
plt.title("Distribution curves for the saving model", size=18)
plt.xlabel("Money")
plt.ylabel("Frequency")
plt.legend(['$\lambda = 0$','$\lambda = 0.25$', '$\lambda = 0.5$', '$\lambda = 0.9$'])
plt.show()
