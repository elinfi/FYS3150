import numpy as np
import matplotlib.pyplot as plt

def hist_diff(filename):

    data = np.loadtxt(filename)

    max = np.max(data)
    step = 30
    bins = int(max/step)
    m_val = np.linspace(0, max, bins)

    hist_arr = np.zeros((len(data[0]), len(m_val)-1))
    diff = np.zeros(len(data[:-2,1]))
    hist_arr[0] = np.histogram(data[0], bins = m_val)[0]
    for i in range(2, len(data[0])):

        hist_arr[i,:] = np.histogram(data[:i], bins=m_val)[0]
        diff[i-2]=np.linalg.norm(((hist_arr[i]/i)-(hist_arr[i-1]/(i-1))))
    print(diff[600])

    plt.plot(diff)


hist_diff("5c_hist_1e7_1e3_N500_m500_l0_a0_g0.txt")
hist_diff("5c_hist_1e7_1e3_N500_m500_l025_a0_g0.txt")
hist_diff("5c_hist_1e7_1e3_N500_m500_l05_a0_g0.txt")
hist_diff("5c_hist_1e7_1e3_N500_m500_l09_a0_g0.txt")
plt.title('Relative difference in Monte Carlo steps', size=18)
plt.xlabel('Monte Carlo steps', size=14)
plt.ylabel('Relative difference', size=14)
plt.legend(["$\lambda = 0.9", "$\lambda = 0.25$", "$\lambda = 0.5$", "$\lambda = 0.9$"])
plt.show()


hist_diff("5d_hist_1e7_1e4_N1000_m500_l09_a2_g0.txt")
plt.title('Relative difference in Monte Carlo steps', size=18)
plt.xlabel('Monte Carlo steps', size=14)
plt.ylabel('Relative difference', size=14)
plt.legend(["$\lambda = 0.9, \\alpha = 2, N = 1000$"])
plt.show()
