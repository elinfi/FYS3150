# This Python file uses the following encoding: utf-8

# if__name__ == "__main__":
#     pass

import numpy as np
import matplotlib.pyplot as plt

file = open("test.txt", "r")
N = 10

solution = np.array(N+1)

h = 1.0/(N + 1)
print(h)
#i = 0
#for line in file:
#    line.rstrip('\n')
#    print(line)
#    solution[i] = int(line)
#    i += 1

line = file.readline().split(",")
for i in range(0, len(line)):
    line[i] = float(line[i])
print(line)


x = h*np.linspace(0, 1, N)
plt.plot(x, line, label="N = "+str(N))
plt.xlabel(x)
plt.ylabel(solution)
plt.legend()
plt.show()
