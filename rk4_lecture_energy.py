import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("n_body_energy.txt")
t = data[:, 0]
E = data[:, 1]

plt.plot(t, E)
plt.xlabel("Time (s)")
plt.ylabel("Total Energy (J)")
plt.title("Conservation de l'énergie")
plt.grid()
plt.show()
