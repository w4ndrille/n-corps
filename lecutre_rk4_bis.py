import numpy as np
import matplotlib.pyplot as plt


x_numérique = np.loadtxt('write_results_q10.txt')

# Solution analytique
x0 = 1
omega0 = 0.5
x_analytique = x0 * np.cos(omega0 * x_numérique[:,0])

x = x_numérique[:,1] - x_analytique

# Tracé de la figure
plt.plot(x_numérique[:,0], x_numérique[:,1], label='Solution numérique')
plt.plot(x_numérique[:,0], x_analytique, label='Solution analytique', linestyle='--')
plt.legend()
plt.grid(True)

plt.savefig('q09.png')
plt.show()
