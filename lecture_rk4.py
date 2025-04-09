import numpy as np
import matplotlib.pyplot as plt

# Load the simulation data
data = np.loadtxt('n_body_simulation.txt')

# Extract time and positions
time = data[:500, 0]
num_particles = (data.shape[1] - 1) // 2  # Number of particles
positions = data[:, 1:]  # Exclude the time column

# Plot trajectories of each particle

for i in range(num_particles):
    x = positions[:500, 2 * i]     # x-coordinates of particle i
    y = positions[:500, 2 * i + 1] # y-coordinates of particle i
    # t = time[:]
    plt.plot(x, y, label=f'Corps {i + 1}', marker='o')



# Add labels, legend, and grid
plt.xlabel('position x (en mètres)')
plt.ylabel('position y (en mètres)')
plt.title('Trajectoires des corps')
plt.legend()
plt.grid(True)

# Save and show the plot
plt.savefig('trajectories_3bodies.png')
plt.show()