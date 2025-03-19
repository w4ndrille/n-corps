import numpy as np
import matplotlib.pyplot as plt

# Load the simulation data
data = np.loadtxt('n_body_simulation.txt')

# Extract time and positions
time = data[:, 0]
num_particles = (data.shape[1] - 1) // 2  # Number of particles
positions = data[:, 1:]  # Exclude the time column

# Plot trajectories of each particle
plt.figure(figsize=(8, 6))
for i in range(num_particles):
    x = positions[:, 2 * i]     # x-coordinates of particle i
    y = positions[:, 2 * i + 1] # y-coordinates of particle i
    plt.plot(x, y, label=f'Particle {i + 1}')

# Add labels, legend, and grid
plt.xlabel('x')
plt.ylabel('y')
plt.title('Trajectories of Bodies')
plt.legend()
plt.grid(True)

# Save and show the plot
plt.savefig('trajectories.png')
plt.show()