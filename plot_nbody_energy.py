import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

# donnees 
data = np.loadtxt("n_body_energy.txt")
t = data[:, 0]
E = data[:, 1]

# Affichage de l'energie totale
plt.figure(figsize=(10, 5))
plt.plot(t, E, label="Energie Totale", lw=2)
plt.xlabel("Temps (s)")
plt.ylabel("Energie (J)")
plt.title("Conservation de l'energie dans une simulation N-corps (RK4)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Charger les positions
try:
    sim_data = np.loadtxt("n_body_simulation.txt")
    t_pos = sim_data[:, 0]
    x1 = sim_data[:, 1]
    y1 = sim_data[:, 2]
    x2 = sim_data[:, 3]
    y2 = sim_data[:, 4]

    # Affichage des trajectoires
    plt.figure(figsize=(6, 6))
    plt.plot(x1, y1, label="Corps 1")
    plt.plot(x2, y2, label="Corps 2")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.title("Trajectoires des corps")
    plt.axis('equal')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Animation des trajectoires
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(min(np.min(x1), np.min(x2)), max(np.max(x1), np.max(x2)))
    ax.set_ylim(min(np.min(y1), np.min(y2)), max(np.max(y1), np.max(y2)))
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title("Animation des trajectoires")
    ax.grid(True)

    line1, = ax.plot([], [], 'b-', label="Corps 1")
    line2, = ax.plot([], [], 'r-', label="Corps 2")
    point1, = ax.plot([], [], 'bo')
    point2, = ax.plot([], [], 'ro')
    ax.legend()

    def init():
        line1.set_data([], [])
        line2.set_data([], [])
        point1.set_data([], [])
        point2.set_data([], [])
        return line1, line2, point1, point2

    def update(frame):
        line1.set_data(x1[:frame], y1[:frame])
        line2.set_data(x2[:frame], y2[:frame])
        point1.set_data([x1[frame]], [y1[frame]])
        point2.set_data([x2[frame]], [y2[frame]])
        return line1, line2, point1, point2

    ani = animation.FuncAnimation(fig, update, frames=len(t_pos), init_func=init, blit=True, interval=10)
    plt.tight_layout()
    plt.show()

except Exception as e:
    print("\n[INFO] Impossible d'afficher les trajectoires:", e)

