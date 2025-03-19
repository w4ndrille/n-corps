# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Paramètres du cercle
radius = 5

# Création de la figure
fig, ax = plt.subplots()
ax.set_xlim(-radius - 1, radius + 1)
ax.set_ylim(-radius - 1, radius + 1)
ax.set_aspect('equal')

# Tracé du cercle
theta = np.linspace(0, 2*np.pi, 100)
x_circle = radius * np.cos(theta)
y_circle = radius * np.sin(theta)
ax.plot(x_circle, y_circle, 'b-')  # Cercle fixe

# Création du point mobile
point, = ax.plot([], [], 'ro', markersize=8)

# Fonction d'animation
def animate(i):
    angle = 2 * np.pi * i / 100  # Rotation progressive
    x = radius * np.cos(angle)
    y = radius * np.sin(angle)
    point.set_data([x], [y])  # Mettre les valeurs sous forme de liste
    return point,

# Création de l'animation
anim = FuncAnimation(fig, animate, frames=361, interval=20, blit=True)

plt.show()
