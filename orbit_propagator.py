import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.integrate import RK45
from matplotlib import animation
from functools import partial
from typing import List

class StaticBody:
    def __init__(self, mass, position) -> None:
        self.mass = mass
        self.position = position
    
def newton_gravitation(t, y, bodies: List[StaticBody]):
    # f = g*m1*m2/r**2
    # a = f/m2 = g*m1/r**2 = -mu/r**2
    a = 0
    for body in bodies:
        position = y[:3] - body.position
        mu = body.mass
        mu = 3.986004418e5
        a += -mu*position/np.linalg.norm(position)**3
    return np.hstack((y[3:], a))

def animate_func(idx, x, bodies: List[StaticBody]):
    ax.clear()
    x.step()
    global position
    position = np.append(position, [x.y], axis=0)

    for body in bodies:
        u = np.linspace(0, 2 * np.pi, 15)
        v = np.linspace(0, np.pi, 15)
        x = body.position[0] + re * np.outer(np.cos(u), np.sin(v))
        y = body.position[1] + re * np.outer(np.sin(u), np.sin(v))
        z = body.position[2] + re * np.outer(np.ones(np.size(u)), np.cos(v))
        ax.plot_surface(x, y, z, alpha=0.5)

    ax.plot(position[:, 0], position[:, 1], position[:, 2], 'red', zorder=3)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.grid()

def orbital_elements_to_state_vector(h, e, i, O, w, theta):
    mu = 3.986004418e5
    r = (h**2/mu)*(1/(1 + e*np.cos(theta)))*np.array([np.cos(theta), np.sin(theta), 0])
    v = (mu/h)*np.array([-np.sin(theta), e + np.cos(theta), 0])
    Q = np.array([[np.cos(O)*np.cos(w) - np.sin(O)*np.sin(w)*np.cos(i), -np.cos(O)*np.sin(w) - np.sin(O)*np.cos(i)*np.cos(w), np.sin(O)*np.sin(i)],
                  [np.sin(O)*np.cos(w) + np.cos(O)*np.cos(i)*np.sin(w), -np.sin(O)*np.sin(w) + np.cos(O)*np.cos(i)*np.cos(w), -np.cos(O)*np.sin(i)],
                  [np.sin(i)*np.sin(w), np.sin(i)*np.cos(w), np.cos(i)]])
    r = np.dot(Q, r)
    v = np.dot(Q, v)
    return r, v


if __name__ == "__main__":
    # Earth's radius is 6371000 m
    bodies = [StaticBody(3.986004418e5, np.array([0, 0, 0])), StaticBody(3.986004418e5, np.array([-30000, 10000, 0]))]
    re = 6371
    r, v = orbital_elements_to_state_vector(80000, 0.6, np.pi*0/180, np.pi*0/180, np.pi*0/180, np.pi*0/180)
    x = RK45(partial(newton_gravitation, bodies=bodies), 0, np.append(r, v), 10000000, max_step=100)
    global position
    position = np.array([x.y])

    fig = plt.figure()
    ax = plt.axes(projection="3d")
    line_ani = animation.FuncAnimation(fig, partial(animate_func, x=x, bodies=bodies), interval=100)
    plt.show()