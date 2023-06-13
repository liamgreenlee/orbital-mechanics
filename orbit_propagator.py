import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.integrate import RK45

def newton_gravitation(t, y):
    # f = g*m1*m2/r**2
    # a = f/m2 = g*m1/r**2 = -mu/r**2
    position = y[:3]
    mu = 3.986004418e14
    a = -mu*position/np.linalg.norm(position)**3
    return np.hstack((y[3:], a))
    

if __name__ == "__main__":
    # Earth's radius is 6371000 m
    fig = plt.figure()
    ax = plt.axes(projection="3d")

    re = 6371000
    x = RK45(newton_gravitation, 0, np.array([2*re, 0, 0, 1200, 7000, 2200]), 10000000, max_step=1000)
    plt.plot(x.y[0], x.y[1], x.y[2], 'o', zorder=4, label='start')
    position = np.array([x.y])
    while x.t < 1000000:
        x.step()
        position = np.append(position, [x.y], axis=0)
    plt.plot(x.y[0], x.y[1], x.y[2], 'o', zorder=4, label='stop')

    u = np.linspace(0, 2 * np.pi, 15)
    v = np.linspace(0, np.pi, 15)
    x = re * np.outer(np.cos(u), np.sin(v))
    y = re * np.outer(np.sin(u), np.sin(v))
    z = re * np.outer(np.ones(np.size(u)), np.cos(v))

    ax.plot_surface(x, y, z, alpha=0.5)
    ax.plot(position[:, 0], position[:, 1], position[:, 2], 'red',zorder=3)
    ax.set_aspect("equal")
    plt.legend()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.grid()
    plt.show()

