import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.integrate import RK45
from matplotlib import animation
from functools import partial

def newton_gravitation(t, y):
    # f = g*m1*m2/r**2
    # a = f/m2 = g*m1/r**2 = -mu/r**2
    position = y[:3]
    mu = 3.986004418e5
    a = -mu*position/np.linalg.norm(position)**3
    return np.hstack((y[3:], a))

def animate_func(idx, x):
    ax.clear()
    x.step()
    global position
    position = np.append(position, [x.y], axis=0)

    u = np.linspace(0, 2 * np.pi, 15)
    v = np.linspace(0, np.pi, 15)
    x = re * np.outer(np.cos(u), np.sin(v))
    y = re * np.outer(np.sin(u), np.sin(v))
    z = re * np.outer(np.ones(np.size(u)), np.cos(v))

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
    re = 6371
    r, v = orbital_elements_to_state_vector(80000, 0.5, np.pi*30/180, np.pi*40/180, np.pi*60/180, np.pi*30/180)
    x = RK45(newton_gravitation, 0, np.append(r, v), 10000000, max_step=100)
    global position
    position = np.array([x.y])

    fig = plt.figure()
    ax = plt.axes(projection="3d")
    line_ani = animation.FuncAnimation(fig, partial(animate_func, x=x), interval=100)
    plt.show()