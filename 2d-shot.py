"""
file: 2d-shot.py

Author: Caitlin Welch
Date created: August 24, 2020
Date modified: August 29, 2020

Brief: Uses RK4 ODE solver to plot the 2D trajectory of a Wilson basketball
"""

import numpy as np
import pylab as plt
from scipy.integrate import odeint






def f(r):

    C_d = 0.47  # drag coefficient for a sphere
    rho = 1.225  # air density in kg/m^3
    C = 0.7493  # circumference in m, from basketball's circumference of 29.5 inches
    A = C**2/(4*np.pi)  # cross-sectional area in m^2
    g = 9.81  # gravitational constant on Earth
    m = 22/35.274  # gives mass in kg of 22oz basketball

    coef = -C_d*rho*A  # only used so I don't have to keep writing this out
    x = r[0]
    y = r[1]
    v_x = r[2]
    v_y = r[3]

    fx = v_x
    fy = v_y
    f_v_x = coef*fx*np.sqrt(fx**2+fy**2)/(2*m)
    f_v_y = coef*fy*np.sqrt(fx**2+fy**2)/(2*m)-g

    return np.array([fx, fy, f_v_x, f_v_y], float)


def RK4(init):
    x_points = []
    y_points = []
    v_x_points = []
    v_y_points = []

    r = np.array(init, float)
    h = 0.01
    while r[1] >= 0:
        x_points.append(r[0])
        y_points.append(r[1])
        v_x_points.append(r[2])
        v_y_points.append(r[3])

        k1 = h*f(r)
        k2 = h*f(r+0.5*k1)
        k3 = h*f(r+0.5*k2)
        k4 = h*f(r+k3)
        r += (k1+2*k2+2*k3+k4)/6

    return x_points, y_points, v_x_points, v_y_points


def main():

    """
    Constants for linspace
    """
    a = 0
    b = 1000
    N = 10000
    T = np.linspace(a, b, N)
    h = (b-a)/N

    """
    Change each of these initial values for different shots when doing the Monte Carlo simulation
    """
    angle = 38  # float(input("Enter angle shot at in degrees: "))
    v0 = 9.8  # m/s # float(input("Enter initial velocity (m/s): "))
    theta = (np.pi/180) * angle  # float(input("Enter angle shot at in degrees: "))
    start_height = 1.8
    init = [0, start_height, v0*np.cos(theta), v0*np.sin(theta)]
    x_points, y_points, v_x_points, v_y_points = RK4(init)

    # if using odeint
    #solution = odeint(f, init, T)

    # shortens the array for only positive y values
    """
    positive_y = []
    for ys in solution[:,1]:
        if ys >= 0:
            positive_y.append(ys)
    """

    # only plots until the basketball hits the ground
    #plt.plot(solution[:len(positive_y), 0], positive_y)
    plt.plot(x_points,y_points)
    plt.show()


if __name__ == "__main__":
    main()
