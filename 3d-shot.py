"""
file: 3d-shot.py

Author: Caitlin Welch
Date created: August 31, 2020
Date modified: August 31, 2020

Brief: Uses RK4 ODE solver to plot the 3D trajectory of a Wilson basketball
"""

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt


def f(r):
    C_d = 0.47  # drag coefficient for a sphere
    rho = 1.225  # air density in kg/m^3
    C = 0.7493  # circumference in m, from basketball's circumference of 29.5 inches
    radius = C/(2*np.pi)
    A = C**2/(4*np.pi)  # cross-sectional area in m^2
    g = 9.81  # gravitational constant on Earth
    m = 22/35.274  # gives mass in kg of 22oz basketball

    coef = -C_d*rho*A  # only used so I don't have to keep writing this out
    x = r[0]
    y = r[1]
    z = r[2]
    v_x = r[3]
    v_y = r[4]
    v_z = r[5]

    fx = v_x
    fy = v_y
    fz = v_z

    f_v_x = coef*fx*np.sqrt(fx**2+fy**2+fz**2)/(2*m)
    f_v_y = coef*fy*np.sqrt(fx**2+fy**2+fz**2)/(2*m)
    f_v_z = coef*fz*np.sqrt(fx**2+fy**2+fz**2)/(2*m)-g

    return np.array([fx, fy, fz, f_v_x, f_v_y, f_v_z], float)


def RK4(init):
    x_points = []
    y_points = []
    z_points = []
    v_x_points = []
    v_y_points = []
    v_z_points = []

    r = np.array(init, float)
    h = 0.01

    while r[2] >= 0:
        x_points.append(r[0])
        y_points.append(r[1])
        z_points.append(r[2])
        v_x_points.append(r[3])
        v_y_points.append(r[4])
        v_z_points.append(r[5])

        k1 = h*f(r)
        k2 = h*f(r+0.5*k1)
        k3 = h*f(r+0.5*k2)
        k4 = h*f(r+k3)
        r += (k1+2*k2+2*k3+k4)/6

    return x_points, y_points, z_points, v_x_points, v_y_points, v_z_points

def main():
    """
    Constants for linspace
    """
    a = 0
    b = 1000
    N = 10000
    T = np.linspace(a,b,N)
    h = (b-a)/N

    """
    Change each of these initial values for different shots when doing the Monte Carlo simulation
    """
    angle = 38  # float(input("Enter angle shot at in degrees: "))
    theta = (np.pi/180) * angle  # float(input("Enter angle shot at in degrees: "))
    v0 = 9.8  # m/s # float(input("Enter initial velocity (m/s): "))
    start_height = 1.8 # m - this for someone around 6ft tall
    # postive x values represent a shot coming from the right of the hoop, negative x values represent shots coming from the left of the hoop
    start_x = 3
    start_y = 4 # all y values will be positive - distance from hoop, perpendicular to endline
    phi = np.arctan(start_y/start_x)

    init = [start_x, start_y, start_height, v0*np.cos(theta)*np.cos(phi), v0*np.cos(theta)*np.sin(phi), v0*np.sin(theta)]
    x_points, y_points, z_points, v_x_points, v_y_points, v_z_points = RK4(init)

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(x_points, y_points, z_points)
    plt.show()


if __name__ == "__main__":
    main()
