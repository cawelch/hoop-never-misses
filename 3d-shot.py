"""
file: 3d-shot.py

Author: Caitlin Welch
Date created: August 31, 2020
Date modified: August 31, 2020

Brief: Uses RK4 ODE solver to plot the 3D trajectory of a Wilson basketball
"""

#from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

eps = 1e-2 #margin of error when checking if two quantities are equal

"""
Function used to set up the second order differential equations for the
x, y and z positions and velocities of the ball.

Parameters: r - array of the x and y positons and x, y and z velocities
Returns: array of the differential equation for the x, y and z positions and velocities.
"""
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


"""
Function used to solve the differential equations, using the Fourth order
Runge-Kutta method.

Parameters: init - array of the initial x, y and z positions and velocities
Returns: solved differential equations of the x, y and z positions and velocities
"""
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


"""
Creates new arrays of the x, y and z positions and velocities of the basketball
until the time that it hits the backboard.

Parameters:
x_points, y_points, z_points, v_x_points, v_y_points v_z_points - solved arrays
    for the x, y and positions and velocities
T - time array for which the shot is being modelled
distance_backboard - float value, which will be changed for different shots

Returns: new_x, new_y, new_z, new_v_x, new_v_y new_v_z - arrays of the original
        x and y positions or velocities of the shot until it hits the backboard
"""
def elastic(x_points, y_points, z_points, v_x_points, v_y_points, v_z_points, T):
    t = T[1]-T[0]
    length = len(x_points)
    new_x = []
    new_y = []
    new_z = []
    new_v_x = []
    new_v_y = []
    new_v_z = []

    for i in range(length):
        if not ((z_points[i] >= 3.048 and z_points[i] <= 4.1148) and
            y_points[i] <= eps and (x_points[i] <= 0.9144 and x_points[i] >= -0.9144)):

            new_x.append(x_points[i])
            new_y.append(y_points[i])
            new_z.append(z_points[i])
            new_v_x.append(v_x_points[i])
            new_v_y.append(v_y_points[i])
            new_v_z.append(v_z_points[i])

        else:
            return new_x, new_y, new_z, new_v_x, new_v_y, new_v_z


"""
Determines whether a basket was made for the given shot.

Parameters: x_points, y_points, z_points - x, y and z positions of the ball at all times

Returns: boolean value indicating whether the ball went through the hoop's dimensions
"""
def in_basket(x_points, y_points, z_points):
    C = 0.7493  # circumference in m, from basketball's circumference of 29.5 inches
    radius = C/(2*np.pi)
    shot_made = False

    for i in range(len(x_points)):
        if np.absolute(z_points[i]-3.048) <= eps:
            if ((y_points[i]-radius) >= 0.15 and
                (y_points[i]+radius) <= 0.6):
                if ((x_points[i]-radius) >= -0.2286) and (x_points[i]+radius <= 0.2286):
                    shot_made = True
                    return shot_made

    return False


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
    angle = 80  # float(input("Enter angle shot at in degrees: "))
    theta = (np.pi/180) * angle  # float(input("Enter angle shot at in degrees: "))
    v0 = 9.8  # m/s # float(input("Enter initial velocity (m/s): "))
    start_height = 1.8 # m - this for someone around 6ft tall
    # postive x values represent a shot coming from the right of the hoop, negative x values represent shots coming from the left of the hoop
    start_x = 0
    start_y = 2.5 # all y values will be positive - distance from hoop, perpendicular to endline
    phi = np.arctan(start_x/start_y)

    init = [start_x, start_y, start_height, v0*np.cos(theta)*np.sin(phi), v0*np.cos(theta)*np.cos(phi), v0*np.sin(theta)]
    x_points, y_points, z_points, v_x_points, v_y_points, v_z_points = RK4(init)

    try:
        x_before, y_before, z_before, v_x_before, v_y_before, v_z_before = elastic(x_points, y_points, z_points, v_x_points, v_y_points, v_z_points,T)
        backboard_hit = [x_before[-1],y_before[-1],z_before[-1],v_x_before[-1],-1*v_y_before[-1],v_z_before[-1]]
        x_after, y_after, z_after, v_x_after, v_y_after, v_z_after = RK4(backboard_hit)
        x_backboard = x_before + x_after
        y_backboard = y_before + y_after
        z_backboard = z_before + z_after
        v_x_backboard = v_x_before + v_x_after
        v_y_backboard = v_y_before + v_y_after
        v_z_backboard = v_z_before + v_z_after

        print(in_basket(x_backboard,y_backboard,z_backboard))
        #PLOT 3D BACKBOARD SHOTS HERE

    except:
        print("Ball does not come in contact with the backboard.")
        print(in_basket(x_points,y_points,z_points))
        #PLOT 3D SHOTS HERE


    """
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(x_points, y_points, z_points)
    plt.show()
    """

if __name__ == "__main__":
    main()
