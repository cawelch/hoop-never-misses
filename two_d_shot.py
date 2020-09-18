"""
file: 2d-shot.py

Author: Caitlin Welch
Date created: August 24, 2020
Date modified: August 31, 2020

Brief: Uses RK4 ODE solver to plot the 2D trajectory of a Wilson basketball
"""

import numpy as np
import pylab as plt
from scipy.integrate import odeint

eps_backboard = 1e-2 # margin of error when checking if two quantities are equal
eps_rim = 1e-1


"""
Function used to set up the second order differential equations for the
y and z positions and velocities of the ball.

Parameters: r - array of the y and z positons and y and z velocities
Returns: array of the differential equation for the y and z positions and velocities
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
    y = r[0]
    z = r[1]
    v_y = r[2]
    v_z = r[3]

    fy = -v_y
    fz = v_z

    f_v_y = -coef*fy*np.sqrt(fy**2+fz**2)/(2*m)
    f_v_z = coef*fz*np.sqrt(fy**2+fz**2)/(2*m)-g

    return np.array([fy, fz, f_v_y, f_v_z], float)


"""
Function used to solve the differential equations, using the Fourth order
Runge-Kutta method.

Parameters: init - array of the initial y and z positions and velocities
Returns: solved differential equations of the y and z positions and velocities
"""
def RK4(init):
    y_points = []
    z_points = []
    v_y_points = []
    v_z_points = []

    r = np.array(init, float)
    h = 0.01
    while r[1] >= 0:
        y_points.append(r[0])
        z_points.append(r[1])
        v_y_points.append(r[2])
        v_z_points.append(r[3])

        k1 = h*f(r)
        k2 = h*f(r+0.5*k1)
        k3 = h*f(r+0.5*k2)
        k4 = h*f(r+k3)
        r += (k1+2*k2+2*k3+k4)/6

    return y_points, z_points, v_y_points, v_z_points


"""
Creates new arrays of the y and z positions and velocities of the basketball
until the time that it hits the backboard.

Parameters:
y_points, z_points, v_y_points, v_z_points - solved arrays for the y and z
                                            positions and velocities
T - time array for which the shot is being modelled
distance_backboard - float value, which will be changed for different shots

Returns: new_y, new_z, new_v_y, new_v_z - arrays of the original y and z
        positions or velocities of the shot until it hits the backboard
"""
def elastic(y_points, z_points, v_y_points, v_z_points, T):
    C = 0.7493  # circumference in m, from basketball's circumference of 29.5 inches
    radius = C/(2*np.pi)
    t = T[1]-T[0]
    length = len(y_points) #this could have been any of the four arrays since they're the same length
    new_y = []
    new_z = []
    new_v_y = []
    new_v_z = []

    for i in range(length):
        y = 0
        z_max = 4.1148
        z_min = 3.048

        #checks if the ball comes within the tolerance of the backboard in the
        #y directino
        if not (np.absolute(y_points[i]) <= eps_backboard):
            #checks if any point on the circumference of the circle hits the z range
            #of the backboard
            try:
                if not (z_points[i]+np.sqrt(radius**2-(y-y_points[i])**2) >= z_min and
                    z_points[i]+np.sqrt(radius**2-(y-y_points[i])**2) <= z_max):
                    
                    new_y.append(y_points[i])
                    new_z.append(z_points[i])
                    new_v_y.append(v_y_points[i])
                    new_v_z.append(v_z_points[i])
            except:
                #if the square root can't be evaluated, then the y value is too
                #large in comparison to the radius, so it won't hit the backboard
                new_y.append(y_points[i])
                new_z.append(z_points[i])
                new_v_y.append(v_y_points[i])
                new_v_z.append(v_z_points[i])

        else:
            return new_y, new_z, new_v_y, new_v_z

"""
Determines whether a basket was made for the given shot.

Parameters:
y_points, z_points - y and z positions of the ball at all times
distance_backboard - starting distance from the backboard

Returns: boolean value indicating whether the ball went through the hoop's dimensions
"""
def in_basket(y_points, z_points):
    C = 0.7493  # circumference in m, from basketball's circumference of 29.5 inches
    radius = C/(2*np.pi)
    shot_made = False

    for i in range(len(y_points)):
        #here I have whether the centre of the ball goes through the hoop
        if np.absolute(z_points[i]-3.048) <= eps_rim:
            if ((y_points[i]+radius) <= (0.6) and (y_points[i]-radius) >= (0.15)):
                shot_made = True
                return shot_made

    return False


def backboard(y,z,v_y,v_z,T):
    C = 0.7493  # circumference in m, from basketball's circumference of 29.5 inches
    radius = C/(2*np.pi)
    try:
        y_before, z_before, v_y_before, v_z_before = elastic(y, z, v_y, v_z, T)
        backboard_hit = [y_before[-1],z_before[-1],-1*v_y_before[-1],v_z_before[-1]]
        y_after, z_after, v_y_after, v_z_after = RK4(backboard_hit)
        y = y_before + y_after
        z = z_before + z_after
        v_y = v_y_before + v_y_after
        v_z = v_z_before + v_z_after
    except:
        print("Ball does not come in contact with the backboard.")

    print(in_basket(y,z))

    """
    Plots the trajectory, backboard, tolerance levels for backboard, hoop and
    area of the hoop that the centre of the basketball can go through for the
    whole basketball to make it.
    """
    plt.plot(y,z,'.')
    plt.plot(np.linspace(0.6,0.15,100),[3.048]*100)
    plt.plot([0.15+radius]*100,np.linspace(3.148,2.948,100))
    plt.plot([0.6-radius]*100,np.linspace(3.148,2.948,100))

    plt.plot([0]*100,np.linspace(3.048,4.1148,100))
    plt.plot([-.01]*100,np.linspace(3.048,4.1148,100))
    plt.plot([.01]*100,np.linspace(3.048,4.1148,100))

    plt.axis("equal")
    plt.show()


def main():
    """
    Constants for linspace
    """
    a = 0
    b = 1000
    N = 10000000
    T = np.linspace(a, b, N)
    h = (b-a)/N

    """
    Change each of these initial values for different shots when doing the Monte Carlo simulation
    """
    angle = 80  # float(input("Enter angle shot at in degrees: "))
    theta = (np.pi/180) * angle  # float(input("Enter angle shot at in degrees: "))
    v0 = 9.8  # m/s # float(input("Enter initial velocity (m/s): "))
    start_height = 1.8
    start_y = 2.5 #arbitrary distance from backboard in m


    init = [start_y, start_height, v0*np.cos(theta), v0*np.sin(theta)]
    y_points, z_points, v_y_points, v_z_points = RK4(init)

    backboard(y_points,z_points,v_y_points,v_z_points,T)




    # if using odeint
    #solution = odeint(f, init, T)

    # shortens the array for only positive y values
    """
    positive_z = []
    for zs in solution[:,1]:
        if zs >= 0:
            positive_z.append(ys)
    """



    # only plots until the basketball hits the ground
    #plt.plot(solution[:len(positive_z), 0], positive_z)
    #plt.plot(y_points, z_points)


if __name__ == "__main__":
    main()
