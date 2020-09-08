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

eps = 1e-1 # margin of error when checking if two quantities are equal


"""
Function used to set up the second order differential equations for the
x and y positions and velocities of the ball.

Parameters: r - array of the x and y positons and x and y velocities
Returns: array of the differential equation for the x and y positions and velocities
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
    v_x = r[2]
    v_y = r[3]

    fx = v_x
    fy = v_y

    f_v_x = coef*fx*np.sqrt(fx**2+fy**2)/(2*m)
    f_v_y = coef*fy*np.sqrt(fx**2+fy**2)/(2*m)-g

    return np.array([fx, fy, f_v_x, f_v_y], float)


"""
Function used to solve the differential equations, using the Fourth order
Runge-Kutta method.

Parameters: init - array of the initial x and y positions and velocities
Returns: solved differential equations of the x and y positions and velocities
"""
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


"""
Creates new arrays of the x and y positions and velocities of the basketball
until the time that it hits the backboard.

Parameters:
x_points, y_points, v_x_points, v_y_points - solved arrays for the x and y
                                            positions and velocities
T - time array for which the shot is being modelled
distance_backboard - float value, which will be changed for different shots

Returns: new_x, new_y, new_v_x, new_v_y - arrays of the original x and y
        positions or velocities of the shot until it hits the backboard
"""
def elastic(x_points, y_points, v_x_points, v_y_points, T, distance_backboard):
    C = 0.7493  # circumference in m, from basketball's circumference of 29.5 inches
    radius = C/(2*np.pi)
    t = T[1]-T[0]
    length = len(x_points) #this could have been any of the four arrays since they're the same length
    new_x = []
    new_y = []
    new_v_x = []
    new_v_y = []

    for i in range(length):
        if not (np.absolute(x_points[i]+radius-distance_backboard) <= eps and
            y_points[i] >= 3.048 and y_points[i] <= 4.1148):

            new_x.append(x_points[i])
            new_y.append(y_points[i])
            new_v_x.append(v_x_points[i])
            new_v_y.append(v_y_points[i])

        else:
            return new_x, new_y, new_v_x, new_v_y

"""
Determines whether a basket was made for the given shot.

Parameters:
x_points, y_points - x and y positions of the ball at all times
distance_backboard - starting distance from the backboard

Returns: boolean value indicating whether the ball went through the hoop's dimensions
"""
def in_basket(x_points, y_points, distance_backboard):
    C = 0.7493  # circumference in m, from basketball's circumference of 29.5 inches
    radius = C/(2*np.pi)
    shot_made = False

    for i in range(len(x_points)):
        if np.absolute(y_points[i]-3.048) <= eps:
            if ((x_points[i]-radius) >= (distance_backboard-0.6) and
                (x_points[i]+radius) <= (distance_backboard-0.15)):
                shot_made = True
                return shot_made

    return False


def backboard(x,y,v_x,v_y,T,distance_backboard):
    try:
        x_before, y_before, v_x_before, v_y_before = elastic(x, y, v_x, v_y, T, distance_backboard)
        backboard_hit = [x_before[-1],y_before[-1],-1*v_x_before[-1],v_y_before[-1]]
        x_after, y_after, v_x_after, v_y_after = RK4(backboard_hit)
        x = x_before + x_after
        y = y_before + y_after
        v_x = v_x_before + v_x_after
        v_y = v_y_before + v_y_after

    except:
        print("Ball does not come in contact with the backboard.")

    print(in_basket(x,y,distance_backboard))
    plt.plot(x,y)
    plt.plot(np.linspace(distance_backboard-0.6,distance_backboard-0.15,100),[3.048]*100)
    plt.plot([distance_backboard]*100,np.linspace(3.048,4.1148,100))
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
    distance_backboard = 2.5 #arbitrary distance from backboard in m


    init = [0, start_height, v0*np.cos(theta), v0*np.sin(theta)]
    x_points, y_points, v_x_points, v_y_points = RK4(init)

    backboard(x_points,y_points,v_x_points,v_y_points,T,distance_backboard)




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
    #plt.plot(x_points, y_points)


if __name__ == "__main__":
    main()
