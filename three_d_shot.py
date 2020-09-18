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


eps_backboard = 1e-2 #margin of error when checking if two quantities are equal
eps_rim = 1e-1

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

    fx = -v_x
    fy = -v_y
    fz = v_z

    f_v_x = -coef*fx*np.sqrt(fx**2+fy**2+fz**2)/(2*m)
    f_v_y = -coef*fy*np.sqrt(fx**2+fy**2+fz**2)/(2*m)
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
for the motion of the basketball when the backboard is in play. Uses the values
of each position and velocity right when the ball hits the backboard to solve
for the positions and velocities of the ball after it hits the backboard.

Parameters:
x_points, y_points, z_points, v_x_points, v_y_points v_z_points - solved arrays
    for the x, y and positions and velocities

Returns: x,y,z,v_x,v_y,v_z - arrays of each respective variable including both
        before and after the ball hits the backboard
"""
def backboard(x_points, y_points, z_points, v_x_points, v_y_points, v_z_points):
    C = 0.7493  # circumference in m, from basketball's circumference of 29.5 inches
    radius = C/(2*np.pi)
    hits_backboard = False

    index = -1
    for i in range(len(y_points)):
        #how to extend lines 114-115 to 3D code
        if np.absolute(y_points[i]-radius) <= eps:
            index = i

    if z_points[index] >= 3.048 and z_points[index] <= 4.1148 and x_points[index] <= 0.9144 and x_points[index] >= -0.9144:
        hits_backboard = True

    if hits_backboard:
        print("hit")
        x_before = x_points[:index]
        y_before = y_points[:index]
        z_before = z_points[:index]
        v_x_before = v_x_points[:index]
        v_y_before = v_y_points[:index]
        v_z_before = v_z_points[:index]

        backboard_hit = [x_before[-1],y_before[-1],z_before[-1],v_x_before[-1],-1*v_y_before[-1],v_z_before[-1]]
        x_after,y_after,z_after,v_x_after,v_y_after,v_z_after = RK4(backboard_hit)

        x = x_before + x_after
        y = y_before + y_after
        z = z_before + z_after
        v_x = v_x_before + v_x_after
        v_y = v_y_before + v_y_after
        v_z = v_z_before + v_z_after
    else:
        print("Ball does not come in contact with the backboard.")
        x = x_points
        y = y_points
        z = z_points
        v_x = v_x_points
        v_y = v_y_points
        v_z = v_z_points

    return x,y,z,v_x,v_y,v_z


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
                if ((x_points[i]-radius) >= -0.225) and (x_points[i]+radius <= 0.225):
                    shot_made = True
                    return shot_made

    return False


def main():
    """
    Constants for linspace
    """
    a = 0
    b = 1000
    N = 1000000
    T = np.linspace(a,b,N)
    h = (b-a)/N

    """
    Change each of these initial values for different shots when doing the Monte Carlo simulation
    """
    angle = 80  # float(input("Enter angle shot at in degrees: "))
    theta = (np.pi/180) * angle  # float(input("Enter angle shot at in degrees: "))
    v0 = 9.8  # m/s # float(input("Enter initial velocity (m/s): "))
    start_height = 1.8 # m - this is for someone around 6ft tall
    # postive x values represent a shot coming from the right of the hoop, negative x values represent shots coming from the left of the hoop
    start_x = 0
    start_y = 2.5 # all y values will be positive - distance from hoop, perpendicular to endline
    phi = np.arctan(start_x/start_y)

    init = [start_x, start_y, start_height, v0*np.cos(theta)*np.sin(phi), v0*np.cos(theta)*np.cos(phi), v0*np.sin(theta)]
    x_points, y_points, z_points, v_x_points, v_y_points, v_z_points = RK4(init)

    x_backboard,y_backboard,z_backboard,v_x_backboard,v_y_backboard,v_z_backboard = backboard(x_points,y_points,z_points,v_x_points,v_y_points,v_z_points)
    print(in_basket(x_backboard,y_backboard,z_backboard))

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(x_backboard, y_backboard, z_backboard)

    ax.plot3D(np.linspace(-0.9144,0.9144,1000),[0]*1000,[3.048]*1000,'r')
    ax.plot3D(np.linspace(-0.9144,0.9144,1000),[0]*1000,[4.1148]*1000,'r')
    ax.plot3D([-0.9144]*1000,[0]*1000,np.linspace(3.048,4.1148,1000),'r')
    ax.plot3D([0.9144]*1000,[0]*1000,np.linspace(3.048,4.1148,1000),'r')
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    circle_x = np.linspace(-0.225,0.225,1000)
    circle_y1 = []
    circle_y2 = []
    for x in circle_x:
        circle_y1.append(np.sqrt(.225**2-x**2)+.375)
        circle_y2.append(-np.sqrt(.225**2-x**2)+.375)
    ax.plot3D(circle_x,circle_y1,[3.048]*1000,'g')
    ax.plot3D(circle_x,circle_y2,[3.048]*1000,'g')
    plt.show()


if __name__ == "__main__":
    main()
