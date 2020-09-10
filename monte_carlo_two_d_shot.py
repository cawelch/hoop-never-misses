"""
file: monte-carlo-2d-shot.py

Author: Caitlin Welch
Date created: September 6, 2020
Date modified: September 9, 2020

Brief: First attempt at running a Monte Carlo simulation for a 2D shot to form
        a continuous curved backboard
"""

import two_d_shot
import numpy as np
import pylab as plt

eps = 1e-1 #margin of error when checking if two quantities are equal


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

    fy = v_y
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
Solves for the y and z coordinates and velocities of the ball when it comes in
contact with the angled backboard.

Parameters: phi - angle the backboard makes with the vertical z axis
            y_points, z_points, v_y_points, v_z_points - the y and z positions
            and velocities of the ball for its whole trajectory
            T - array of time indices

Returns: y_before, z_before, v_y_before, v_z_before - the y and z positions and
        velocities of the ball up until it hits the backboard
"""
def backboard_hit_location(phi,y_points,z_points,v_y_points,v_z_points,T):
    C = 0.7493 # circumference in m, from basketball's circumference of 29.5 inches
    radius = C/(2*np.pi)
    t = T[1]-T[0]
    length = len(y_points) #this could have been any of the four arrays since they're the same length
    height_backboard = 1.0668
    hits_backboard = False

    index = -1
    for i in range(length):
        if y_points[i]-radius >= 0 and y_points[i]+radius <= height_backboard*np.sin(phi):
            if z_points[i]-radius >= 3.048 and z_points[i]+radius <= 3.048+height_backboard+np.cos(phi):
                index = i
                hits_backboard = True

    if hits_backboard:
        y_before = y_points[:index]
        z_before = z_points[:index]
        v_y_before = v_y_points[:index]
        v_z_before = v_z_points[:index]

        v = np.sqrt((v_y_before[-1])**2+(v_z_before[-1])**2)
        theta = np.arctan(v_z_before[-1]/v_y_before[-1])

        backboard_hit = [y_before[-1],z_before[-1],v*np.cos(np.pi-theta-2*phi),v*np.sin(np.pi-theta-2*phi)]
        y_after,z_after,v_y_after,v_z_after = RK4(backboard_hit)

        y = y_before + y_after
        z = z_before + z_after
        v_y = v_y_before + v_y_after
        v_z = v_z_before + v_z_after

    else:
        y = y_points
        z = z_points
        v_y = v_y_points
        v_z = v_z_points


    return hits_backboard,y,z,v_y,v_z


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
        if np.absolute(z_points[i]-3.048) <= eps:
            if ((y_points[i]+radius) <= 0.6 and
                (y_points[i]-radius) >= 0.15):
                shot_made = True
                return shot_made

    return False


"""
Runs a Monte Carlo Simulation to determine the optimal angle of the backboard
for which the most random shots go in the basket.

Parameters:

Returns:
"""
def monte_carlo():
    """
    Constants for linspace
    """
    a = 0
    b = 1000
    N = 1000000
    T = np.linspace(a,b,N)
    h = (b-a)/N

    number_of_experiments = 100
    backboard_angles = np.arange(-90,90,0.1)
    start_ys = np.arange(0,10,)
    start_zs = np.arange(1,2,0.1)
    v0s = np.arange(4,9,0.1)
    start_angles = np.arange(40,85,0.1)
    start_thetas = []
    for angle in start_angles:
        start_thetas.append((np.pi/180) * angle)
    percent_in = []

    best_angle = 0
    most_in = -1
    y = []
    z = []
    v_y = []
    v_z = []

    for phi in backboard_angles:
        shots_hit_backboard = number_of_experiments
        shots_made = 0
        for n in range(number_of_experiments):
            random_y = np.random.choice(start_ys)
            random_z = np.random.choice(start_zs)
            random_v0 = np.random.choice(v0s)
            random_theta = np.random.choice(start_thetas)

            init = [random_y, random_z, random_v0*np.cos(random_theta), random_v0*np.sin(random_theta)]
            y_points, z_points, v_y_points, v_z_points = RK4(init)
            hits_backboard, y_backboard, z_backboard, v_y_backboard, v_z_backboard = backboard_hit_location(phi,y_points,z_points,v_y_points,v_z_points,T)

            if hits_backboard:
                shot_in_basket = in_basket(y_backboard,z_backboard)
                if shot_in_basket:
                    shots_made += 1
            else:
                shots_hit_backboard -= 1
        try:
            percent_in.append(shots_made/shots_hit_backboard)
            if shots_made/shots_hit_backboard >= most_in:
                best_angle = phi
                most_in = shots_made/shots_hit_backboard
                y = random_y
                z = random_z
                v_y = random_v_y
                v_z = random_v_z
        except:
            percent_in.append(-1)
        
    return best_angle, most_in, y, z, v_y, v_z

def main():
    height_backboard = 1.0668

    phi,num_in,y_points,z_points,v_y_points,v_z_points =  monte_carlo()
    plt.plot(y_points,z_points)
    plt.plot(np.linspace(0,height_backboard*np.sin(phi),1000),np.linspace(3.048,3.048+height_backboard*np.cos(phi),1000))
    plt.show()

if __name__ == "__main__":
    main()
