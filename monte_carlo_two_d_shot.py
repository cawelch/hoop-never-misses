"""
file: monte-carlo-2d-shot.py

Author: Caitlin Welch
Date created: September 6, 2020
Date modified: September 10, 2020

Brief: First attempt at running a Monte Carlo simulation for a 2D shot to form
        a continuous curved backboard
"""

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
def backboard_hit_location(phi,y_points,z_points,v_y_points,v_z_points):
    C = 0.7493 # circumference in m, from basketball's circumference of 29.5 inches
    radius = C/(2*np.pi)
    length = len(y_points) #this could have been any of the four arrays since they're the same length
    height_backboard = 1.0668
    hits_backboard = False

    index = -1
    for i in range(length):
        if phi <= np.pi/2 and phi >= 0:
            if y_points[i]-radius >= 0 and y_points[i]+radius <= height_backboard*np.cos(phi):
                if z_points[i]-radius >= 3.048 and z_points[i]+radius <= 3.048+height_backboard+np.sin(phi):
                    index = i
                    hits_backboard = True
        elif phi > np.pi/2 and phi <= np.pi:
            if y_points[i]-radius >= -height_backboard*np.cos(np.pi-phi) and y_points[i]+radius <= 0:
                if z_points[i]-radius >= 3.048 and z_points[i]+radius <= 3.048+height_backboard*np.sin(np.pi-phi):
                    index = i
                    hits_backboard = True

    #Checks whether the ball hits the wrong side of the backboard
    if phi <= np.pi/2 and hits_backboard:
        print("1")
        #If the ball is in the downward part of its trajectory and coming down at a more vertical angle than phi
        if v_z_points[index] < 0 and np.arctan(z_points[index]/y_points[index]) >= phi:
            hits_backboard = False
            print("2")
    elif phi >= np.pi/2 and hits_backboard:
        #If the ball is in the upward part of its trajectory and coming up at an angle
        if v_z_points[index] > 0 and np.arctan(z_points[index]/y_points[index]) <= phi:
            hits_backboard = False

    if hits_backboard:
        y_before = y_points[:index]
        z_before = z_points[:index]
        v_y_before = v_y_points[:index]
        v_z_before = v_z_points[:index]

        v = np.sqrt((v_y_before[-1])**2+(v_z_before[-1])**2)
        theta = np.arctan(v_z_before[-1]/v_y_before[-1])


        if phi > np.pi/2 and phi <= np.pi:
            if v_z_before[-1] > 0:
                backboard_hit = [y_before[-1],z_before[-1],-v*np.cos(2*phi-theta-np.pi/2),-v*np.sin(2*phi-theta-np.pi/2)]
            else:
                backboard_hit = [y_before[-1],z_before[-1],-v*np.cos(2*phi-theta-np.pi/2),v*np.sin(2*phi-theta-np.pi/2)]
        elif phi <= np.pi/2 and phi >= 0:
            if v_z_before[-1] > 0:
                backboard_hit = [y_before[-1],z_before[-1],-v*np.sin(np.pi/2+theta-2*phi),-v*np.cos(np.pi/2+theta-2*phi)]
            else:
                backboard_hit = [y_before[-1],z_before[-1],-v*np.sin(np.pi/2+theta-2*phi),v*np.cos(np.pi/2+theta-2*phi)]
        y_after, z_after, v_y_after, v_z_after = RK4(backboard_hit)

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
    Different initial conditions for the shot to randomly choose from.
    """
    number_of_shots = 100
    backboard_angles = np.linspace(np.pi/6,5*np.pi/6,1000)
    start_ys = np.arange(0,10,0.1)
    start_zs = np.arange(1,2,0.1)
    v0s = np.arange(4,9,0.1)
    start_angles = np.arange(40,85,0.1)
    start_thetas = []
    for angle in start_angles:
        start_thetas.append((np.pi/180) * angle)
    percent_in = []
    y = []
    z = []
    v_y = []
    v_z = []

    """
    Iterate through the potential angles the backboard can be at.
    """
    for i in range(len(backboard_angles)):
        y.append([])
        z.append([])
        v_y.append([])
        v_z.append([])
        phi = backboard_angles[i]
        """
        Variables that will help in finding the optimum angle for a bunch of different shots.
        """
        shots_hit_backboard = number_of_shots
        shots_made = 0
        """
        Iterate through (eventually a large number of) random shots.
        """
        for n in range(number_of_shots):
            random_y = np.random.choice(start_ys)
            random_z = np.random.choice(start_zs)
            random_v0 = np.random.choice(v0s)
            random_theta = np.random.choice(start_thetas)

            """
            Run the RK4 function and backboard location for our random initial shot conditions.
            """
            init = [random_y, random_z, random_v0*np.cos(random_theta), random_v0*np.sin(random_theta)]
            y_points, z_points, v_y_points, v_z_points = RK4(init)
            hits_backboard, y_backboard, z_backboard, v_y_backboard, v_z_backboard = backboard_hit_location(phi,y_points,z_points,v_y_points,v_z_points)

            """
            We only care about when the ball hits the backboard.
            """
            if hits_backboard:
                """
                Determines whether or not the ball goes in the basket after the rebound off the backboard.
                """
                shot_in_basket = in_basket(y_backboard,z_backboard)
                if shot_in_basket:
                    shots_made += 1 #keep track of number of shots that work for this angle
                    """
                    Keep track of the index where the ball hits the backboard if it goes in the
                    hoop and the y and z positions and velocities of the ball for that shot.
                    """
                    y[i].append(y_points)
                    z[i].append(z_points)
                    v_y[i].append(v_y_points)
                    v_z[i].append(v_z_points)
            else:
                shots_hit_backboard -= 1
        try:
            percent_in.append(shots_made/shots_hit_backboard)
        except:
            percent_in.append(-1)
    highest_percent = np.max(percent_in)
    index = percent_in.index(highest_percent)
    best_phi = backboard_angles[index]
    print(y)

    return best_phi,highest_percent,y[index],z[index],v_y[index],v_z[index]

def main():
    height_backboard = 1.0668
    phi,highest_percent,y,z,v_y,v_z = monte_carlo()
    for i in range(len(y)):
        plt.plot(y[i],z[i])
    if phi <= np.pi/2 and phi >= 0:
        plt.plot(np.linspace(0,height_backboard*np.cos(phi),1000),np.linspace(3.048,3.048+height_backboard*np.sin(phi),1000))
    elif phi >= np.pi/2 and phi <= np.pi:
        plt.plot(np.linspace(0,-height_backboard*np.cos(np.pi-phi),1000),np.linspace(3.048,3.048+height_backboard*np.sin(np.pi-phi),1000))
    plt.axis('scaled')
    plt.show()
    print(phi,highest_percent)


if __name__ == "__main__":
    main()
