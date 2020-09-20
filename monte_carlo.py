"""
file: monte_carlo.py

Author: Caitlin Welch
Date created: September 20, 2020
Date modified: September 20, 2020

Brief: Attempt at running a Monte Carlo simulation for a 2D shot to form a
        continuous curved bacboard.

Current progress: solving for percent of shots in for one backboard angle
"""

import numpy as np
import pylab as plt

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
For given y and z coordinates, determines whether the ball is currently hitting the backboard.

Parameters: phi - angle of the backboard with the horizontal
            y, z - current y and z position of the ball

Returns: boolean value indicating whether the ball is hitting the backboard
"""
def hit_backboard(phi,y,z):
    C = 0.7493 # circumference in m, from basketball's circumference of 29.5 inches
    radius = C/(2*np.pi)
    height_backboard = 1.0668
    eps_backboard = 1e-2

    theta = np.linspace(0,2*np.pi,500)
    circle_ys = radius*np.cos(theta)
    circle_zs = radius*np.sin(theta)

    for i in range(len(circle_ys)):
        if phi <= np.pi/2:
            if np.absolute(y-circle_ys[i]) >= eps_backboard and y+circle_ys[i] <= height_backboard*np.cos(phi):
                if z-circle_zs[i] >= 3.048 and z+circle_zs[i] <= 3.048+height_backboard*np.sin(phi):
                    return True
        else:
            if y-circle_ys[i] >= height_backboard*np.cos(phi) and y+circle_ys[i] <= eps_backboard:
                if z-circle_zs[i] >= 3.048 and z+z_points[i] <= 3.048+height_backboard*np.sin(phi):
                    return True

        return False

"""
For a shot that hits the backboard, determines the y and z positions and velocities
of the shot immediately after its collision with the backboard.

Parameters: phi - angle of teh backboard with the horizontal
            y, z, v_y, v_z - y and z positions and velocities of the ball upon
            hitting the backboard

Returns: y and z positions and velocities of the ball immediately after its
        collision with the backboard
"""
def elastic_bounce(phi,y,z,v_y,v_z):
    v = np.sqrt((v_y)**2+(v_z)**2)
    theta = np.arctan(v_z/v_y)

    y_after = y
    z_after = z

    if phi >= np.pi/2:
        v_y_after = -v*np.cos(2*phi-theta-np.pi/2)
        if v_z > 0:
            v_z_after = -v*np.sin(2*phi-theta-np.pi/2)
        else:
            v_z_after = v*np.sin(2*phi-theta-np.pi/2)
    else:
        v_y_after = -v*np.cos(np.pi/2+theta-2*phi)
        if v_z > 0:
            v_z_after = -v*np.cos(np.pi/2+theta-2*phi)
        else:
            v_z_after = v*np.cos(np.pi/2+theta-2*phi)

    return y_after,z_after,v_y_after,v_z_after

"""
Function used to solve the differential equations, using the Fourth order
Runge-Kutta method.

Parameters: init - array of the initial y and z positions and velocities
            phi - angle that the backboard makes with the horizontal
Returns: solved differential equations of the y and z positions and velocities
"""
def RK4(init,phi):
    y_points = []
    z_points = []
    v_y_points = []
    v_z_points = []
    backboard = False

    r = np.array(init, float)
    h = 0.01
    while r[1] >= 3.048 or r[3] > 0:
        if hit_backboard(phi,r[0],r[1]):
            backboard = True
            y_bounce, z_bounce, v_y_bounce, v_z_bounce = elastic_bounce(phi,r[0],r[1],r[2],r[3])
            break
        else:
            y_points.append(r[0])
            z_points.append(r[1])
            v_y_points.append(r[2])
            v_z_points.append(r[3])

            k1 = h*f(r)
            k2 = h*f(r+0.5*k1)
            k3 = h*f(r+0.5*k2)
            k4 = h*f(r+k3)
            r += (k1+2*k2+2*k3+k4)/6

    if backboard:
        r = np.array([y_bounce,z_bounce,v_y_bounce,v_z_bounce],float)
        while r[1] >= 3.048 or r[3] > 0:
            y_points.append(r[0])
            z_points.append(r[1])
            v_y_points.append(r[2])
            v_z_points.append(r[3])

            k1 = h*f(r)
            k2 = h*f(r+0.5*k1)
            k3 = h*f(r+0.5*k2)
            k4 = h*f(r+k3)
            r += (k1+2*k2+2*k3+k4)/6


    if backboard:
        return y_points, z_points, v_y_points, v_z_points
    else:
        return [],[],[],[]


"""
Determines whether, given the y and z position of the ball at a given time, the
ball is in the basket.

Parameters: y,z - y and z position of the ball

Returns: boolean value indicating whether or not the ball is in the basket
"""
def in_basket(y,z):
    C = 0.7493  # circumference in m, from basketball's circumference of 29.5 inches
    radius = C/(2*np.pi)
    eps_rim = 1e-1

    if np.absolute(z-3.048) <= eps_rim:
        if y+radius <= 0.6 and y-radius >= 0.15:
            return True

    return False


"""
For a given angle of the backboard, determines how many randomly generated
shots go in teh backboard.

Parameters: phi - angle that the backboard makes with the horizontal

Returns: float value of the percent of shots that go in the hoop
"""
def percent_in(phi):
    num_shots = 10000
    shots_made = 0
    start_ys = np.arange(0,10,0.1)
    start_zs = np.arange(1,2,0.1)
    v0s = np.arange(4,9,0.1)
    start_angles = np.arange(105,140,0.1)
    start_thetas = np.pi/180*start_angles

    for i in range(num_shots):
        random_y = np.random.choice(start_ys)
        random_z = np.random.choice(start_zs)
        random_v0 = np.random.choice(v0s)
        random_theta = np.random.choice(start_thetas)

        init = [random_y,random_z,random_v0*np.cos(random_theta),random_v0*np.sin(random_theta)]
        y_points,z_points,v_y_points,v_z_points = RK4(init,phi)

        for j in range(len(y_points)):
            if in_basket(y_points[j],z_points[j]):
                shots_made += 1
                plt.plot(y_points,z_points)

    pct = np.float(shots_made)/np.float(num_shots)
    return pct


def main():
    height_backboard = 1.0668
    C = 0.7493 # circumference in m, from basketball's circumference of 29.5 inches
    radius = C/(2*np.pi)
    phi = np.pi/2
    print(percent_in(phi))
    if phi <= np.pi/2 and phi >= 0:
        plt.plot(np.linspace(0,height_backboard*np.cos(phi),1000),np.linspace(3.048,3.048+height_backboard*np.sin(phi),1000))
    elif phi >= np.pi/2 and phi <= np.pi:
        plt.plot(np.linspace(0,height_backboard*np.cos(phi),1000),np.linspace(3.048,3.048+height_backboard*np.sin(phi),1000))
    plt.plot(np.linspace(0.6,0.15,100),[3.048]*100)
    plt.plot([0.15+radius]*100,np.linspace(3.148,2.948,100))
    plt.plot([0.6-radius]*100,np.linspace(3.148,2.948,100))
    plt.axis('scaled')
    plt.show()


if __name__ == "__main__":
    main()
