"""
file: monte_carlo.py

Author: Caitlin Welch
Date created: September 20, 2020
Date modified: September 24, 2020

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
    slope = np.tan(phi)

    """
    z_range = np.linspace(3.048,3.048+height_backboard*np.sin(phi),1000)
    y_range = (z_range-3.048)/np.tan(phi)
    """

    theta = np.linspace(0,2*np.pi,500)
    circle_ys = y+radius*np.cos(theta)
    circle_zs = z+radius*np.sin(theta)
    #plt.plot(circle_ys,circle_zs)

    for i in range(len(circle_ys)):
        #plt.plot(circle_ys[i],circle_zs[i],'.')
        if phi <= np.pi/2:
            if np.absolute(circle_ys[i]) >= eps_backboard and circle_ys[i] <= height_backboard*np.cos(phi):
                if circle_zs[i] >= 3.048 and circle_zs[i] <= 3.048+height_backboard*np.sin(phi):
                    #print('within rectangle')
                    if np.absolute(circle_zs[i]-3.048-slope*circle_ys[i]) <= eps_backboard:
                        #print('hit')
                        return True
        else:
            if circle_ys[i] >= height_backboard*np.cos(phi) and circle_ys[i] <= eps_backboard:
                if circle_zs[i] >= 3.048 and circle_zs[i] <= 3.048+height_backboard*np.sin(phi):
                    #print('within rectangle')
                    if np.absolute(circle_zs[i]-3.048-slope*circle_ys[i]) <= eps_backboard:
                        #print('hit')
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

    if phi < np.pi/2:
        #print(v_z)
        if v_z > 0:
            alpha = np.pi/2-2*phi-theta
            #alpha = np.pi/2-2*phi+theta
            #print(alpha)
            v_y_after = v*np.cos(alpha)
            v_z_after = v*np.sin(alpha)
        else:
            alpha = np.pi/2-2*phi+theta
            v_y_after = v*np.cos(alpha)
            v_z_after = v*np.sin(alpha)
    elif phi == np.pi/2:
        v_y_after = -v_y
        v_z_after = v_z
    else:
        if v_z > 0:
            alpha = -np.pi/2+2*phi+theta
            v_z_after = v*np.sin(alpha)
            v_y_after = v*np.cos(alpha)
        else:
            alpha = np.pi-2*phi+theta
            v_y_after = v*np.cos(alpha)
            v_z_after = -v*np.sin(alpha)


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
    height_backboard = 1.0668

    r = np.array(init, float)
    h = 0.01
    while (r[1] >= 3.048 or r[3] > 0) and r[0] >= height_backboard*np.cos(phi)-0.1:
        #print("y=",r[0])
        y_points.append(r[0])
        z_points.append(r[1])
        v_y_points.append(r[2])
        v_z_points.append(r[3])

        if (not backboard) and hit_backboard(phi,r[0],r[1]):
            backboard = True
            y_bounce, z_bounce, v_y_bounce, v_z_bounce = elastic_bounce(phi,r[0],r[1],r[2],r[3])
            r = np.array([y_bounce,z_bounce,v_y_bounce,v_z_bounce],float)
        else:
            k1 = h*f(r)
            k2 = h*f(r+0.5*k1)
            k3 = h*f(r+0.5*k2)
            k4 = h*f(r+k3)
            r += (k1+2*k2+2*k3+k4)/6


    #if backboard:
    return y_points, z_points, v_y_points, v_z_points
    #else:
    #    return [],[],[],[]


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
shots go in the backboard.

Parameters: phi - angle that the backboard makes with the horizontal

Returns: float value of the percent of shots that go in the hoop
"""
def percent_in(phi):
    num_shots = 10
    shots_made = 0
    smallest_y = 4.1148/(np.tan(np.arctan(1.0668/.6)))
    start_ys =  np.array([2.5]) #np.arange(smallest_y,10,0.1)
    start_zs =  np.array([1.8]) #np.arange(1,2,0.1)
    v0s =  np.array([7.5]) #np.arange(4,9,0.1)
    start_angles =  np.array([125]) #np.arange(105,140,0.1)
    start_thetas = np.pi/180*start_angles

    for i in range(num_shots):
        random_y =  np.random.choice(start_ys)
        random_z =  np.random.choice(start_zs)
        random_v0 = np.random.choice(v0s)
        random_theta = np.random.choice(start_thetas)

        init = [random_y,random_z,random_v0*np.cos(random_theta),random_v0*np.sin(random_theta)]
        y_points,z_points,v_y_points,v_z_points = RK4(init,phi)
        plt.plot(y_points,z_points)

        for j in range(len(y_points)):
            if in_basket(y_points[j],z_points[j]):
                shots_made += 1
                break
                #plt.plot(y_points,z_points)

    pct = np.float(shots_made)/np.float(num_shots)
    return pct


"""
Finds the best angle for the backboard, from a given array of backboard angles.

Parameters: phi_array - array of backboard angles

Returns: best_phi - idea backboard angle from the array
        best_pct - percentage of random shots that go in for the best_phi
                    backboard angle (note best_pct is the greatest of all percentages)
"""
def best_angle(phi_array):
    height_backboard = 1.0668
    best_pct = 0.0
    best_phi = np.pi/2 #assume the best angle is the angle that the backboard is in regulation and change this if there's a better angle

    for phi in phi_array:
        plt.plot(np.linspace(0,height_backboard*np.cos(phi),1000),np.linspace(3.048,3.048+height_backboard*np.sin(phi),1000))
        if percent_in(phi) > best_pct:
            best_pct = percent_in(phi)
            best_phi = phi

    return best_phi, best_pct


def main():
    height_backboard = 1.0668
    C = 0.7493 # circumference in m, from basketball's circumference of 29.5 inches
    radius = C/(2*np.pi)
    phi_array = np.array([np.pi/2,np.pi/3])
    best_phi, best_pct = best_angle(phi_array)
    print(best_phi,best_pct)

    """
    Plotting for the backboard and hoop.
    """
    #plt.plot(np.linspace(0,height_backboard*np.cos(best_phi),1000),np.linspace(3.048,3.048+height_backboard*np.sin(best_phi),1000))
    plt.plot(np.linspace(0.6,0.15,100),[3.048]*100)
    plt.plot([0.15+radius]*100,np.linspace(3.148,2.948,100))
    plt.plot([0.6-radius]*100,np.linspace(3.148,2.948,100))
    plt.axis('scaled')
    plt.show()


if __name__ == "__main__":
    main()
