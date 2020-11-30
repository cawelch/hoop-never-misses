"""
File: monte_carlo_two_d.py
Author: Caitlin Welch
Date modified: November 30, 2020

Brief: Finds the optimal backboard for a given number of points on the backboard,
        a given number of random backboards to choose from and a given number
        of initial shots to simulate. The optimal backboard is the backboard
        that has the greatest percent of shots that go in the basket for the
        random shots.
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

Parameters: phi_array - array of angles of the backboard with the horizontal
            y, z - current y and z position of the ball
            backboard_y, backboard_z - arrays of the y and z coordinates of the
                                        random backboard points

Returns: boolean value indicating whether the ball is hitting the backboard
"""
def hit_backboard(phi_array,y,z,backboard_y,backboard_z):
    C = 0.7493 # circumference in m, from basketball's circumference of 29.5 inches
    radius = C/(2*np.pi)
    height_backboard = 1.0668
    eps_backboard = 1e-4
    eps = 1e-1
    #print(phi_array)

    theta = np.linspace(0,2*np.pi,100)
    circle_ys = y+radius*np.cos(theta)
    circle_zs = z+radius*np.sin(theta)
    for j in range(len(phi_array)):
        phi = phi_array[j]
        slope = np.tan(phi_array[j])
        for i in range(len(circle_ys)):
            on_slope = np.absolute(circle_zs[i]-backboard_z[j]-slope*(circle_ys[i]-backboard_y[j]))
            if phi <= np.pi/2:
                if circle_ys[i] >= backboard_y[j] and circle_ys[i] <= backboard_y[j+1] and circle_zs[i] >= backboard_z[j] and circle_zs[i] <= backboard_z[j+1] and on_slope <= eps:
                    hit_phi = phi
                    return True, hit_phi
            else:
                if circle_ys[i] <= backboard_y[j] and circle_ys[i] >= backboard_y[j+1] and circle_zs[i] >= backboard_z[j] and circle_zs[i] <= backboard_z[j+1] and on_slope <= eps:
                    hit_phi = phi
                    return True, hit_phi

    hit_phi = -1
    return False,hit_phi


"""
For a shot that hits the backboard, determines the y and z positions and velocities
of the shot immediately after its collision with the backboard.

Parameters: phi - angle of the backboard with the horizontal
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
        if v_z > 0:
            alpha = np.pi/2-2*phi-theta
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
            v_y_after = -v*np.cos(alpha)
            v_z_after = v*np.sin(alpha)
        else:
            alpha = np.pi-2*phi+theta
            v_y_after = v*np.cos(alpha)
            v_z_after = -v*np.sin(alpha)


    return y_after,z_after,v_y_after,v_z_after


"""
Function used to solve the differential equations, using the Fourth order
Runge-Kutta method.

Parameters: init - array of the initial y and z positions and velocities
            phi_array - array of angles that the backboard makes with the horizontal
            backboard_y, backboard_z - arrays of the y and z coordinates of the
                                        random backboard points
Returns: solved differential equations of the y and z positions and velocities
"""
def RK4(init,phi_array,backboard_y,backboard_z):
    y_points = []
    z_points = []
    v_y_points = []
    v_z_points = []
    backboard = False
    stop = False
    basket = False
    stop_index = 0
    height_backboard = 1.0668
    min_y = -0.5

    r = np.array(init, float)
    h = 0.01
    while (r[1] >= 3.048 or r[3] > 0) and r[0] >= min_y:
        y_points.append(r[0])
        z_points.append(r[1])
        v_y_points.append(r[2])
        v_z_points.append(r[3])

        backboard, hit_phi = hit_backboard(phi_array,r[0],r[1],backboard_y,backboard_z)
        if (not stop) and backboard:
            stop = True
            stop_index = len(y_points)
            y_bounce, z_bounce, v_y_bounce, v_z_bounce = elastic_bounce(hit_phi,r[0],r[1],r[2],r[3])
            r = np.array([y_bounce,z_bounce,v_y_bounce,v_z_bounce],float)

        else:
            k1 = h*f(r)
            k2 = h*f(r+0.5*k1)
            k3 = h*f(r+0.5*k2)
            k4 = h*f(r+k3)
            r += (k1+2*k2+2*k3+k4)/6

        if in_basket(r[0],r[1]):
            basket = True

    if stop:
        plt.plot(y_points[:stop_index],z_points[:stop_index],'b')
        plt.plot(y_points[stop_index:],z_points[stop_index:],'r')
        return y_points, z_points, v_y_points, v_z_points,basket
    else:
        return [],[],[],[],basket


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
Chooses num_points random coordinates for the backboard. The y coordinates are
random within some given range, and the z coordinates are evenly spaced along the
range of 3.058 to 4.1148m. Creates a random shape backboard.

Parameters: num_points - the number of points that the backboard will be composed of

Returns: backboard_y, backboard_z - the y and z coordinates of the random backboard
"""
def move_points(num_points):
    backboard_y = np.zeros(num_points) #[0,-0.1,0]
    backboard_z = np.zeros(num_points) #[3.048,3.5814,4.1148]
    backboard_z[0] = 3.048

    for i in range(1,num_points):
        backboard_z[i]=(backboard_z[i-1]+(4.1148-3.048)/(num_points-1))
        backboard_y[i]=(np.random.uniform(-0.5,0.5))
    plt.plot(backboard_y,backboard_z)
    return backboard_y,backboard_z


"""
Simulates num_shots shots at num_backboards random backboards composed of
num_backboard_points. Keeps a tally of the number of shots that go in the basket
for each backboard.

Parameters - num_backboard_points - integer number of points that make up the backboard
             num_backboards - integer number of random backboards for which we will compare the percent of shots in
             num_shots - integer number of random shots thrown at each backboard

Returns - percent_shots_in - array of the percent of shots that go in the basket
                             for each random backboard
          backboard_y, backboard_z - nested array of the y and z coordinates for
                                     each random backboard
"""
def percent_in(num_backboard_points, num_backboards, num_shots):
    backboard_y = np.zeros((num_backboards,num_backboard_points))
    backboard_z = np.zeros((num_backboards,num_backboard_points))
    shots_in = np.zeros(num_backboards)

    for i in range(num_backboards):
        back_y, back_z = move_points(num_backboard_points)
        backboard_y[i] = back_y
        backboard_z[i] = back_z
        del_y = []
        del_z = []
        for m in range(len(backboard_y[i])-1):
            del_y.append(backboard_y[i][m+1]-backboard_y[i][m])
            del_z.append(backboard_z[i][m+1]-backboard_z[i][m])

        if 0 in del_y:
            for k in range(len(del_y)):
                phi_array = []
                if del_y[k] == 0:
                    phi_array.append(np.pi/2)
                else:
                    phi_array.append(np.arctan(np.array(del_z[k])/np.array(del_y[k])))
        else:
            phi_array = np.arctan(np.array(del_z)/np.array(del_y))

        for j in range(len(phi_array)):
            if phi_array[j] < 0:
                phi_array[j] += np.pi
    index = 0
    for k in range(num_shots):
        random_y = 7.5 #np.random.uniform(0.5,7)
        random_z = 1.8 #np.random.uniform(1.5,2)
        random_v0 = 9.8 #np.random.uniform(6,10)
        random_theta = np.random.uniform(np.pi/2,np.pi)
        init = [random_y,random_z,random_v0*np.cos(random_theta),random_v0*np.sin(random_theta)]

        for i in range(num_backboards):
            y_points,z_points,v_y_points,v_z_points,basket = RK4(init,phi_array,backboard_y[i],backboard_z[i])

            if basket:
                shots_in[i] += 1

    percent_shots_in = np.array(shots_in)/np.float(num_shots)
    return percent_shots_in,backboard_y,backboard_z


"""
Finds the best backboard out of a set of randomly generated backboards.

Parameters - num_backboard_points - integer number of points that make up the backboard
             num_backboards - integer number of random backboards for which we will compare the percent of shots in
             num_shots - integer number of random shots thrown at each backboard

Returns: max_percent - percent of shots that go in the basket for the optimal
                       backboard of the randomly generated backboards
         backboard_y, backboard_z - y and z coordinates of the optimal backboard
                                    of the randomly generated backboards
"""
def best_backboard(num_backboard_points, num_backboards, num_shots):
    percent_shots_in,backboard_y_all,backboard_z_all = percent_in(num_backboard_points, num_backboards, num_shots)
    max_index = np.argmax(percent_shots_in)
    max_percent = percent_shots_in[max_index]
    backboard_y = backboard_y_all[max_index]
    backboard_z = backboard_z_all[max_index]

    return max_percent,backboard_y,backboard_z


def main():
    num_backboard_points = 10
    num_backboards = 1000
    num_shots = 1000
    print(best_backboard(num_backboard_points, num_backboards, num_shots))
    plt.plot(np.linspace(0.6,0.15,100),[3.048]*100)
    plt.savefig("plot-output.png")
    plt.show()

if __name__ == "__main__":
    main()
