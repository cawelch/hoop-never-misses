"""
file: monte_carlo.py

edited for three point backboard

Author: Caitlin Welch
Date created: September 20, 2020
Date modified: October 1, 2020

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
def hit_backboard(phi,y,z,fixed_y,fixed_z,next_z): #WILL HAVE TO ADJUST THIS FUNCTION TO ACCOUNT FOR SEVERAL BACKBOARD SEGMENTS
    C = 0.7493 # circumference in m, from basketball's circumference of 29.5 inches
    radius = C/(2*np.pi)
    height_backboard = 1.0668
    eps_backboard = 1e-2
    eps = 1e-1
    slope = np.tan(phi)

    theta = np.linspace(0,2*np.pi,100)
    circle_ys = y+radius*np.cos(theta)
    circle_zs = z+radius*np.sin(theta)
    #plt.plot(circle_ys,circle_zs)

    for i in range(len(circle_ys)):
        plt.plot(circle_ys[i],circle_zs[i],'.')
        if phi < np.pi/2:
            if np.absolute(circle_ys[i]) >= eps_backboard and circle_ys[i] <= (next_z-fixed_z)/slope:
                if circle_zs[i]-radius >= fixed_z and circle_zs[i] <= next_z:
                    #print('within rectangle')
                    #print(np.absolute(circle_zs[i]-3.048-slope*circle_ys[i]))
                    if np.absolute(circle_zs[i]-fixed_z-slope*circle_ys[i]) <= eps:
                        #print('hit')
                        return True
        elif phi == np.pi/2:
            if np.absolute(circle_ys[i]) <= eps_backboard:
                if circle_zs[i]-radius >= fixed_z and circle_zs[i] <= next_z:
                    #print('hit')
                    return True
        else:
            if circle_ys[i] >= (next_z-fixed_z)/slope and circle_ys[i] <= eps_backboard:
                if circle_zs[i]-radius >= fixed_z and circle_zs[i] <= next_z:
                    #print('within rectangle')
                    if np.absolute(circle_zs[i]-fixed_z-slope*circle_ys[i]) <= eps:
                        #print('hit')
                        return True

    return False

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
            phi - angle that the backboard makes with the horizontal
Returns: solved differential equations of the y and z positions and velocities
"""
def RK4(init,phi,fixed_y,fixed_z,next_z):
    y_points = []
    z_points = []
    v_y_points = []
    v_z_points = []
    backboard = False
    height_backboard = 1.0668

    r = np.array(init, float)
    h = 0.01
    while (r[1] >= 3.048 or r[3] > 0) and r[0] >= (next_z-fixed_z)/np.tan(phi)-0.1:
        #print("y=",r[0])
        y_points.append(r[0])
        z_points.append(r[1])
        v_y_points.append(r[2])
        v_z_points.append(r[3])

        if (not backboard) and hit_backboard(phi,r[0],r[1],fixed_y,fixed_z,next_z):
            backboard = True
            y_bounce, z_bounce, v_y_bounce, v_z_bounce = elastic_bounce(phi,r[0],r[1],r[2],r[3])
            r = np.array([y_bounce,z_bounce,v_y_bounce,v_z_bounce],float)
        else:
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
shots go in the backboard.

Parameters: phi_array - array of angles that the backboard makes with the horizontal

Returns: float value of the percent of shots that go in the hoop
"""
def percent_in(phi_array,fixed_y,fixed_z,next_z):
    num_shots = 1
    length = len(phi_array)
    num_backboard = np.array([0]*length)
    shots_made = np.array([0]*length)
    pct = np.array([0]*length)
    smallest_y = 4.1148/(np.tan(np.arctan(1.0668/.6)))

    while np.amin(num_backboard) <= 1000:
        random_y =  np.random.uniform(smallest_y,10)
        random_z =  np.random.uniform(1.5,2)
        random_v0 = np.random.uniform(4,9)
        random_theta = np.random.uniform(0,np.pi)

        init = [random_y,random_z,random_v0*np.cos(random_theta),random_v0*np.sin(random_theta)]
<<<<<<< HEAD
        y_points,z_points,v_y_points,v_z_points = RK4(init,phi)
            #plt.plot(y_points,z_points)

        if y_points != []:
            num_backboard += 1

        for k in range(len(y_points)):
            if in_basket(y_points[k],z_points[k]):
                shots_made += 1
                    #plt.plot(y_points,z_points)
                break

    print(num_shots,num_backboard,shots_made)


=======
        for j in range(length):
            y_points,z_points,v_y_points,v_z_points = RK4(init,phi_array[j],fixed_y,fixed_z,next_z)
            plt.plot(y_points,z_points)
            #print(num_backboard,shots_made)
            if y_points != []:
                num_backboard[j] += 1

            for k in range(len(y_points)):
                if in_basket(y_points[k],z_points[k]):
                    shots_made[j] += 1
                        #plt.plot(y_points,z_points)
                    break

    print(num_shots,num_backboard,shots_made)
>>>>>>> three_points
    try:
        #percent is number of shots made for shots that hit the backboard. If
        #no shots hit the backboard, set the percent to 0.
        pct = shots_made/num_backboard
        #print("HERE",pct)
    except:
        for i in range(length):
            if num_backboard[i] == 0:
                pct[i] = 0
            else:
                pct[i] = np.float(shots_made[i])/np.float(num_backboard[i])
            #print(pct)
    #print(pct,type(pct))
    return pct


"""
Finds the best angle for the backboard, from a given array of backboard angles.

Parameters: phi_array - array of backboard angles

Returns: best_phi - idea backboard angle from the array
        best_pct - percentage of random shots that go in for the best_phi
                    backboard angle (note best_pct is the greatest of all percentages)
"""
def best_angle(min_phi,max_phi,fixed_y,fixed_z,next_z):
    height_backboard = 1.0668
    best_pct = 0.0
<<<<<<< HEAD
    num_angles = 100
    pct = []
    phi_array = []

    for i in range(num_angles):
        phi_array.append(np.random.uniform(min_phi,max_phi))
 #plt.plot(np.linspace(0,height_backboard*np.cos(phi),1000),np.linspace(3.048,3.048+height_backboard*np.sin(phi),1000))
        pct.append(percent_in(phi_array[i]))
        #print(pct)
        #print("Percent in: ",pct,"Angle: ",phi)
    
    #    if pct >= best_pct:
     #       best_pct = pct
      #      best_phi = phi
    
=======
    num_angles = 10
    phi_array = []

    for i in range(num_angles):
        phi = np.random.uniform(min_phi,max_phi)
        phi_array.append(phi)
        """
        pct = percent_in(phi,fixed_y,fixed_z,next_z)
        print("Percent in: ",pct,"Angle: ",phi)
        if pct >= best_pct:
            best_pct = pct
            best_phi = phi
        #plt.plot(np.linspace(0,height_backboard*np.cos(phi),1000),np.linspace(3.048,3.048+height_backboard*np.sin(phi),1000))
        """

    pct = percent_in(phi_array,fixed_y,fixed_z,next_z)
    #print("1",pct,type(pct))
    #print(pct)
>>>>>>> three_points
    best_index = np.argmax(pct)
    best_pct = pct[best_index]
    best_phi = phi_array[best_index]
    '''
    if 0 in num_backboard:
        pct = []
        for num in num_backboard:
            if num == 0:
                pct.append(0)
            else:
                pct.append(np.float(shots_made)/np.float(num_backboard))
    else:
        pct = np.float(shots_made)/np.float(num_backboard)
    '''


<<<<<<< HEAD
=======


>>>>>>> three_points
    return best_phi


"""
Optimization function, which narrows the range of angles tested but with the same
number of points, to try to find convergence for the backboard angle with the
highest percent of shots made.

Returns: central_phi - the angle at which the highest percentage of shots in
        converges to
"""
def optimize_angle(fixed_y,fixed_z,next_z):
    height_backboard = 1.0668
    central_phi = [np.pi/2]
    n = 2 #we start by testing points plus/minus pi/2 of the central_phi
    tolerance = 1e-2
    diff = np.pi/2 #start at value of pi/2 because we know the while loop condition will be true
    index = 0

    while diff >= tolerance:
        print("Tolerance: ",index)
        min_phi = central_phi[index]-np.pi/n
        max_phi = central_phi[index]+np.pi/n
        central_phi.append(best_angle(min_phi,max_phi,fixed_y,fixed_z,next_z))
        plt.plot(np.linspace(0,height_backboard*np.cos(central_phi),1000),np.linspace(3.048,3.048+height_backboard*np.sin(central_phi),1000))
        index += 1
        n *= 2
        diff = np.absolute(central_phi[index]-central_phi[index-1])
        print(np.absolute(diff-tolerance))

    print(central_phi)
    return central_phi[-1]

"""
def move_points():
    num_points = 3
    y = [0]
    z = [3.048]

    for j in range(1,num_points):
        z.append(z[j-1]+(4.1148-3.048)/(num_points-1))
        phi = optimize_angle(y[j-1],z[j-1],z[j])
        next_y = z[j]/np.tan(phi)
        if next_y >= 2 or next_y <= -2:
            phi = optimize_angle(y[j-1],z[j-1],z[j])
            next_y = z[j]/np.tan(phi)
        y.append(next_y)
        plt.plot(np.linspace(y[j-1],y[j],100),np.linspace(z[j-1],z[j],100))

    return y,z
"""


def main():
    height_backboard = 1.0668
    C = 0.7493 # circumference in m, from basketball's circumference of 29.5 inches
    radius = C/(2*np.pi)

    #print(move_points())

    print("Optimized angle is: ",optimize_angle(0,3.048,4.1148))

    """
    phi_array = np.array([np.pi/4,np.pi/3,np.pi/2,2*np.pi/3,3*np.pi/4])
    best_phi, best_pct = best_angle(phi_array)
    #print("Best angle: ",best_phi,"Percent of shots made: ",best_pct)
    """
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
