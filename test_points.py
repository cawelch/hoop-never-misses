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
def hit_backboard(phi_array,y,z):
    C = 0.7493 # circumference in m, from basketball's circumference of 29.5 inches
    radius = C/(2*np.pi)
    height_backboard = 1.0668
    eps_backboard = 1e-2
    eps = 1e-1

    theta = np.linspace(0,2*np.pi,100)
    circle_ys = y+radius*np.cos(theta)
    circle_zs = z+radius*np.sin(theta)
    #plt.plot(circle_ys,circle_zs)

    for phi in phi_array:
        slope = np.tan(phi)
        for i in range(len(circle_ys)):
            #plt.plot(circle_ys[i],circle_zs[i],'.')
            if phi < np.pi/2:
                if np.absolute(circle_ys[i]) >= eps_backboard and circle_ys[i] <= height_backboard*np.cos(phi):
                    if circle_zs[i]-radius >= 3.048 and circle_zs[i] <= 3.048+height_backboard*np.sin(phi):
                        #print('within rectangle')
                        #print(np.absolute(circle_zs[i]-3.048-slope*circle_ys[i]))
                        if np.absolute(circle_zs[i]-3.048-slope*circle_ys[i]) <= eps:
                            #print('hit')
                            hit_phi = phi
                            return True,hit_phi
            elif phi == np.pi/2:
                if np.absolute(circle_ys[i]) <= eps_backboard:
                    if circle_zs[i]-radius >= 3.048 and circle_zs[i] <= 3.048+height_backboard:
                        #print('hit')
                        hit_phi = phi
                        return True,hit_phi
            else:
                if circle_ys[i] >= height_backboard*np.cos(phi) and circle_ys[i] <= eps_backboard:
                    if circle_zs[i]-radius >= 3.048 and circle_zs[i] <= 3.048+height_backboard*np.sin(phi):
                        #print('within rectangle')
                        if np.absolute(circle_zs[i]-3.048-slope*circle_ys[i]) <= eps:
                            #print('hit')
                            hit_phi = phi
                            return True,hit_phi

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
def RK4(init,phi_array):
    y_points = []
    z_points = []
    v_y_points = []
    v_z_points = []
    backboard = False
    stop = False
    basket = False
    height_backboard = 1.0668
    min_y = -0.5 #height_backboard*np.cos(phi)-0.1

    r = np.array(init, float)
    h = 0.01
    #print(r[1],r[3],r[0],height_backboard*np.cos(phi)-0.1)
    while (r[1] >= 3.048 or r[3] > 0) and r[0] >= min_y:
        #print("y=",r[0])
        #print("here")
        y_points.append(r[0])
        z_points.append(r[1])
        v_y_points.append(r[2])
        v_z_points.append(r[3])

        backboard, hit_phi = hit_backboard(phi_array,r[0],r[1])
        if (not stop) and backboard:
            backboard = True
            y_bounce, z_bounce, v_y_bounce, v_z_bounce = elastic_bounce(hit_phi,r[0],r[1],r[2],r[3])
            r = np.array([y_bounce,z_bounce,v_y_bounce,v_z_bounce],float)
            #add in basket here
            if in_basket(r[0],r[1]):
                basket = True

        else:
            k1 = h*f(r)
            k2 = h*f(r+0.5*k1)
            k3 = h*f(r+0.5*k2)
            k4 = h*f(r+k3)
            r += (k1+2*k2+2*k3+k4)/6

    """
    if backboard:
        return y_points, z_points, v_y_points, v_z_points,basket
    else:
        return [],[],[],[],basket
    """

    return y_points,z_points,v_y_points,v_z_points,basket


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


def move_points(num_points):
    backboard_y = [0]
    backboard_z = [3.048]

    for i in range(1,num_points):
        backboard_z.append(backboard_z[i-1]+(4.1148-3.048)/(num_points-1))
        backboard_y.append(np.random.uniform(-0.5,0.5))

    plt.plot(backboard_y,backboard_z)
    return backboard_y,backboard_z



def main():
    backboard_y,backboard_z = move_points(3)
    plt.plot(backboard_y,backboard_z)

    random_y = 8    #np.random.uniform(smallest_y,7)
    random_z = 1.8 #np.random.uniform(1.5,2)
    random_v0 = 9.8 #np.random.uniform(6,10)
    random_theta = 2*np.pi/3 #np.random.uniform(np.pi/2,np.pi)
    init = [random_y,random_z,random_v0*np.cos(random_theta),random_v0*np.sin(random_theta)]
    phi_array = np.arctan(np.array(backboard_z),np.array(backboard_y))
    y_points,z_points,v_y_points,v_z_points,basket = RK4(init,phi_array)
    plt.plot(y_points,z_points)

    plt.show()

if __name__ == "__main__":
    main()
