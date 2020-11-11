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
def hit_backboard(phi_array,y,z,backboard_y,backboard_z):
    C = 0.7493 # circumference in m, from basketball's circumference of 29.5 inches
    radius = C/(2*np.pi)
    height_backboard = 1.0668
    eps_backboard = 1e-4
    eps = 1e-1

    theta = np.linspace(0,2*np.pi,100)
    circle_ys = y+radius*np.cos(theta)
    circle_zs = z+radius*np.sin(theta)
    #print(circle_ys,circle_zs)
    #plt.plot(circle_ys,circle_zs)
    #print(phi_array)
    for j in range(len(phi_array)-1):
        #print(len(phi_array)-1)
        #print("j",j)
        phi = phi_array[j]
        #print(phi)
        slope = np.tan(phi_array[j])
        for i in range(len(circle_ys)):
            on_slope = np.absolute(circle_zs[i]-backboard_z[j]-slope*(circle_ys[i]-backboard_y[j]))
            #print("here")
            if phi <= np.pi/2:
                if circle_ys[i] >= backboard_y[j] and circle_ys[i] <= backboard_y[j+1] and circle_zs[i] >= backboard_z[j] and circle_zs[i] <= backboard_z[j+1] and on_slope <= eps:
                    hit_phi = phi
                    #print("hit")
                    return True, hit_phi
            else:
                if circle_ys[i] <= backboard_y[j] and circle_ys[i] >= backboard_y[j+1] and circle_zs[i] >= backboard_z[j] and circle_zs[i] <= backboard_z[j+1] and on_slope <= eps:
                    hit_phi = phi
                    #print("hit")
                    return True, hit_phi
            #print("i after",i)

    hit_phi = -1
    #print("no")
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
def RK4(init,phi_array,backboard_y,backboard_z):
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
    #print("hi")
    while (r[1] >= 3.048 or r[3] > 0) and r[0] >= min_y:
        #print("while")
        #plt.plot(r[0],r[1])
        y_points.append(r[0])
        z_points.append(r[1])
        v_y_points.append(r[2])
        v_z_points.append(r[3])

        backboard, hit_phi = hit_backboard(phi_array,r[0],r[1],backboard_y,backboard_z)
        #print(backboard)
        if (not stop) and backboard:
            stop = True
            y_bounce, z_bounce, v_y_bounce, v_z_bounce = elastic_bounce(hit_phi,r[0],r[1],r[2],r[3])
            r = np.array([y_bounce,z_bounce,v_y_bounce,v_z_bounce],float)
            if in_basket(r[0],r[1]):
                basket = True

        else:
            k1 = h*f(r)
            k2 = h*f(r+0.5*k1)
            k3 = h*f(r+0.5*k2)
            k4 = h*f(r+k3)
            r += (k1+2*k2+2*k3+k4)/6

    if stop:
        #print("hit")
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


def move_points(num_points):
    backboard_y = np.zeros(num_points) #[0,-0.1,0]
    backboard_z = np.zeros(num_points) #[3.048,3.5814,4.1148]
    backboard_z[0] = 3.048

    for i in range(1,num_points):
        backboard_z[i]=(backboard_z[i-1]+(4.1148-3.048)/(num_points-1))
        backboard_y[i]=(np.random.uniform(-0.5,0.5))
    #print(backboard_y,backboard_z)
    plt.plot(backboard_y,backboard_z)
    return backboard_y,backboard_z


def percent_in():
    num_backboards = 2
    num_backboard_points = 3
    num_shots = 1000
    backboard_y = np.zeros((num_backboards,num_backboard_points))
    backboard_z = np.zeros((num_backboards,num_backboard_points))
    shots_in = np.zeros(num_backboards)

    for i in range(num_backboards):
        back_y, back_z = move_points(num_backboard_points)
        backboard_y[i] = back_y
        backboard_z[i] = back_z
        #print(backboard_y,backboard_z)
        #backboard_y[i],backboard_z[i] = move_points(num_backboard_points)
        #plt.plot(backboard_y[i],backboard_z[i])
        del_y = []
        del_z = []
        #print(backboard_y)
        #Creates an array of the change in y and z for each point in the backboard
        for m in range(len(backboard_y[i])-1):
            del_y.append(backboard_y[i][m+1]-backboard_y[i][m])
            del_z.append(backboard_z[i][m+1]-backboard_z[i][m])
        phi_array = np.arctan(np.array(del_z)/np.array(del_y))
        for j in range(len(phi_array)):
            if phi_array[j] < 0:
                phi_array[j] += np.pi

    while np.amax(shots_in)==0:
        random_y = 7.5 #np.random.uniform(0.5,7) #7.5
        random_z = 1.8 #np.random.uniform(1.5,2) #1.8
        random_v0 = 9.8 #np.random.uniform(6,10) #9.8
        random_theta = np.random.uniform(np.pi/2,np.pi) #2*np.pi/3
        init = [random_y,random_z,random_v0*np.cos(random_theta),random_v0*np.sin(random_theta)]

        for i in range(num_backboards):
            y_points,z_points,v_y_points,v_z_points,basket = RK4(init,phi_array,backboard_y[i],backboard_z[i])
            plt.plot(y_points,z_points)

            if basket:
                shots_in[i] += 1
            #print(shots_in[i],type(shots_in[i]),type(num_shots))

    #print(shots_in)
    percent_shots_in = np.array(shots_in)/np.float(num_shots)
    print(percent_shots_in)
    
    return percent_shots_in,backboard_y,backboard_z

    
def best_backboard():
    percent_shots_in,backboard_y_all,backboard_z_all = percent_in()
    max_index = np.argmax(percent_shots_in)
    backboard_y = backboard_y_all[max_index]
    backboard_z = backboard_z_all[max_index]
    
    print(backboard_y,backboard_z)
    

def main():
    best_backboard()
    
    plt.savefig("plot-output.png")
    plt.show()

if __name__ == "__main__":
    main()
