"""
file: 2d-shot.py

Author: Caitlin Welch
Date created: August 24, 2020
Date modified: Auguest 25, 2020

Brief: Uses RK4 ODE solver to plot the 2D trajectory of a Wilson basketball
"""

import numpy as np
import pylab as plt
from scipy.integrate import odeint

C_d = 0.47 #drag coefficient for a sphere
rho = 1.225 #air density in kg/m^3
C = 0.7493 #circumference in m, from basketball's circumference of 29.5 inches
A = C**2/(4*np.pi) #cross-sectional area in m^2
g = 9.81 #gravitational constant on Earth
m = 22/35.274 #gives mass in kg of 22oz basketball

coef = -C_d*rho*A #only used so I don't have to keep writing this out
    
def f(r,t):
    x = r[0]
    y = r[1]
    v_x = r[2]
    v_y = r[3]
    
    fx = v_x
    fy = v_y
    f_v_x = coef*fx*np.sqrt(fx**2+fy**2)/(2*m)
    f_v_y = coef*fy*np.sqrt(fx**2+fy**2)/(2*m)-g
    
    return np.array([fx,fy,f_v_x,f_v_y],float)


def main():
    v0 = float(input("Enter initial velocity (m/s): "))
    theta = (np.pi/180) * float(input("Enter angle shot at in degrees: "))
    init = [0,0,v0*np.cos(theta),v0*np.sin(theta)]
    
    t = np.linspace(0,100,1000)
    
    solution = odeint(f,init,t)
    
    plt.plot(solution[:,0],solution[:,1])
    
    
if __name__ == "__main__":
    main()
