"""
file: monte-carlo-2d-shot.py

Author: Caitlin Welch
Date created: September 6, 2020
Date modified: September 6, 2020

Brief: First attempt at running a Monte Carlo simulation for a 2D birdseye shot
        to form a continuous curved backboard
"""

import two_d_shot_birdseye
import numpy as np

def backboard(x_min,x_max,y_min,y_max,phi,x_points,y_points,v_x_points,v_y_points,T):



def change_backboard(x_points,y_points,v_x_points,v_y_points):
    phis = np.linspace(0,90,1000) #angle the backboard makes with the horizontal
    shots_and_makes = {} #this dictionary will hold all the phi angles as keys and number of shots that are made for that angle as values

    for phi in phis:
        hits_backboard,x_backboard,y_backboard,v_x_backboard,v_y_backboard = backboard()

        if hits_backboard:



def monte_carlo():
    start_x = np.linspace(-7.5,7.5,1000)
    start_y = np.linspace(0,10,1000)
    angle = np.linspace(40,85,1000)
    theta = []
    for a in angle:
        theta.append((np.pi/180)*a)
    v0 = np.linspace(4,10,1000)




def main():
    """
    Constants for linspace
    """
    a = 0
    b = 100
    N = 1000
    T = np.linspace(a,b,N)
    h = (b-a)/N

    start_x =
    start_y =
    angle =
    theta = (np.pi/180) * angle
    v0 =

    init = [start_x, start_y, v0*np.cos(theta)*np.sin(phi), v0*np.cos(theta)*np.cos(phi)]
    x_points, y_points, v_x_points, v_y_points = two_d_shot_birdseye.RK4(init,T)



if __name__ == "__main__":
    main()
