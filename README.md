# Hoop That Never Misses
This repository contains files for finding the ideal basketball backboard shape so that the ball will always go into the basket. While the code in this repository does not provide the final result, the files provide an outline for the Monte Carlo optimization process required for finding the backboard shape.

## File Structure
* Python files: basic scripts that plot the trajectory of the basketball and run the Monte Carlo optimization for the backboard shape.
    * monte_carlo_two_d.py - for a given number of points used to compose the backboard, a given number of random shots to test and a given number of random backboards composed of the specified number of points to compare, returns the backboard that results in the most shots going in the basket.
    * Backboard Plot
        * plot_backboard.py - plots a backboard given the coordinates for the backboard points. monte_carlo_two_d.py returns the y and z coordinates of the optimal backboard, so plot_backboard.py is used to plot those backboard points to show the backboard shape.
    * Plotting Shots
        * two_d_shot.py - solves for and plots the trajectory of a basketball shot, including hitting the backboard in an elastic collision and checking whether the ball goes in the hoop. The two dimensions used in this file are y and z - the distance out from the backboard towards the half court line and height.
        * two_d_shot_birdseye.py - same as two_d_shot.py, but the two dimensions used in this file are x and y - the distance from the center of the hoop in terms of the width of the court and the distance out from the backboard towards the half court line.
        * three_d_shot.py - same as two_d_shot.py but includes all three dimensions.
* Shell scripts:
    * monte_carlo.sh - SLURM shell script used to run monte_carlo_two_d.py 
* png files: found in Weekly Presentations > Photos. Plots saved from running each of the python files for the trajectories of the basketball shots and the different backboard shapes.
* Powerpoints: weekly presentations for research meetings and final research presentation.

## Dependencies
All python code written for this project is in python3. The packages used in the code are **numpy** and **pylab**.

To effectively use the files, you can clone the repository and download **numpy** and **pylab** through your command line:
`pip install numpy pylab`

## Running the Code
As noted in *File Structure*, *monte_carlo.sh* will be used to run the main code file, *monte_carlo_two_d.py*. The file should be run using SLURM, to account for the computationally intensive code and long run times. 

In order to run the code, first ensure you are in the correct directory in Jupyter Hub. Then, use the following command in the Jupyter command line:
`sbatch monte_carlo.sh`

To change the number of simulations, you must edit the *monte_carlo_two_d_.py* directly. 

## Important Information Pertaining to the Project
* 