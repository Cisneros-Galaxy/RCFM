# Luminous Convolution Model

About
-----
This repository holds projects in Python for my work on the Luminous Convolution Model for spiral galaxies. The main program for the time being is lcm_graphs.py, which reads files containing the rotation curve data of 50 galaxy samples and plots them using raw and scaled radii values. It can then be visualized to compare galaxy rotation curves side by side. 

Dependencies
-----
* Python 3 (otherwise interpolate.BSpline will NOT work)
* scipy 
* matplotlib

Development Road Map
-----
* Add interpolation of scaled data to produce functional data (Added 1/18/19)
* Identify landmarks (first peak) of each rotation curve
* Use warping function to align the peaks of each rotation curve to further aid in visual comparison of data.
