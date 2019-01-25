import matplotlib.pyplot as plt
import numpy as np
import os
import math
from os import listdir
import scipy.interpolate as interpolate



DATA_DIR = './Data/'

# Import galaxy data
def import_data():
    data_list = []
    for item in os.listdir(DATA_DIR):
        if not item.startswith('.'):
            data_list.append(item)
    return data_list


# Load data from files into matrix [[[radii], [velocities]]]
def get_data(ycol):
    data_matrix = []
    data_list = import_data()
    for file in data_list:
        with open(DATA_DIR + file, 'r') as f:
            next(f)
            x, y = [], []
            for line in f:
                values = [float(s) for s in line.split(   )]
                # Remove gas, bulge data for galaxies with no recorded values
                # if values[ycol] != 0:
                x.append(values[0])
                y.append(values[ycol])
            data_matrix.append([x, y])
    return data_matrix


# Plot raw galaxy data
def plot_raw_data(ycol, xlabel, ylabel, title):
    data_matrix = get_data(ycol)
    plt.clf()
    for i in range(len(data_matrix)):
        x = data_matrix[i][0]
        y = data_matrix[i][1]
        plt.plot(x, y, marker='o', markersize=2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()


# Scale radius such that 0<r<1. The list scaled_radii preserves groupings of
# radii by galaxy. Index [0]-[50] corresponds to the radii of the
# respective file number
def scale_radii(ycol):
    scaled_radii = []
    ycol_data_matrix = get_data(ycol)
    for i in range(len(ycol_data_matrix)):
        radii_values = ycol_data_matrix[i][0]
        max_radius = max(radii_values)
        scaled_galaxy_radii = []
        for j in range(len(radii_values)):
            scaled_radius = radii_values[j] / max_radius
            scaled_galaxy_radii.append(scaled_radius)
        scaled_radii.append(scaled_galaxy_radii)
    return scaled_radii


# Calculates luminous velocity for each galaxy using
# v_lum_val = sqrt(vbulge^2 + vdisk^2 +vgas^2)
def calc_luminous_velocity(ycol1, ycol2, ycol3):
    v_lum = []
    v_disk_matrix = get_data(ycol1)
    v_gas_matrix = get_data(ycol2)
    v_bulge_matrix = get_data(ycol3)
    for i in range(len(v_disk_matrix)):
        v_lum_list = []
        for j in range(len(v_disk_matrix[i][1])):
            v_disk_sqrd = v_disk_matrix[i][1][j]**2
            v_gas_sqrd = v_gas_matrix[i][1][j]**2
            v_bulge_sqrd = v_bulge_matrix[i][1][j]**2
            v_lum_val = math.sqrt(v_disk_sqrd + v_gas_sqrd + v_bulge_sqrd)
            v_lum_list.append(v_lum_val)
        v_lum.append(v_lum_list)
    return v_lum


# Function to plot v_lum data vs. scaled r values
def plot_scaled_data():
    scal_rad = scale_radii(4)
    v_lum_vals = calc_luminous_velocity(4, 5, 6)
    v = []
    for i in range(len(scal_rad)):
        r = scal_rad[i]
        v = v_lum_vals[i]
        plt.plot(r, v, marker='o', markersize=2)
    plt.xlabel('Radius (scaled)')
    plt.ylabel('Luminous Velocity (km/s)')
    plt.title('Scaled Data')
    plt.show()


def Bspline_interpolate_v_lum():
    plt.clf()

    # Produce scaled radii and luminous velocity data
    scal_rad = scale_radii(4)
    v_lum = calc_luminous_velocity(4, 5, 6)

    for i in range(len(scal_rad)):
        r = np.array(scal_rad[i])
        v = np.array(v_lum[i])
        N = len(r)
        rmin, rmax = r.min(), r.max()
        rr = np.linspace(rmin, rmax, N)

        # Creat tuple (t, c, k) containing vector of knots,
        # B-spline coefficients, and degree of spline
        t, c, k = interpolate.splrep(r, v, s=0, k=4)


        # Create B-spline object using (t, c, k)
        spline = interpolate.BSpline(t, c, k, extrapolate = False)

        # Plot B-spline fit (functional data)
        plt.plot(r, spline(rr), label = 'BSpline')

    plt.xlabel('Radius (scaled)')
    plt.ylabel('Luminous Velocity (km/s)')
    plt.title('Functional Data')
    plt.show()


# Function to visualize raw data. ycol corresponds to velocity data:
# observed, disk, gas, or bulge.
def visualize_raw_data(ycol, y_label, title):
    x_label = 'Radius (kpc)'
    plot_raw_data(ycol, x_label, y_label, title)


def main():
    # Close all figure windows to free memory before initialization
    plt.close('all')

    visualize_raw_data(2, 'Observed Velocity (km/s)', 'Raw Data: Observed Velocity')
    visualize_raw_data(4, 'Disk Velocity (km/s)', 'Raw Data (Disk Velocity)')
    visualize_raw_data(5, 'Gas Velocity (km/s)', 'Raw Data (Gas Velocity)')
    visualize_raw_data(6, 'Bulge Velocity (km/s)', 'Raw Data (Bulge Velocity)')

    # Scaled luminous velocity data
    plot_scaled_data()

    # Convert raw luminous velocity data to functional data 
    # using B spline interpolation
    Bspline_interpolate_v_lum()


main()
