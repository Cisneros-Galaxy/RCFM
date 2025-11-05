# -*- coding: utf-8 -*-
"""
Created on Tue Nov  4 09:13:09 2025

@author: canto
"""

import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
from scipy.optimize import curve_fit

# Helper functions from DataAid.py and DataImport.py
import DataAid
import DataImporter

# Numerically stable class of functions from Neros_v2.py
import Neros

# Load Galaxy Data
sparcTset = DataAid.GetGalaxyData("data/Sparc/TrainingSet/")

# Load Milky Way Model Data
xueSofueGalaxies = DataAid.GetGalaxyData("data/XueSofue/")

# Create array of Milky Way radius and vlum tuples from model data
MWXueSofue = np.array(xueSofueGalaxies['MW_lum'])

# Create Neros instance to perform calculations with the supplied Milky Way model as comparison
# Change Milky Way model by changing the variable in the parentheses
# i.e. neros_fns = Neros_v2.Neros(MWModelVariable)
neros_fns = Neros.Neros(MWXueSofue)
MW_name = "XueSofueGaia2" # Change this if you change the MW model in neros_fns!

MW_rad = neros_fns.mw_rad
MW_vLum = neros_fns.mw_vLum
#MW_phi = neros_fns.mw_phi

MW_vLum_interp_func = neros_fns.mw_vLum_interp

# This designates which galaxy sample to fit
#galaxies = sparcGalaxies
galaxies = sparcTset

galaxyName = 'NGC2976_rotmod'

# Extract out the needed galaxy components
galaxy = np.array(galaxies[galaxyName])
galaxy_rad = galaxy[:,0]
galaxy_vObs = galaxy[:,1]
galaxy_error = galaxy[:,2]
galaxy_gas = galaxy[:,3]
galaxy_disk = galaxy[:,4]
galaxy_bulge = galaxy[:,5]

neros_fns.fit(galaxy_rad, galaxy_gas, galaxy_disk, galaxy_bulge, galaxy_vObs, galaxy_error)
fit_results = neros_fns.get_fit_results(galaxy_rad)