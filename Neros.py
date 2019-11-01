from scipy.integrate import cumtrapz as cumtrapz
import math
import numpy as np
# Variables needed for other functions
c = 3 * (10**5) # km/s

# VlumSquared - Given Vgas, Vdisk, and Vbulge this will return VSquared
# Params
#  Vgas
#  Vdisk
#  Vbulge
def vLumSquared(vGas, vDisk, vBulge):
  return vGas*vGas + vDisk*vDisk + vBulge*vBulge

def phi(radii, vLum):
  # Square the vlum
  vLumSquared = np.square(vLum)
  # Create an array of just the radii
  x = radii
  # Empty array of the y values that will be added
  y = []
  # Loop over x values and apply the function to get y values
  for i in range(len(x)):
    y.append(vLumSquared[i]/(x[i]*(c*c)))
  x = np.insert(x, 0, 0)
  y = np.insert(y, 0, 0)
  # Calculate the integration
  phi = cumtrapz(y, x)

  return phi

# Kappa - Given phiMilkyWay and phiOtherGalaxy, calculate kappa
# Params
#  phiMW - array of phi values for Milky Way
#  phiOther - array of phi values for other galaxy
def kappa(MW_phi, Other_phi):
  return Other_phi / MW_phi

# Beta - Calculate the beta value for use in E Tsi
# Params
#  VlumOther - Vlum data for another galaxy
def beta(vLumOther):
  return vLumOther / c

# eTsiFlat - Calculate eTsiFlat given Vlum
# Params
#  VlumOther - Vlum data for another galaxy
def eTsiFlat(beta):
  return np.sqrt((1 + beta)/(1 - beta))

# eTsiCurve - Calculate eTsiCurve given Vlum
# Params
#  VlumOther - Vlum data for another galaxy
def eTsiCurve(MW_phi, Other_phi):
  return np.sqrt((1 - 2*MW_phi)/(1 - 2*Other_phi))

def v1(eTsiCurve):
  num = 2
  den = eTsiCurve + 1/eTsiCurve
  return 1 - (num/den)

def v2(eTsiFlat, eTsiCurve):
  num = eTsiFlat + eTsiCurve
  den = eTsiFlat - eTsiCurve
  return num/den

def vLcm(radii, MW_vLum, Other_vLum):
  MW_phi = phi(radii, MW_vLum)
  Other_phi = phi(radii, Other_vLum)
  b = beta(Other_vLum)
  etflat = eTsiFlat(b)
  etCurve = eTsiCurve(MW_phi, Other_phi)
  k = kappa(MW_phi, Other_phi)

  return c*c*k*k*v1(etCurve)*v2(etflat, etCurve)

def vNeros(Other_Vlum, vLCM, freeParam):
  return np.sqrt(np.square(Other_Vlum) + (vLCM*freeParam) )

def chiSquared(observed, expected):
  chiSq = 0
  for i in len(observed):
    chiSquared = chiSquared + (((observed - expected)**2) / expected)
    
  return chiSq