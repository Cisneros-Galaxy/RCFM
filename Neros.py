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
def VlumSquared(Vgas, Vdisk, Vbulge):
  return Vgas*Vgas + Vdisk*Vdisk + Vbulge*Vbulge

def Phi(radii, Vlum):
  # Square the vlum
  VlumSquared = np.square(Vlum)
  # Create an array of just the radii
  x = radii
  # Empty array of the y values that will be added
  y = []
  # Loop over x values and apply the function to get y values
  for i in range(len(x)):
    y.append(VlumSquared[i]/(x[i]*(c*c)))
  x = np.insert(x, 0, 0)
  y = np.insert(y, 0, 0)
  # Calculate the integration
  phi = cumtrapz(y, x)

  # Printing to check
  print("[X VALUES TO INTEGRATE]")
  print(len(x))
  print(x)
  print("[Y VALUES TO INTEGRATE]")
  print(len(y))
  print(y)
  print("[PHI VALUES]")
  print(len(phi))
  print(phi)

  return phi

# Kappa - Given phiMilkyWay and phiOtherGalaxy, calculate kappa
# Params
#  phiMW - array of phi values for Milky Way
#  phiOther - array of phi values for other galaxy
def kappa(phiMW, phiOther):
  return phiOther / phiMW

# Beta - Calculate the beta value for use in E Tsi
# Params
#  VlumOther - Vlum data for another galaxy
def beta(VlumOther):
  return VlumOther / c

# eTsiFlat - Calculate eTsiFlat given Vlum
# Params
#  VlumOther - Vlum data for another galaxy
def eTsiFlat(beta):
  return ((beta+1)/(-beta+1))**(1.0/2.0)

# eTsiCurve - Calculate eTsiCurve given Vlum
# Params
#  VlumOther - Vlum data for another galaxy
def eTsiCurve(MW_phi, Other_phi):
  return np.sqrt((-2*MW_phi+1)/(-2*Other_phi + 1))

def v1(eTsiCurve):
  num = 2
  den = eTsiCurve + 1/eTsiCurve
  return 1 - (num/den)

def v2(eTsiFlat, eTsiCurve):
  num = eTsiFlat + eTsiCurve
  den = eTsiFlat - eTsiCurve
  return num/den

def Vlcm(radii, MW_Vlum, Other_Vlum):
  print("[CALCULATING PHI FOR MW]")
  MW_phi = Phi(radii, MW_Vlum)
  print(MW_phi)
  print("[CALCULATING PHI FOR OTHER]")
  Other_phi = Phi(radii, Other_Vlum)
  print(Other_phi)
  b = beta(Other_Vlum)
  etflat = eTsiFlat(b)
  etCurve = eTsiCurve(MW_phi, Other_phi)
  k = kappa(MW_phi, Other_phi)

  return c*c*k*v1(etCurve)*v2(etflat, etCurve)