from scipy.integrate import quad as integrate
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

def Phi(radii, VlumSquared):
  phi = []
  radii = np.insert(radii, 0, 0)
  for i in range(len(radii)-1):

    # Get min radius. Will be 0 for the first time, then the previous r for every other
    rMin = radii[i]
        
    # Max radius will always be the radius at current index
    rMax = radii[i+1]
    vS = VlumSquared[i]

    # New Phi data points are cumulative, so we need to add 
    #   the previous point (if it exists)
    newPhiDataPoint = 0
    if i != 0:
      newPhiDataPoint = phi[i-1]
    
    # Calculate the integral at this point
    integral, err = integrate(lambda r, vS: (vS)/(r*(c*c)), rMin, rMax, vS)

    # Add it to the newPoint (essentially adding to the sum)
    newPhiDataPoint += integral

    print("---Integrating----")
    print(rMin, " -> ", rMax)
    print("With constant v^2: ", vS)
    print("Value: ", integral)

    # print(rMin, rMax, vS, vS/(rMin*(c*c)), newPhiDataPoint)

    # Append the new point to the array!
    phi.append(newPhiDataPoint)
    
  # phi = phi[1:]
  # print(phi)
  return np.array(phi)

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
def eTsiCurve(phiMW, phiOther):
  return np.sqrt((-2*phiMW+1)/(-2*phiOther + 1))

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
  MW_phi = Phi(radii, MW_Vlum*MW_Vlum)
  print(MW_phi)
  print("[CALCULATING PHI FOR OTHER]")
  Other_phi = Phi(radii, Other_Vlum*Other_Vlum)
  print(Other_phi)
  b = beta(Other_Vlum)
  etflat = eTsiFlat(b)
  etCurve = eTsiCurve(MW_phi, Other_phi)
  k = kappa(MW_phi, Other_phi)

  return c*c*k*v1(etCurve)*v2(etflat, etCurve)