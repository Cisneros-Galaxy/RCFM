from scipy.integrate import quad as integrate
import math

# Variables needed for other functions
c = 3 * (10**5) # km/s

# VlumSquared - Given Vgas, Vdisk, and Vbulge this will return VSquared
# Params
#  Vgas
#  Vdisk
#  Vbulge
def VlumSquared(Vgas, Vdisk, Vbulge):
  return Vgas*Vgas + Vdisk*Vdisk + Vbulge*Vbulge

# CalcPhi - Given radii and vSquared, this will return Phi at each radius
# Params
#  radii    - an array of all the radii to calculate at
#  VlumSquared - an array of VlumSquared at each radius
def Phi(radii, VlumSquared):
  phi = []
  phi.append(0)
  for i in range(len(radii)):
    if i == 0:
        continue

    # Get min radius. Will be 0 for the first time, then the previous r for every other
    rMin = radii[i-1]
        
    # Max radius will always be the radius at current index
    rMax = radii[i]
    vS = VlumSquared[i]

    # New Phi data points are cumulative, so we need to add 
    #   the previous point (if it exists)
    newPhiDataPoint = phi[i-1]
    
    # Calculate the integral at this point
    integral, err = integrate(lambda r, vS: (vS)/(r*(c*c)), rMin, rMax, vS)

    # Add it to the newPoint (essentially adding to the sum)
    newPhiDataPoint += integral

    # Append the new point to the array!
    phi.append(newPhiDataPoint)
    
  phi = phi[1:]
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
def eTsiCurve(phiMW, phiOther):
  return ((-2*phiMW+1)/(-2*phiOther + 1))**(1.0/2.0)

def v1(eTsiCurve):
  num = 2
  den = eTsiCurve + 1/eTsiCurve
  return (num/den) - 1

def v2(eTsiFlat, eTsiCurve):
  num = eTsiFlat + eTsiCurve
  den = eTsiFlat - eTsiCurve
  return num/den

def Vlcm(radii, MW_phi, Other_phi, MW_Vlum, Other_Vlum):
  MW_phi = Phi(radii, MW_Vlum)
  Other_phi = Phi(radii, Other_Vlum)
  b = beta(Other_Vlum)
  etflat = eTsiFlat(b)
  etCurve = eTsiCurve(MW_phi, Other_phi)

  return c*c*kappa(MW_phi, Other_phi)*v1(etCurve)*v2(etflat, etCurve)