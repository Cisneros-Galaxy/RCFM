from scipy.integrate import quad as integrate

# This is the function for calc Phi.
# NOTE: This function is to be used with Scypi.integrate.quad - when doing this you can specify
#   arguments to pass in, but it will always pass in the arguments AFTER x (or r in our case)
def phiFunc(r, v):
    return (v)/(r*(c*c))

# Variables needed for other functions
c = 3 * (10**5) # km/s

# CalcPhi - Given radii and vSquared, this will return Phi at each radius
# Params
#  radii    - an array of all the radii to calculate at
#  vSquared - an array of vSquared at each radius
def CalcPhi(radii, vSquared):
  phi = []
  for i in range(len(radii)):
    # Get min radius. Will be 0 for the first time, then the previous r for every other
    rMin = 0
    if i != 0:
        rMin = radii[i-1]
        
    # Max radius will always be the radius at current index
    rMax = radii[i]
    vS = vSquared[i]

    # New Phi data points are cumulative, so we need to add 
    #   the previous point (if it exists)
    newPhiDataPoint = 0
    if i != 0:
        newPhiDataPoint += phi[i-1]
    
    # Calculate the integral at this point
    integral, err = integrate(phiFunc, rMin, rMax, (vS))
    # Add it to the newPoint (essentially adding to the sum)
    newPhiDataPoint += integral
    # Append the new point to the array!
    phi.append(newPhiDataPoint)
    
  return phi