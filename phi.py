from scipy.integrate import quad as integrate

# Variables needed for other functions
c = 3 * (10**5) # km/s

# This is the function for calc Phi.
# NOTE: This function is to be used with Scypi.integrate.quad - when doing this you can specify
#   arguments to pass in, but it will always pass in the arguments AFTER x (or r in our case)
def phiFunc(r, vSquared):
    return (vSquared)/(r*(c*c))

# CalcPhi - Given radii and vSquared, this will return Phi at each radius
# Params
#  radii    - an array of all the radii to calculate at
#  vSquared - an array of vSquared at each radius
def CalcWithRadiiAndVSquared(radii, vSquared):
  phi = []
  phi.append(0)
  for i in range(len(radii)):
    if i == 0:
        continue

    # Get min radius. Will be 0 for the first time, then the previous r for every other
    rMin = radii[i-1]
        
    # Max radius will always be the radius at current index
    rMax = radii[i]
    vS = vSquared[i]

    # New Phi data points are cumulative, so we need to add 
    #   the previous point (if it exists)
    newPhiDataPoint = phi[i-1]
    
    # Calculate the integral at this point
    integral, err = integrate(phiFunc, rMin, rMax, vS)

    # Add it to the newPoint (essentially adding to the sum)
    newPhiDataPoint += integral

    # Append the new point to the array!
    phi.append(newPhiDataPoint)
    
  phi = phi[1:]
  return phi