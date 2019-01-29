from scipy.integrate import cumtrapz as cumtrapz

# Variables needed for other functions
c = 3 * (10**5) # km/s

# This is the function for calc Phi.
# NOTE: This function is to be used with Scypi.integrate.quad - when doing this you can specify
#   arguments to pass in, but it will always pass in the arguments AFTER x (or r in our case)
def force(r, vSquared):
    return (2*vSquared)/(r*(c*c))

# CalcPhi - Given radii and vSquared, this will return Phi at each radius
# Params
#  radii    - an array of all the radii to calculate at
#  vSquared - an array of vSquared at each radius
def CalcWithRadiiAndVSquared(radii, vSquared):
  # Make sure radii and vsquared are same size
  if (len(radii) != len(vSquared)):
    raise ValueError('CalcWithRadiiAndVSquared: Radii and vSquared parameters must be the same size.')

  # Calculate phi for every radius
  phi = []
  for i in range(len(radii)):
    phi.append(force(radii[i], vSquared[i]))
  
  # Calculate cumulative integration
  integration = cumtrapz(phi, radii)

  return integration