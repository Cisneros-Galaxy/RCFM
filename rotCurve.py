import math

c = 3*(10**8)
overc = 1/c

def vRot(radius, vSquared):
  vRot = math.sqrt()

def kappa(phi_mw, phi_other, phi_o):
  return (phi_other-phi_o)/(phi_mw-phi_o)

def vlum(vgas, vdisk, vbul):
  return math.sqrt(vgas*vgas + vdisk+vdisk + vbul*vbul)

def v1(phi_mw, phi_other):
  return (1-(2/(math.sqrt((1 - phi_mw)/(1 - phi_other) ) + math.sqrt((1-phi_other)/(1 -
phi_mw)))))

def v2():
  return (math.sqrt((1+vlum(r_other_galaxy)*overc)/(1.-vlum(r)*overc))+sqrt((1 -phi(r_mw_i))/(1 -
phi(r_other_galaxy)) ))/(sqrt((1.+vlum(r_other_galaxy)*overc)/(1.-vlum(r_other_galaxy)*overc))-
sqrt((1 -phi(r_mw_i))/(1 -phi(r_other_galaxy)) ) )

def vLCM(phi_mw, phi_other, phi_o, alpha):
  return c*c*kappa()