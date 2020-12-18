# This is an attempt at a more numerically stable, cleaned-up
# version of the Neros calculations. 
# The big changes are that it uses the calculations in NumericalStability.ipynb
# and uses and interpolated MW phi to avoid data issues
# I'm going to implement it as a class, instead
# of a namespace of functions, which should make using it easier
# Particularly, since we'll be using the same Milky Way data
# to do a family of fits, that will be a class member
# Usage will be documented in the class

import numpy as np
from scipy.interpolate import interp1d
import scipy.integrate


c = 3 * (10**5) # km/s

class Neros:
    """The Neros Model
    
    This class implements the actual computations of
    the Neros model. Since the model depends on the Milky Way
    mass curve, each instance of this class will need to be 
    supplied with one.  Meaning that if you want to do multiple
    Milky Way models, each will need its own instance.
    
    Create an instance by calling
    Neros(milky_way_data)
    
    The data can be in the form of a NumPy array (of shape N x 2),
    a Pandas DataFrame with two columns, or a list of tuples.
    Look at the DataImporter.py file for help reading from files 
    if needed. The data must come in the order
    radius, vLum (it ignores Pandas column names).
    
    The Milky Way data can be changed later by calling setMilkyWay
    """
    
    def __init__(self, milky_way_data):
        #This will check if the Milky Way data
        #is properly formatted, and raise an exception
        #if it isn't
        self.setMilkyWay(milky_way_data)

    def setMilkyWay(self, milky_way_data):
        """Changes the internal Milky Way
        
        Data must be two columns in the order radius, vLum, 
        either as a NumPy array, a Pandas DataFrame, or a
        list of tuples
        
        This automatically sets up the appropriate interpolators"""
        
        data = np.array(milky_way_data)
        
        if len(data.shape) != 2 or data.shape[0] < 2 or data.shape[1] != 2:
            raise ValueError("Milky Way data was not in the form of two columns")
        
        #Retain the data, just in case, but we'll mostly be doing interpolation
        self.milky_way_data = data
        self.mwRad = data[:,0]
        self.mwVLum = data[:,1]
        self.mw_vLum_interp = interp1d(data[:,0], data[:,1], kind='cubic')
        mw_phi = self.phi(data[:,0], data[:,1])
        self.mw_phi_interp = interp1d(data[:,0], mw_phi)
        
    def vLumSquared(self, vGas, vDisk, vBulge, disk_scale=1, bulge_scale=1):
        """Calculates total luminous velocity from the sum of the squares of
           the contributions from gas, disk, and bulge"""
        return vGas*vGas + disk_scale*vDisk*vDisk + bulge_scale*vBulge*vBulge
    
    def curve_fit_fn(self, galaxyData, alpha, vLumFreeParam, disk_scale, bulge_scale):
        """Formerly known as 'simple'.
        This is used as the fitting function in scipy.curve_fit
        to find alpha and vLumFreeParam, and it calculates vNeros 
        using those values. 
        
        The parameters are
        :galaxyData: A 4-D Numpy Array of radii and vGas, vDisk, VBulge data
        :alpha: fitting parameter from scipy.curve_fit
        :vLumFreeParam: fitting parameter from scipy.curve_fit"""
        
        # Parse out data for the galaxy
        rad,vGas,vDisk,vBulge = galaxyData
        # Apply the vLum free param to the data
        vLum_scaled = self.vLumSquared(vGas,vDisk,vBulge,disk_scale,bulge_scale)
        # Calc and return vNeros
        return self.vNeros(rad, vLum_scaled, alpha)
        
    def phi(self, radius, vlum):
        """Computes potential. phi = integrate vlum^2/r/c^2"""
        
        x = radius
        y = np.square(vlum) / (x*c*c)
        x = np.concatenate([np.array([0]), x])
        y = np.concatenate([np.array([0]), y])
        #May want to switch to Simpson's rule: scipy.integrate.simps
        return scipy.integrate.cumtrapz(y,x)
    
    def vNeros(self, galaxy_rad, galaxy_vLum, alpha, phi_zero=3e-11):
        """This computes the predicted vObs
        
        The parameters are
        :galaxy_rad: A 1-D NumPy array or Pandas DataSeries of radii
        :galaxy_vLum: A 1-D NumPy array or Pandas DataSeries of vLums
        :alpha: From the equation vObs^2 = vLum^2 + alpha*vLCM^2
        :phi_zero: The zero point for the phi integration, used in kappa calculation
        
        This calls sqrt internally, so it can fail for some values of alpha."""
        
        return np.sqrt(self.vNerosSquared(galaxy_rad, galaxy_vLum, alpha, phi_zero=3e-11))
        
    def vNerosSquared(self, galaxy_rad, galaxy_vLum, alpha, phi_zero=3e-11):
        """This computes the predicted vObs^2, basically by calling vLCM and applying alpha
        
        The parameters are
        :galaxy_rad: A 1-D NumPy array or Pandas DataSeries of radii
        :galaxy_vLum: A 1-D NumPy array or Pandas DataSeries of vLums
        :alpha: From the equation vObs^2 = vLum^2 + alpha*vLCM^2
        :phi_zero: The zero point for the phi integration, used in kappa calculation
        
        This avoids the issues inherent in calling sqrt, but can produce non-physical vObs^2"""
        
        vLCM = self.vLCM(galaxy_rad, galaxy_vLum, phi_zero)
        # this formula allows vNerosSquared to be negative for some fits which is problematic
        # chisquared becomes NaN
        return galaxy_vLum*galaxy_vLum + alpha*vLCM
    
    def vLCM(self, galaxy_rad, galaxy_vLum, phi_zero=3e-11):
        """This computes the vLCM - the actual model
        
        The parameters are
        :galaxy_rad: A 1-D NumPy array or Pandas DataSeries of radii
        :galaxy_vLum: A 1-D NumPy array or Pandas DataSeries of vLums
        :phi_zero: The zero point for the phi integration, used in kappa calculation"""
        
        MW_phi = self.mw_phi_interp(galaxy_rad)
        galaxy_phi = self.phi(galaxy_rad, galaxy_vLum)
        k = self.kappa(MW_phi, galaxy_phi, phi_zero)
        v1 = self.v1(MW_phi, galaxy_phi)
        v2 = self.v2(MW_phi, galaxy_phi, galaxy_vLum)
        vLCM = c * c * k * k * v1 * v2
        
        return vLCM
    
    def kappa(self, MW_phi, other_phi, phi_zero=3e-11):
        """kappa(r) in the paper, just phi_gal(r)/phi_mw(r)"""
        
        return (other_phi - phi_zero) / (MW_phi - phi_zero)
    
    def v1(self, MW_phi, other_phi):
        etc = self._eTsiCurveMinusOne(MW_phi, other_phi)
        num = etc**2
        den = (1 + etc)**2 + 1
        return num/den
    
    def v2(self, MW_phi, other_phi, other_vlum):
        etFlat = self._eTsiFlatMinusOne(other_vlum)
        etCurve = self._eTsiCurveMinusOne(MW_phi, other_phi)
        num = 2 + etFlat + etCurve
        den = etFlat - etCurve
        return num/den

    
    def chiSquared(self, observed, expected, error):
        """This computes chi squared"""
        chiSquared = 0

        for i in range(len(observed)):
            chiSquared = chiSquared + (((observed[i] - expected[i])**2) / (error[i]**2))
        return chiSquared / len(observed)
    
    def _eTsiFlatMinusOne(self, other_vlum):
        """This computes eTsiFlat - 1, compared to the old code, for numerical stability"""
        
        beta = other_vlum / c
        numerator = 2*beta / (1 - beta)
        denominator = np.sqrt((1+beta) / (1-beta)) + 1
        return numerator / denominator
    
    def _eTsiCurveMinusOne(self, MW_phi, other_phi):
        """This computes eTsiFlat - 1, compared to the old code, for numerical stability"""

        numerator = (2*other_phi - 2*MW_phi) / (1 - 2*other_phi)
        denominator = np.sqrt((1 - 2*MW_phi) / (1 - 2*other_phi)) + 1
        return numerator / denominator
