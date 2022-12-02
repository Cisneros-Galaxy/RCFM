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
from scipy.optimize import curve_fit


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
        self.mw_rad = data[:,0]
        self.mw_vLum = data[:,1]
        self.mw_vLum_interp = interp1d(data[:,0], data[:,1], kind='cubic')
        self.mw_phi = self.phi(data[:,0], data[:,1])
        self.mw_phi_interp = interp1d(data[:,0], self.mw_phi)


    def vLumSquared(self, vGas, vDisk, vBulge, disk_scale=1, bulge_scale=1):
        """Calculates total luminous velocity from the sum of the squares of
           the contributions from gas, disk, and bulge"""
        return vGas**2 + (disk_scale*vDisk)**2 + (bulge_scale*vBulge)**2


    def fit(self, rad, vGas, vDisk, vBulge, vObs, vObsError):
        """Fits a galaxy using the LCM model, core functionality of this class
        
        This takes a lot of the things that were done ad-hoc in Model.ipynb and
        puts them together to simplify the process.  Results will be stored 
        internally to simplify later operations. Ideally, the usage simplifies
        to simply calling this fit, then asking for chi_squared, vLum, etc
        
        Parameters:
        :rad: The radii of the galaxy being fit, as a numpy array
        :vGas: Inferred gas mass as a velocity, galaxy_vGas, as a numpy array
        :vDisk: Inferred disk mass as a velocity, galaxy_vDisk, as a numpy arry
        :vBulge: Inferred bulge mass as a velocity, galaxy_vBulge, as a numpy array
        :vObs: Observed galaxy rotation velocity, as a numpy array
        :vObsError: Measurement error on vObs, as a numpy array"""
        
        # First we need to clip the galaxy data so it doesn't extend beyond
        # the range of our Milky Way data, we may improve this method later
        # We're storing these for later computations, they'll get overwritten
        # every time we call fit
        valid_rad = rad <= self.mw_rad[-1]
        self.rad = rad[valid_rad]
        self.vGas = vGas[valid_rad]
        self.vDisk = vDisk[valid_rad]
        self.vBulge = vBulge[valid_rad]
        self.vObs = vObs[valid_rad]
        self.vObsError = vObsError[valid_rad]
        
        fit_vals, cov = curve_fit(self.curve_fit_fn,(self.rad, self.vGas, self.vDisk, self.vBulge),
                          self.vObs, sigma=self.vObsError, maxfev=10000)
        
        fit_parameter_names  = ['alpha', 'disk_scale', 'bulge_scale']
        self.best_fit_values = dict(zip(fit_parameter_names, fit_vals))

    
    def curve_fit_fn(self, galaxyData, alpha, disk_scale, bulge_scale):
        """Formerly known as 'simple'.
        This is used as the fitting function in scipy.curve_fit
        to find alpha and vLumFreeParam, and it calculates vNeros 
        using those values. 
        
        The parameters are
        :galaxyData: A 4-D Numpy Array of radii and vGas, vDisk, vBulge data
        :alpha: fitting parameter from scipy.curve_fit, now defined as sqrt(old_alpha)
        :disk_scale: scaling of vDisk
        :bulge_scale: scaling of vBulge"""
        
        # Parse out data for the galaxy
        rad,vGas,vDisk,vBulge = galaxyData

        # Apply the vLum free param to the data
        vLum_scaled = np.sqrt(self.vLumSquared(vGas,vDisk,vBulge,disk_scale,bulge_scale))

        # Calc and return vNeros
        return self.vNeros(rad, vLum_scaled, alpha)


    def get_fit_results(self, galaxy_rad, old_alpha=True):
        """Returns the numerical fit results: chi^2 and best fit parameters
        
        Parameters:
        :old_alpha: Whether to return new_alpha^2 to match old format"""
        
        if not hasattr(self, 'best_fit_values'):
            raise RuntimeError("Please call fit before trying to get the fit results")
        
        chi_squared = self.get_chi_squared()
        alpha = self.best_fit_values['alpha']
        disk_scale = self.best_fit_values['disk_scale']
        bulge_scale = self.best_fit_values['bulge_scale']
        phi_zero = self.get_phi_zero(galaxy_rad)
        
        if old_alpha:
            alpha = alpha**2
            disk_scale = abs(disk_scale)
            bulge_scale = abs(bulge_scale)

        return {'chi_squared': chi_squared, 'alpha': alpha, 'disk_scale': disk_scale, 'bulge_scale': bulge_scale, 'phi_zero': phi_zero}


    def get_chi_squared(self):
        if not (hasattr(self, 'vObs') and hasattr(self, 'vObsError')):
            raise RuntimeError("Please call fit before trying to get the fit results")

        prediction = self.get_vNeros()
        return self.chiSquared(prediction, self.vObs, self.vObsError)


    def get_vLum_scaled(self):
        if not hasattr(self, 'best_fit_values'):
            raise RuntimeError("Please call fit before trying to get the fit results")
        if not (hasattr(self, 'vGas') and hasattr(self, 'vDisk') and hasattr(self, 'vBulge')):
            raise RuntimeError("Please call fit before trying to get the fit results")
        
        disk_scale = self.best_fit_values['disk_scale']
        bulge_scale = self.best_fit_values['bulge_scale']
        vLum_scaled = np.sqrt(self.vLumSquared(self.vGas, self.vDisk, self.vBulge, disk_scale, bulge_scale))
        return vLum_scaled


    def get_vNeros(self):
        if not hasattr(self, 'best_fit_values'):
            raise RuntimeError("Please call fit before trying to get the fit results")
        alpha = self.best_fit_values['alpha']
        vLum_scaled = self.get_vLum_scaled()
        prediction = self.vNeros(self.rad, vLum_scaled, alpha)
        return prediction


    def get_rad(self):
        """Helper function to get the trimmed galaxy radii
        Alteratively, could just call (this object).rad"""
        return self.rad


    def get_vObs(self):
        """Helper function to get the trimmed observation values
        Alteratively, could just call (this object).vObs"""
        return self.vObs


    def get_vObsError(self):
        """Helper function to get the trimmed observation value errors
        Alteratively, could just call (this object).vObsError"""
        return self.vObsError
    
    def get_phi_zero(self, galaxy_rad):
        valid_rad = valid_rad = galaxy_rad <= self.mw_rad[-1]
        trimmed_phi = self.mw_phi[:len(valid_rad)]
        return trimmed_phi[-1]


    def phi(self, radius, vlum):
        """Computes potential. phi = integrate vlum^2/r/c^2"""
        
        x = radius
        y = np.square(vlum) / (x*c*c)
        x = np.concatenate([np.array([0]), x])
        y = np.concatenate([np.array([0]), y])
        #May want to switch to Simpson's rule: scipy.integrate.simps
        return scipy.integrate.cumtrapz(y,x)


    def vNeros(self, galaxy_rad, galaxy_vLum, alpha):
        """This computes the predicted vObs
        
        The parameters are
        :galaxy_rad: A 1-D NumPy array or Pandas DataSeries of radii
        :galaxy_vLum: A 1-D NumPy array or Pandas DataSeries of vLums
        :alpha: From the equation vObs^2 = vLum^2 + alpha*vLCM^2
        
        This calls sqrt internally, so it can fail for some values of alpha."""
        return np.sqrt(self.vNerosSquared(galaxy_rad, galaxy_vLum, alpha))


    def vNerosSquared(self, galaxy_rad, galaxy_vLum, alpha):
        """This computes the predicted vObs^2, basically by calling vLCM and applying alpha
        
        The parameters are
        :galaxy_rad: A 1-D NumPy array or Pandas DataSeries of radii
        :galaxy_vLum: A 1-D NumPy array or Pandas DataSeries of vLums
        :alpha: From the equation vObs^2 = vLum^2 + alpha*vLCM^2
        
        This avoids the issues inherent in calling sqrt, but can produce non-physical vObs^2
        
        Note that alpha is now defined differently, so that it is sqrt(old_alpha)"""
        
        vLCM = self.vLCM(galaxy_rad, galaxy_vLum)
        # this formula allows vNerosSquared to be negative for some fits which is problematic
        # chisquared becomes NaN
        return galaxy_vLum**2 + (alpha**2)*vLCM


    def vLCM(self, galaxy_rad, galaxy_vLum):
        """This computes the vLCM - the actual model
        
        The parameters are
        :galaxy_rad: A 1-D NumPy array or Pandas DataSeries of radii
        :galaxy_vLum: A 1-D NumPy array or Pandas DataSeries of vLums"""
        #phi_zero = 2.0e-04
        valid_rad = galaxy_rad <= self.mw_rad[-1]
        # Which MW_phi calculation should we use?
        # Right now, this way gives better chi_squared
        MW_phi = self.mw_phi_interp(galaxy_rad)
        #MW_phi = self.mw_phi
        trimmed_phi = MW_phi[:len(valid_rad)]
        phi_zero = trimmed_phi[-1]
        galaxy_phi = self.phi(galaxy_rad, galaxy_vLum)
        k = self.kappa(trimmed_phi, galaxy_phi, phi_zero)
        v1 = self.v1(trimmed_phi, galaxy_phi, phi_zero)
        v2 = self.v2(trimmed_phi, galaxy_phi, galaxy_vLum, phi_zero)
        vLCM = c * c * k * k * v1 * v2
        #vLCM = c * c * v2 * v2
        
        return vLCM


    def kappa(self, MW_phi, other_phi, phi_zero):
        """kappa(r) in the paper, just phi_gal(r)/phi_mw(r)"""
        
        return   other_phi  /  MW_phi 
        
        

    def v1(self, MW_phi, other_phi, phi_zero):
        etc = self._eTsiCurveMinusOne(MW_phi, other_phi, phi_zero)
        #sinh: uncomment following two lines
        num = (etc+1)**2 - 1
        den = 2*(1 + etc)
        return num/den
       

    def v2(self, MW_phi, other_phi, other_vlum, phi_zero):
        etFlat = self._eTsiFlatMinusOne(other_vlum)
        etCurve = self._eTsiCurveMinusOne(MW_phi, other_phi, phi_zero)
        #COSH:NN uncomment following two linesversion 1(e^2z =e^c * e^f)
        num = (etFlat +1)*(etCurve+1)  +1
        den = 2*np.sqrt((etFlat + 1)*(etCurve +1))
        return num/den
        


    def chiSquared(self, model, expected, error):
        """This computes chi squared"""
        chiSquared = 0

        for i in range(len(model)):
            chiSquared = chiSquared + (((model[i] - expected[i])**2) / (error[i]**2))
        return chiSquared / len(model)


    def _eTsiFlatMinusOne(self, other_vlum):
        """This computes eTsiFlat - 1, compared to the old code, for numerical stability"""
        
        beta = other_vlum / c
        numerator = 2*beta / (1 - beta)
        denominator = np.sqrt((1+beta) / (1-beta)) + 1
        return numerator / denominator


    def _eTsiCurveMinusOne(self, MW_phi, other_phi, phi_zero):
        """This computes eTsiCurve - 1, compared to the old code, for numerical stability"""
        #NOTE: this is the correct frame order and correct signs for numerical stability calculation
        numerator = ( 2*(MW_phi )-2*(other_phi )) / (1 - 2*(MW_phi ))
        denominator = np.sqrt((1 - 2*(other_phi)) / (1 - 2*(MW_phi ))) + 1
        return numerator / denominator