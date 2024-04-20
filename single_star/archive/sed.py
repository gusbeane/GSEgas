import numpy as np
import platform
import h5py as h5
from scipy.interpolate import RegularGridInterpolator, interp1d
from numba import njit

def auto_convert_to_ndarray(func):
    def wrapper(*args, **kwargs):
        # ignore first arg if we are a class method
        arg0 = args[1] if isinstance(args[0], object) else args[0]
        args_in = [np.array([arg]) if isinstance(arg, (float, int)) 
                   else np.asarray(arg) if isinstance(arg, list)
                   else arg for arg in args]
        result = func(*args_in, **kwargs)
        # Convert ndarray result to float if the input was a float
        # print(args[0])
        return float(result[0]) if isinstance(arg0, (float, int)) else result
    return wrapper

class IMFClass(object):
    def __init__(self, mmin=0.1, mmax=100, A=0.852464, B=0.237912,
                 mc=0.079, sigma=0.69, mhinge=1):
        
        self.mmin = mmin
        self.mmax = mmax
        self.A = A
        self.B = B
        self.mc = mc
        self.sigma = sigma
        self.mhinge = mhinge
    
        self._test_norm()

    def _test_norm(self):
        m = np.logspace(np.log10(self.mmin), np.log10(self.mmax), 1000)
        imf = self(m)
        norm = np.trapz(m*imf, m)
        assert np.isclose(norm, 1.0)

    @auto_convert_to_ndarray
    def __call__(self, m):
        imf = np.zeros(m.shape)

        m = np.clip(m, self.mmin, self.mmax)

        imf_lower = self.A * m**(-1.) * np.exp(-np.log10(m/self.mc)**2 / (2*self.sigma**2))
        imf_upper = self.B * m**(-2.3)

        imf[m < self.mhinge] = imf_lower[m < self.mhinge]
        imf[m >= self.mhinge] = imf_upper[m >= self.mhinge]

        imf[m < self.mmin] = 0.0
        imf[m > self.mmax] = 0.0

        return imf

@auto_convert_to_ndarray
def _my_interp(x, xp, yp, deriv=False):
    # x is assumed to be an array of length (N,)
    # xp is assumed to be an array of length (N, M)
    # yp is assumed to be an array of length (M,)

    # if x is a scalar, make it an array of length (1,)

    # note, this assumes that a _different_ xp is used for each element of x
    # if you have a single xp for many different x, you should use np.interp
    # also, this assumes that each row of xp is sorted in ascending order
    # and will return yp[0] if x < xp[0] and yp[-1] if x > xp[-1]

    # if deriv=True, then we instead return the derivative of the interpolation
    # this is equivalent to returning (yp[i+1] - yp[i]) / (xp[i+1] - xp[i]) if
    # xp[i] < x < xp[i+1]

    # return_float = False
    # if isinstance(x, (int, float)):
        # x = np.array([x])
        # xp = np.expand_dims(xp, 0)
        # return_float = True
    
    if xp.ndim == 1:
        xp = np.expand_dims(xp, 0)
    
    # print(np.shape(x), np.shape(xp))

    bins = (x[:,None] >= xp).sum(axis=1) - 1
    bins = np.clip(bins, 0, xp.shape[1]-2)

    xp_i = xp[np.arange(bins.shape[0]),bins]
    xp_ip1 = xp[np.arange(bins.shape[0]),bins+1]

    s = (x - xp_i) / (xp_ip1 - xp_i)

    # if we return the derivative, we can end here, but need to set to zero if beyond bounds
    if deriv:
        y = (yp[bins+1] - yp[bins]) / (xp_ip1 - xp_i)
        # y[np.logical_or(s < 0, s > 1)] = 0.0
        # y[s < 0] = y
        return y

    s = np.clip(s, 0, 1)
    return (1-s)*yp[bins] + s*yp[bins+1]   

class LifeTimesClass(object):
    def __init__(self, fname):
        self.fname = fname
        data = h5.File(self.fname)
        self._M = data['Masses'][:]
        self._Z = data['Metallicities'][:]
        self._Tau = data['Lifetimes'][:]
        data.close()
        
        self._init_interp()
        self._init_inverse_interp()
    
    @auto_convert_to_ndarray
    def __call__(self, z, m):
        z = np.clip(z, self._Zmin, self._Zmax)
        m = np.clip(m, self._Mmin, self._Mmax)
        return self._interp((z, m))

    def _init_interp(self):
        self._Zmin = self._Z[0]
        self._Zmax = self._Z[-1]
        self._Mmin = self._M[0]
        self._Mmax = self._M[-1]

        # TODO: use two consecutive np.interp to reduce reliance on scipy

        self._interp = RegularGridInterpolator((self._Z, self._M), self._Tau, 
                                               bounds_error=True)
    
    def _init_inverse_interp(self):
        self._Tauz_interp = interp1d(self._Z, self._Tau, axis=0,
                                     bounds_error=False, fill_value=(self._Tau[0], self._Tau[-1]))

    @auto_convert_to_ndarray
    def geriatric_mass(self, z, tau, deriv=False):
        # first interpolate in Z to get an array of lifetimes as a function of mass
        Tau_z = self._Tauz_interp(z)

        # need to flip arrays for _my_interp to work, it assumes x is sorted ascending
        x = np.flip(Tau_z)
        y = np.flip(self._M)

        return _my_interp(tau, x, y, deriv=deriv)

class AGB_SNIIClass(object):
    def __init__(self, fname):
        self.fname = fname
        data = h5.File(fname)
        self._Z = data['Metallicities'][:]
        self._M = data['Masses'][:]

        self._Mej = []
        self._Zej = []
        for i in range(data['Number_of_metallicities'][()]):
            z_str = data['Yield_names'][i]
            self._Mej.append(data['Yields'][z_str]['Ejected_mass'][:])
            self._Zej.append(data['Yields'][z_str]['Total_Metals'][:])
        self._Mej = np.array(self._Mej)
        self._Zej = np.array(self._Zej)

        data.close()
    
        self._init_interp()

    def _init_interp(self):
        self._Zmin = self._Z[0]
        self._Zmax = self._Z[-1]
        self._Mmin = self._M[0]
        self._Mmax = self._M[-1]

        self._Mej_interp = RegularGridInterpolator((self._Z, self._M), self._Mej, 
                                                    bounds_error=True)
        self._Zej_interp = RegularGridInterpolator((self._Z, self._M), self._Zej, 
                                                    bounds_error=True)
    
    @auto_convert_to_ndarray
    def fej(self, z, m):
        z = np.clip(z, self._Zmin, self._Zmax)
        m = np.clip(m, self._Mmin, self._Mmax)
        return self._Mej_interp((z, m)) / m
    
    @auto_convert_to_ndarray
    def fZej(self, z, m):
        z = np.clip(z, self._Zmin, self._Zmax)
        m = np.clip(m, self._Mmin, self._Mmax)
        return self._Zej_interp((z, m)) / m

class SNIaClass(object):
    def __init__(self, fname, tau8=4E7, s=1.12, N0=1.3E-3):
        self.fname = fname

        data = h5.File(fname)
        self._Mej = data['Total_Metals'][()]
        self._MZej = data['Total_Metals'][()]
        data.close()

        self.tau8 = tau8
        self.s = s
        self.N0 = N0

    @auto_convert_to_ndarray
    def DTD(self, tau):
        dtd = self.N0 * (tau/self.tau8)**(-self.s) * (self.s-1)/self.tau8
        dtd[tau < self.tau8] = 0.0

        return dtd
    
    def fej(self):
        return self._Mej
    
    def fZej(self):
        return self._MZej

class SEDClass(object):
    def __init__(self, yieldpath):
        self.yieldpath = yieldpath
    
        self.IMF = IMFClass()
        self.LifeTimes = LifeTimesClass(yieldpath + '/LifeTimes.hdf5')
        self.SNII = AGB_SNIIClass(yieldpath + '/SNII.hdf5')
        self.AGB = AGB_SNIIClass(yieldpath + '/AGB.hdf5')
        self.SNIa = SNIaClass(yieldpath + '/SNIa.hdf5')

        self._init_param()
     
    def _init_param(self):
        # TODO: make these user changable
        self.IMF_MinMass_Msun = 0.1
        self.IMF_MaxMass_Msun = 100.0
        self.AGB_MassTransferOn = 1
        self.SNIa_MassTransferOn = 1
        self.SNII_MassTransferOn = 1
        self.SNII_MinMass_Msun = 8.0
        self.SNII_MaxMass_Msun = 100.0
        self.SNIa_Rate_TAU = 0.04
        self.SNIa_Rate_Norm = 1.3e-3

    @auto_convert_to_ndarray
    def ejected_metallicity(self, z, tau):
        mu = self.LifeTimes.geriatric_mass(z, tau)
        dmu_dtau = self.LifeTimes.geriatric_mass(z, tau, deriv=True)
        imf = self.IMF(mu)

        # print(imf)

        # get the mass ejected by AGB and SNII
        fej_AGB = self.AGB.fej(z, mu)
        fZej_AGB = self.AGB.fZej(z, mu)
        fej_SNII = self.SNII.fej(z, mu)
        fZej_SNII = self.SNII.fZej(z, mu)

        # now filter by min SNII mass
        fej_AGB[mu >= self.SNII_MinMass_Msun] = 0.0
        fZej_AGB[mu >= self.SNII_MinMass_Msun] = 0.0
        fej_SNII[mu < self.SNII_MinMass_Msun] = 0.0
        fZej_SNII[mu < self.SNII_MinMass_Msun] = 0.0
    
        # now multiply by mu * IMF and derivative of inverse lifetime function to get time derivative
        factor = - mu * imf * dmu_dtau
        dmej_dtau_AGB = fej_AGB * factor
        dmzej_dtau_AGB = fZej_AGB * factor
        dmej_dtau_SNII = fej_SNII * factor
        dmzej_dtau_SNII = fZej_SNII * factor

        # get the mass ejected by a single SNIa, then multiply by DTD to get time derivative
        mej_SNIa = self.SNIa.fej()
        mzej_SNIa = self.SNIa.fZej()
        dmej_dtau_SNIa = mej_SNIa * self.SNIa.DTD(tau)
        dmzej_dtau_SNIa = mzej_SNIa * self.SNIa.DTD(tau)

        # total mass ejected and total metallicity ejected
        dmej_dtau = dmej_dtau_SNII + dmej_dtau_AGB + dmej_dtau_SNIa
        dmzej_dtau = dmzej_dtau_SNII + dmzej_dtau_AGB + dmzej_dtau_SNIa

        # also include the initial metallicity of the star for SNII and SNIa
        dmzej_dtau += z * (dmej_dtau_SNII + dmej_dtau_AGB)
        zej = np.divide(dmzej_dtau, dmej_dtau, where=dmzej_dtau > 0.0)

        # if time = 0, then this means tail is big bang. so we just eject the tail metallicity
        zej[tau == 0.0] = z[tau == 0.0]
    
        if np.any(zej < z):            
            key = zej < z
            print('Warning: ejected metallicity is less than initial metallicity')
            print('z = ', z[key])
            print('tau = ', tau[key])
            print('zej = ', zej[key])
            print('imf = ', imf[key])
            print('mu = ', mu[key])
            print('dmu_dtau = ', dmu_dtau[key])
            print('fej_AGB = ', fej_AGB[key])
            print('fej_SNII = ', fej_SNII[key])
            print('fzej_AGB = ', fZej_AGB[key])
            print('fzej_SNII = ', fZej_SNII[key])
   

        return zej