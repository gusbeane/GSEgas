import numpy as np
import h5py as h5
from pygsl import spline
import util
from warnings import warn
from tqdm import tqdm
from scipy.optimize import bisect
from scipy.interpolate import RectBivariateSpline

class TNGParameters(object):
    def __init__(self):
        # "numerical" parameters
        self.GFM_N_MASS_BINS = 200
        self.GFM_MIN_METAL = -20.

        # IMF parameters
        self.IMF_MinMass_Msun = 0.1
        self.IMF_MaxMass_Msun = 100.0
        self.IMF_A = 0.852464
        self.IMF_B = 0.237912
        self.IMF_sigma = 0.69
        self.IMF_mc = 0.079
        self.IMF_Mhinge = 1.0
        self.IMF_alpha = -2.3

        # SNII parameters
        self.SNII_MaxMass_Msun = 100.0
        self.SNII_MinMass_Msun = 8.0

        # SNIa parameters
        self.SNIa_Rate_Norm = 1.3e-3
        self.SNIa_Rate_TAU = 0.04
        self.GFM_SNIA_DTD_POWERLAW_INDEX = 1.12

class ChabrierIMF(object):
    def __init__(self, param, renormalize=False):
        self._verify_params(param)
        self._init_imf(renormalize)
        self._init_cumulative_imf()

    def _verify_params(self, param):
        self._required_params = ['IMF_MinMass_Msun', 'IMF_MaxMass_Msun', 'GFM_N_MASS_BINS', 'IMF_A', 'IMF_B', 'IMF_sigma', 'IMF_mc', 'IMF_Mhinge', 'IMF_alpha']
        for p in self._required_params:
            if not hasattr(param, p):
                raise AttributeError("Must specify parameter " + p)
        self.p = param
    
    def _init_imf(self, renormalize):
        imf_dlog10_Msun = (np.log10(self.p.IMF_MaxMass_Msun) - np.log10(self.p.IMF_MinMass_Msun)) / (self.p.GFM_N_MASS_BINS - 1)
        imf_mass_bin_log10 = np.log10(self.p.IMF_MinMass_Msun) + np.arange(self.p.GFM_N_MASS_BINS) * imf_dlog10_Msun
        imf_mass_bin       = np.power(10, imf_mass_bin_log10)

        imf_number = np.zeros_like(imf_mass_bin)

        imf_number_upper = self.p.IMF_B * np.power(imf_mass_bin, self.p.IMF_alpha)

        tmpa = np.power((np.log10(imf_mass_bin) - np.log10(self.p.IMF_mc)), 2.0)
        tmpb = -2.0 * pow(self.p.IMF_sigma, 2.0)
        imf_number_lower = self.p.IMF_A * np.exp(tmpa / tmpb) / imf_mass_bin

        imf_number[imf_mass_bin > self.p.IMF_Mhinge] = imf_number_upper[imf_mass_bin > self.p.IMF_Mhinge]
        imf_number[imf_mass_bin <= self.p.IMF_Mhinge] = imf_number_lower[imf_mass_bin <= self.p.IMF_Mhinge]

        self._imf_dlog10_Msun = imf_dlog10_Msun
        self._imf_mass_bin_log10 = imf_mass_bin_log10
        self._imf_mass_bin = imf_mass_bin
        self._imf_number = imf_number

        # now renormalize. if the norm is not within 1% of unity, user must specify that they would like to renormalize.
        xmin = np.log10(self.p.IMF_MinMass_Msun)
        xmax = np.log10(self.p.IMF_MaxMass_Msun)
        norm = self.integrate_imf_by_mass(xmin, xmax)

        if np.abs(1.0 - norm) > 0.01 and not renormalize:
            err_str = 'IMF normalization is not within 1% of unity: ' + str(norm) + "\n"
            err_str += "You must specify renormalize=True to renormalize, but make sure this is what you want"
            raise ValueError(err_str)
        else:
            self._imf_number /= norm
            self.renormalized = renormalize
    
    def _init_cumulative_imf(self):
        xmin = np.repeat(self._imf_mass_bin_log10[0], self.p.GFM_N_MASS_BINS)
        xmax = self._imf_mass_bin_log10
        self._cumulative_imf = self.integrate_imf_by_mass(xmin, xmax)

    def _get_imf_bins(self, log_min_mass, log_max_mass):
        log_min_mass = np.clip(log_min_mass, self._imf_mass_bin_log10[0], self._imf_mass_bin_log10[-1])
        log_max_mass = np.clip(log_max_mass, self._imf_mass_bin_log10[0], self._imf_mass_bin_log10[-1])

        ilow  = (log_min_mass[:,None] >= self._imf_mass_bin_log10).sum(axis=1) - 1
        ihigh = 1 + (log_max_mass[:,None] > self._imf_mass_bin_log10).sum(axis=1) - 1
        return ilow, ihigh

    def _integrate_imf(self, log_min_mass, log_max_mass, integrand):
        result = util.trap_integrate(self._imf_mass_bin_log10, integrand, log_min_mass, log_max_mass, self._imf_dlog10_Msun)
        result *= np.log(10.0) # convert from ln to log10
        return result

    @util.auto_convert_to_ndarray
    def integrate_imf_by_number(self, log_min_mass, log_max_mass):
        integrand = self._imf_number * self._imf_mass_bin
        return self._integrate_imf(log_min_mass, log_max_mass, integrand)
    
    @util.auto_convert_to_ndarray
    def integrate_imf_by_mass(self, log_min_mass, log_max_mass):
        integrand = self._imf_number * self._imf_mass_bin**2
        # print(integrand)
        return self._integrate_imf(log_min_mass, log_max_mass, integrand)
    
    @util.auto_convert_to_ndarray
    def integrate_imf_by_func(self, log_min_mass, log_max_mass, func):
        integrand = self._imf_number * self._imf_mass_bin * func
        return self._integrate_imf(log_min_mass, log_max_mass, integrand)

    def __call__(self, x):
        # linearly interpolate the imf_number
        ilow, ihigh = self._get_imf_bins(np.log10(x), np.log10(x))
        s = (np.log10(x) - self._imf_mass_bin_log10[ilow]) / self._imf_dlog10_Msun
        return (1-s) * self._imf_number[ilow] + s * self._imf_number[ihigh]

class LifeTimes(object):
    def __init__(self, yieldpath, param):
        data = h5.File(yieldpath + '/LifeTimes.hdf5', mode='r')

        self._verify_params(param)
        self._read_N(data)
        self._read_tables(data)

        data.close()

    def _verify_params(self, param):
        self._required_params = ['IMF_MaxMass_Msun']
        for p in self._required_params:
            if not hasattr(param, p):
                raise AttributeError("Must specify parameter " + p)
        self.p = param

    def _read_N(self, data):
        self.N_MASS = int(data['Number_of_masses'][()])
        self.N_Z = int(data['Number_of_metallicities'][()])
    
    def _read_tables(self, data):
        self.Mass = data['Masses'][:]
        self.Metallicity = data['Metallicities'][:]
        self.Dyingtime = np.log10(data['Lifetimes'][:])
    
    @util.auto_convert_to_ndarray
    def geriatric_mass(self, tau_in_Gyr, z):
        logtau = np.log10(tau_in_Gyr * 1e9)
        
        i_metal = (z[:,None] >= self.Metallicity).sum(axis=1) - 1
        i_metal = np.clip(i_metal, 0, self.Metallicity.shape[0]-2)

        z_i = self.Metallicity[i_metal]
        z_ip1 = self.Metallicity[i_metal+1]

        s_metal = (z - z_i) / (z_ip1 - z_i)
        s_metal = np.clip(s_metal, 0, 1)

        # interpolate in metal
        aux_time = (1.0-s_metal[:,None]) * self.Dyingtime[i_metal] + s_metal[:,None] * self.Dyingtime[i_metal+1]

        i_time = (logtau[:,None] <= aux_time).sum(axis=1) - 1
        i_time = np.clip(i_time, 0, aux_time.shape[1]-2)

        auxtime_i = aux_time[np.arange(len(i_time)),i_time]
        auxtime_ip1 = aux_time[np.arange(len(i_time)),i_time+1]

        s_time = (logtau - auxtime_i) / (auxtime_ip1 - auxtime_i)
        s_time = np.clip(s_time, 0, 1)

        mass = (1-s_time) * self.Mass[i_time] + s_time * self.Mass[i_time+1]
        mass = np.clip(mass, self.p.IMF_MinMass_Msun, self.p.IMF_MaxMass_Msun)

        return mass

class yieldSNII_AGB(object):
    def __init__(self, yieldpath, param, IMF, SNII=False, AGB=False):
        if SNII ^ AGB:
            self.SNII = SNII
            self.AGB = AGB
        else:
            raise ValueError("Must specify either SNII or AGB (but not both)")
        
        self._verify_params(param)

        if self.SNII:
            self.fname = yieldpath + '/SNII.hdf5'
            data = h5.File(self.fname, mode='r')
        else:
            self.fname = yieldpath + '/AGB.hdf5'
            data = h5.File(self.fname, mode='r')

        self._read_N(data)
        self._read_tables(data)

        self._process_tables()
        self._init_splines()

        data.close()

        self.IMF = IMF

    def _verify_params(self, param):
        self._required_params = ['IMF_MinMass_Msun', 'IMF_MaxMass_Msun', 'GFM_N_MASS_BINS', 'GFM_MIN_METAL',
                                 'SNII_MaxMass_Msun', 'SNII_MinMass_Msun']
        for p in self._required_params:
            if not hasattr(param, p):
                raise AttributeError("Must specify parameter " + p)
        self.p = param

        if self.SNII:
            self._MaxMass_Msun = self.p.SNII_MaxMass_Msun
            self._MinMass_Msun = self.p.SNII_MinMass_Msun
        else:
            self._MaxMass_Msun = self.p.SNII_MinMass_Msun
            self._MinMass_Msun = self.p.IMF_MinMass_Msun

        self._logMinMass = np.log10(self._MinMass_Msun)
        self._logMaxMass = np.log10(self._MaxMass_Msun)

    def _read_N(self, data):
        self.N_ELEMENTS = int(data['Number_of_species'][()])
        self.N_MASS = int(data['Number_of_masses'][()])
        self.N_Z = int(data['Number_of_metallicities'][()])
          
    def _read_tables(self, data):
        self.Mass = data['Masses'][:]
        self.ElementName = data['Species_names'][:]
        self.Metallicity = data['Metallicities'][:]

        self.Yield_names = data['Yield_names'][:]
        self.Ejecta = []
        self.TotalMetals = []
        self.Yield = []
        for i in range(self.N_Z):
            name = self.Yield_names[i].decode()
            self.Ejecta.append(data['Yields/' + name + '/Ejected_mass'][:])
            self.TotalMetals.append(data['Yields/' + name + '/Total_Metals'][:])
            self.Yield.append(data['Yields/' + name + '/Yield'][:])
        self.Ejecta = np.array(self.Ejecta)
        self.TotalMetals = np.array(self.TotalMetals)
        self.Yield = np.array(self.Yield)
    
    def _process_tables(self):
        self.Mass = np.log10(self.Mass)
        self.Metallicity[self.Metallicity < 10.**self.p.GFM_MIN_METAL] = 10.**self.p.GFM_MIN_METAL
        self.Metallicity = np.log10(self.Metallicity)

    def _init_splines(self):
        lm_min = np.log10(self.p.IMF_MinMass_Msun)
        lm_max = np.log10(self.p.IMF_MaxMass_Msun)

        dlm = (lm_max - lm_min) / (self.p.GFM_N_MASS_BINS - 1)

        self.yield_mass_bin = dlm * np.arange(self.p.GFM_N_MASS_BINS) + lm_min

        interp_ej = spline.linear(int(self.N_MASS))
        interp_met = spline.linear(int(self.N_MASS))

        self.Ejecta_spline = np.zeros((self.N_Z, self.p.GFM_N_MASS_BINS))
        self.TotalMetals_spline = np.zeros((self.N_Z, self.p.GFM_N_MASS_BINS))
        for i in range(self.N_Z):
            yield_ej_i = self.Ejecta[i] / (np.power(10., self.Mass))
            yield_met_i = self.TotalMetals[i] / (np.power(10., self.Mass))

            mass = self.Mass
            interp_ej.init(mass, yield_ej_i)
            interp_met.init(mass, yield_met_i)

            for k in range(self.p.GFM_N_MASS_BINS):
                if self.yield_mass_bin[k] < mass[0]:
                    result_ej = yield_ej_i[0]
                    result_met = yield_met_i[0]
                elif self.yield_mass_bin[k] > mass[-1]:
                    result_ej = yield_ej_i[-1]
                    result_met = yield_met_i[-1]
                else:
                    result_ej = interp_ej.eval(self.yield_mass_bin[k])
                    result_met = interp_met.eval(self.yield_mass_bin[k])
            
                self.Ejecta_spline[i, k] = result_ej * np.power(10., self.yield_mass_bin[k])
                self.TotalMetals_spline[i, k] = result_met * np.power(10., self.yield_mass_bin[k])

    def _spline_interp(self, log_metallicity, spline):
        i_metal = (log_metallicity[:,None] >= self.Metallicity).sum(axis=1) - 1
        i_metal = np.clip(i_metal, 0, self.Metallicity.shape[0]-2)

        logz_i = self.Metallicity[i_metal]
        logz_ip1 = self.Metallicity[i_metal+1]

        s_metal = (log_metallicity - logz_i) / (logz_ip1 - logz_i)
        s_metal = np.clip(s_metal, 0, 1)

        # interpolate in metal
        spline = (1.0-s_metal[:,None]) * spline[i_metal] + s_metal[:,None] * spline[i_metal+1]

        return spline

    @util.auto_convert_to_ndarray
    def get_total_mass_ejected(self, log_min_mass, log_max_mass, log_metallicity):
        # interpolate the spline in metallicity first
        mej_spline = self._spline_interp(log_metallicity, self.Ejecta_spline)
        met_spline = self._spline_interp(log_metallicity, self.TotalMetals_spline)

        metallicity = 10.**log_metallicity

        log_min_mass = np.clip(log_min_mass, self._logMinMass, self._logMaxMass)
        log_max_mass = np.clip(log_max_mass, self._logMinMass, self._logMaxMass)

        mej = self.IMF.integrate_imf_by_func(log_min_mass, log_max_mass, mej_spline)
        met = self.IMF.integrate_imf_by_func(log_min_mass, log_max_mass, metallicity[:,np.newaxis]*mej_spline + met_spline)

        return mej, met


class yieldSNIa(object):
    def __init__(self, yieldpath, param):
        self._verify_params(param)
        self.fname = yieldpath + '/SNIa.hdf5'
        data = h5.File(self.fname, mode='r')

        self.TotalMetals_spline = data['Total_Metals'][()]

        data.close()
    
    def _verify_params(self, param):
        self._required_params = ['SNIa_Rate_Norm', 'SNIa_Rate_TAU', 'GFM_SNIA_DTD_POWERLAW_INDEX']
        for p in self._required_params:
            if not hasattr(param, p):
                raise AttributeError("Must specify parameter " + p)
        self.p = param

        # renormalize normalization to get correct # of SNIa in age of universe
        # this needs to be changed for comoving runs
        # whyyyyyyyy would you do it this way
        AgeOfUniverse_in_Gyr = 13.7
        fac = (1 - (AgeOfUniverse_in_Gyr/self.p.SNIa_Rate_TAU)**(1-self.p.GFM_SNIA_DTD_POWERLAW_INDEX))
        self._SNIa_Rate_Norm_Norm = self.p.SNIa_Rate_Norm / fac

    def number_of_SNIa(self, tau_in_Gyr, dtau_in_Gyr):
        tau1 = tau_in_Gyr + dtau_in_Gyr
        tau_in_Gyr = np.clip(tau_in_Gyr, self.p.SNIa_Rate_TAU, np.inf)
        tau1 = np.clip(tau1, self.p.SNIa_Rate_TAU, np.inf)
        dtau_in_Gyr = tau1 - tau_in_Gyr

        # print(tau_in_Gyr, dtau_in_Gyr)
        # dtau_in_Gyr[tau_in_Gyr < self.p.SNIa_Rate_TAU] = tau_in_Gyr + dtau_in_Gyr - self.p.SNIa_Rate_TAU
        # dtau_in_Gyr[dtau_in_Gyr < 0.0] = 0.0
        # tau_in_Gyr[tau_in_Gyr < self.p.SNIa_Rate_TAU] = self.p.SNIa_Rate_TAU

        # print(tau_in_Gyr, dtau_in_Gyr)

        integral = (tau_in_Gyr/self.p.SNIa_Rate_TAU)**(1-self.p.GFM_SNIA_DTD_POWERLAW_INDEX)
        integral -= ((tau_in_Gyr+dtau_in_Gyr)/self.p.SNIa_Rate_TAU)**(1-self.p.GFM_SNIA_DTD_POWERLAW_INDEX)
        
        integral[tau_in_Gyr + dtau_in_Gyr <= self.p.SNIa_Rate_TAU] = 0.0
        assert np.all(integral >= 0.0)

        return self._SNIa_Rate_Norm_Norm * integral

    def get_total_mass_ejected(self, tau_in_Gyr, dtau_in_Gyr):
        mej = self.number_of_SNIa(tau_in_Gyr, dtau_in_Gyr) * self.TotalMetals_spline
        met = mej
        return mej, met

class GFM(object):
    def __init__(self, yieldpath, param, IMF=None, LT=None, SNII=None, AGB=None, SNIa=None):
        self.IMF = ChabrierIMF(param) if IMF is None else IMF
        self.LT = LifeTimes(yieldpath, param) if LT is None else LT
        self.SNII = yieldSNII_AGB(yieldpath, param, self.IMF, SNII=True) if SNII is None else SNII
        self.AGB = yieldSNII_AGB(yieldpath, param, self.IMF, AGB=True) if AGB is None else AGB
        self.SNIa = yieldSNIa(yieldpath, param) if SNIa is None else SNIa

        self.p = param
        self._init_cumulative_eject_fraction()
    
    @util.auto_convert_to_ndarray
    def get_total_mass_ejected(self, tau, dtau, z, verbose=True):
        # for numerical stability
        tau = np.clip(tau, 1e-12, None)

        log_min_mass = np.log10(self.LT.geriatric_mass(tau+dtau, z))
        log_max_mass = np.log10(self.LT.geriatric_mass(tau, z))
        log_metallicity = np.log10(z)

        mej_SNII, met_SNII = self.SNII.get_total_mass_ejected(log_min_mass, log_max_mass, log_metallicity)
        mej_AGB, met_AGB = self.AGB.get_total_mass_ejected(log_min_mass, log_max_mass, log_metallicity)
        mej_SNIa, met_SNIa = self.SNIa.get_total_mass_ejected(tau, dtau)

        mej = mej_SNII + mej_AGB + mej_SNIa
        met = met_SNII + met_AGB + met_SNIa

        return mej, met
    
    def _frac_time_min(self, tau, z, frac):
        cum_ej, _ = self.get_total_mass_ejected(1e-12, tau, z, verbose=False)
        # print(cum_ej)
        diff = cum_ej - frac
        # print(tau, frac, cum_ej, diff)
        # print(diff)
        return diff

    @util.auto_convert_to_ndarray
    def get_time_of_frac(self, frac, z):
        N = frac.shape[0]
        tau = np.zeros_like(frac)
        # for i in tqdm(range(N)):
        for i in range(N):
            try:
                tau[i] = bisect(self._frac_time_min, 0.002, 14.0, args=(z[i], frac[i]))
            except:
                tau[i] = 14.0
                # print(frac[i], z[i])
                # assert False

        return tau

    @util.auto_convert_to_ndarray
    def get_total_mass_ejected_at_frac(self, frac0, frac1, z):
        # frac0 < frac1
        tau0 = self.get_time_of_frac(frac0, z)
        tau1 = self.get_time_of_frac(frac1, z)
        # print(tau0, tau1-tau0)
        return self.get_total_mass_ejected(tau0, tau1-tau0, z)

    def get_tracer_metallicity(self, frac, z):
        mej = self._mej_interp(z, frac, grid=False)
        met = self._met_interp(z, frac, grid=False)
        return met / mej

    def _init_cumulative_eject_fraction(self):
        frac_bins = np.linspace(0.0, 0.5, 100)
        z_bins = np.linspace(10.**(self.p.GFM_MIN_METAL), 0.04, 16)

        z_mesh, frac_mesh = np.meshgrid(z_bins, frac_bins, indexing='ij')
        z_flat = z_mesh.flatten()
        frac_flat = frac_mesh.flatten()

        nonzero = frac_flat > 0.0

        mej = np.zeros_like(frac_flat)
        met = np.zeros_like(frac_flat)

        tau_nonzero = self.get_time_of_frac(frac_flat[nonzero], z_flat[nonzero])
        mej_nonzero, met_nonzero = self.get_total_mass_ejected(tau_nonzero-1e-6, 2e-6, z_flat[nonzero])

        mej[nonzero] = mej_nonzero
        met[nonzero] = met_nonzero

        mej = mej.reshape((len(z_bins), len(frac_bins)))
        met = met.reshape((len(z_bins), len(frac_bins)))

        # numerical stability
        assert np.all(mej >= -1e-12)
        mej[mej < 0.0] = 0.0

        self._mej = mej
        self._met = met

        self._frac_bins = frac_bins
        self._z_bins = z_bins

        self._mej_interp = RectBivariateSpline(z_bins, frac_bins, mej)
        self._met_interp = RectBivariateSpline(z_bins, frac_bins, met)
