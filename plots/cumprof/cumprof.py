import sys
import gc
from pathlib import Path
sys.path.append('/n/holylfs05/LABS/hernquist_lab/Users/abeane/GSEgas/plots')

import lowda

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# mpl.use('Agg')
# mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 8})
mpl.rc('text', usetex=True)
mpl.rc('text.latex', preamble=r'\usepackage{amsmath}')
# mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# color palette
tb_c = ['#4e79a7', '#f28e2b', '#e15759', '#76b7b2', '#59a14f',
        '#edc948', '#b07aa1', '#ff9da7', '#9c755f', '#bab0ac']

columnwidth = 242.26653 / 72.27 # converts pts to inches
textwidth = 513.11743 / 72.27
# frac of host CGM turned into stars

def get_n_T(sn):
    UnitLength = sn.parameters.UnitLength_in_cm
    UnitMass = sn.parameters.UnitMass_in_g
    UnitVelocity = sn.parameters.UnitVelocity_in_cm_per_s

    UnitTime = UnitLength / UnitVelocity
    UnitEnergy = UnitMass * UnitVelocity**2

    HYDROGEN_MASSFRAC = 0.76
    GAMMA = 5./3.
    PROTONMASS = 1.67262178e-24
    BOLTZMANN = 1.38065e-16

    InternalEnergy = sn.part0.InternalEnergy.value
    ElectronAbundance = sn.part0.ElectronAbundance
    Density = sn.part0.Density.value
    
    mu = 4 * PROTONMASS / (1 + 3 * HYDROGEN_MASSFRAC + 4 * HYDROGEN_MASSFRAC * ElectronAbundance)
    T = (GAMMA - 1.) * (InternalEnergy / BOLTZMANN) * (UnitEnergy / UnitMass) * mu

    n = Density / mu
    n *= UnitMass/UnitLength**3
    
    return n, T

def fig_tot_stellar_mass(gal, galiso, reset):
    
    def get_cum_star(gal, numscalars=4):
        CumMass = np.zeros((len(gal['idx_list']), numscalars+1), dtype=float)
        for i,idx in enumerate(gal['idx_list']):
            sn = gal['sn_idx'][idx]
            CumMass[i,0] = sn.Time.value
        
            if sn.NumPart_Total[4] == 0:
                continue
        
            is_star = sn.part4.GFM_StellarFormationTime > 0
            key = is_star
        
            for sc in range(numscalars):
                mass = sn.part4.Masses.value * sn.part4.PassiveScalars[:,sc]
                CumMass[i,sc+1] = np.sum(mass[key])

        return CumMass
    
    fname = 'data/CumStarMass.npy'
    fname_iso = 'data/CumStarMass_iso.npy'
    if reset:
        CumMass = get_cum_star(gal)
        CumMass_iso = get_cum_star(galiso)
        np.save(fname, CumMass)
        np.save(fname_iso, CumMass_iso)
    else:
        file_path = Path(fname)
        if not file_path.exists():
            return 1
        CumMass = np.load(fname)
        
        file_path = Path(fname_iso)
        if not file_path.exists():
            return 1
        CumMass_iso = np.load(fname_iso)
    
    fig, ax = plt.subplots(1, 1, figsize=(columnwidth, 0.75*columnwidth))
    
    ax.plot(CumMass[:,0], CumMass[:,2], c=tb_c[0], label=r'$\textrm{merger}$')
    ax.plot(CumMass_iso[:,0], CumMass_iso[:,2], c=tb_c[1], label=r'$\textrm{isolated}$')
    
    ax.set(xlim=(0, 8), xlabel=r'$t\,[\,\textrm{Gyr}\,]$')
    ax.set(ylim=(0, 0.5), ylabel=r'$\textrm{CGM mass in stars}\,[\,10^{10}\,M_{\odot}\,]$')
    
    ax.legend(frameon=False)
    
    fig.tight_layout()
    fig.savefig('CGM_in_stars.pdf')
    
    return 0

def fig_cold_mass(gal, galiso, reset):
    
    def get_cum_cold(gal, rmin=40, rmax=80, logTcut=4.5, numscalars=4):
        CumMass = np.zeros((len(gal['idx_list']), numscalars+1), dtype=float)
        for i,idx in enumerate(gal['idx_list']):
            sn = gal['sn_idx'][idx]
            CumMass[i,0] = sn.Time.value
        
            MW_COM = gal['MW_COM'][idx]
            pos = sn.part0.pos.value - MW_COM
            r = np.linalg.norm(pos, axis=1)
        
            _, T = get_n_T(sn)
        
            key = np.logical_and(r > rmin, r < rmax)
            key = np.logical_and(key, np.log10(T) < logTcut)
        
            for sc in range(numscalars):
                # mass = sn.part4.Masses.value * sn.part4.PassiveScalars[:,sc]
                mass = sn.part0.Masses.value * sn.part0.PassiveScalars[:,sc]
                CumMass[i,sc+1] = np.sum(mass[key])

        return CumMass
    
    fname = 'data/CumColdMass.npy'
    fname_iso = 'data/CumColdMass_iso.npy'
    if reset:
        CumMass = get_cum_cold(gal)
        CumMass_iso = get_cum_cold(galiso)
        np.save(fname, CumMass)
        np.save(fname_iso, CumMass_iso)
    else:
        file_path = Path(fname)
        if not file_path.exists():
            return 1
        CumMass = np.load(fname)
        
        file_path = Path(fname_iso)
        if not file_path.exists():
            return 1
        CumMass_iso = np.load(fname_iso)
    
    fig, ax = plt.subplots(1, 1, figsize=(columnwidth, 0.75*columnwidth))
    
    ax.plot(CumMass[:,0], CumMass[:,2], c=tb_c[0], label=r'$\textrm{merger}$')
    ax.plot(CumMass_iso[:,0], CumMass_iso[:,2], c=tb_c[1], label=r'$\textrm{isolated}$')
    
    ax.set(xlim=(0, 8), xlabel=r'$t\,[\,\textrm{Gyr}\,]$')
    ax.set(ylim=(0, 0.2), ylabel=r'$\textrm{cold CGM mass}\,[\,10^{10}\,M_{\odot}\,]$')
    
    ax.legend(frameon=False)
    
    fig.tight_layout()
    fig.savefig('CGM_in_cold.pdf')
    
    return 0
    
# CumMass = get_cum_mass(name, lvl)
# CumMass_iso = get_cum_mass(nameiso, lvl)

# fig, ax = plt.subplots(1, 1, figsize=(columnwidth, 0.75*columnwidth))

# ax.plot(CumMass[:,0], CumMass[:,2] - CumMass[0,2], c=tb_c[0], label='gaseous satellite (fiducial)')
# ax.plot(CumMass_iso[:,0], CumMass_iso[:,2] - CumMass_iso[0,2], c=tb_c[1], label='isolated')

# ax.set(xlabel=r'$t\,[\,\textrm{Gyr}\,]$')
# ax.set(ylabel=r'$\textrm{accreted mass}\,[\,10^{10}\,M_{\odot}\,]$')

# ax.set(xlim=(0, 8), ylim=(0, 1.0))

# ax.legend(frameon=False)

# fig.tight_layout()
# fig.savefig('cumprof.pdf')

if __name__ == '__main__':
    name = 'MW4_MHG0.25_GSE2_MHG0.5'
    nameiso = 'MW4iso_fg0.2_MHG0.25_RC9'
    lvl = 4

    reset = True if len(sys.argv) == 2 else False
    
    fields=['Coordinates', 'Masses', 'Velocities', 
            'InternalEnergy', 'Density', 'ElectronAbundance',
            'PassiveScalars', 'StarFormationRate', 'GFM_StellarFormationTime']
    
    if reset:
        gal = lowda.load_gal(name, lvl, parttype=[0, 2, 4], d_idx=20, fields=fields)
        galiso = lowda.load_gal(nameiso, lvl, parttype=[0, 2, 4], d_idx=20, fields=fields)
    else:
        gal = galiso = None
    
    if fig_tot_stellar_mass(gal, galiso, reset):
        print('failed on fig_tot_stellar_mass')
    else:
        print('success on fig_tot_stellar_mass')
    
    if fig_cold_mass(gal, galiso, reset):
        print('failed on fig_cold_mass')
    else:
        print('success on fig_cold_mass')
    