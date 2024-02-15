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

name = 'MW4_MHG0.25_GSE2_MHG0.5'
nameiso = 'MW4iso_fg0.2_MHG0.25_RC9'
lvl = 4

def get_cum_mass(name, lvl, rcut=20, numscalars=4):
    fname = 'CumMass_'+name+'-lvl'+str(lvl)+'.npy'
    
    file_path = Path(fname)
    if file_path.exists():
        print('CumMass already exists, exiting:', name, lvl)
        CumMass = np.load(fname)
        return CumMass
    
    gal = lowda.load_gal(name, lvl)
    num_idx = len(gal['idx_list'])
    CumMass = np.zeros((num_idx, numscalars+1), dtype=float)
    
    for i in range(num_idx):
        idx = gal['idx_list'][i]
        sn = gal['sn_idx'][idx]
        
        CumMass[i,0] = sn.Time.value
        
        MW_COM = gal['MW_COM'][idx]
        pos = sn.part0.pos.value - MW_COM
        r = np.linalg.norm(pos, axis=1)
        
        for sc in range(numscalars):
            mass = sn.part0.Masses.value * sn.part0.PassiveScalars[:,sc]
            CumMass[i,sc+1] = np.sum(mass[r < rcut])            

    np.save(fname, CumMass)
    
    del gal
    gc.collect()
            
    return CumMass

CumMass = get_cum_mass(name, lvl)
CumMass_iso = get_cum_mass(nameiso, lvl)

fig, ax = plt.subplots(1, 1, figsize=(columnwidth, 0.75*columnwidth))

ax.plot(CumMass[:,0], CumMass[:,2] - CumMass[0,2], c=tb_c[0], label='gaseous satellite (fiducial)')
ax.plot(CumMass_iso[:,0], CumMass_iso[:,2] - CumMass_iso[0,2], c=tb_c[1], label='isolated')

ax.set(xlabel=r'$t\,[\,\textrm{Gyr}\,]$')
ax.set(ylabel=r'$\textrm{accreted mass}\,[\,10^{10}\,M_{\odot}\,]$')

ax.set(xlim=(0, 8), ylim=(0, 1.0))

ax.legend(frameon=False)

fig.tight_layout()
fig.savefig('cumprof.pdf')
