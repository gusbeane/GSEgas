import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

columnwidth = 242.26653 / 72.27
textwidth = 513.11743 / 72.27

mpl.use('Agg')
mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 8})
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

tb_c = ['#4e79a7', '#f28e2b', '#e15759', '#76b7b2', '#59a14f',
        '#edc948', '#b07aa1', '#ff9da7', '#9c755f', '#bab0ac']

basepath = '/n/holylfs05/LABS/hernquist_lab/Users/abeane/GSEgas/'

dil_path = basepath + 'anlys/dilution/'

MW3iso_corona3 = 'MW3iso_fg0.7_MHG0.25_RC9'
MW3_GSE2_merge2 = 'MW3_MHG0.25_GSE2_MHG0.18_Rcut10'

lvl = 4

def run(snap=800):
    fname = dil_path+'dil_'+MW3_GSE2_merge2+'-lvl'+str(lvl)+'.npy'
    dil = np.load(fname, allow_pickle=True).item()
    
    aveR = dil['aveR'][snap]
    dilution = dil['dilution'][snap]
    metallicity = dil['metallicity'][snap]
    
    fig, ax = plt.subplots(2, 1, figsize=(columnwidth, 1.2*columnwidth), sharex=True)
    
    ax[0].plot(aveR, dilution, c=tb_c[0])
    ax[0].set(ylim=(0, 0.2))
    ax[0].set_ylabel(r'dilution fraction')
    
    ax[1].plot(aveR, np.log10(metallicity/0.0127), c=tb_c[0])
    ax[1].axhline(-0.3, c='k', ls='dashed')
    ax[1].axhline(-0.6, c='k', ls='dashed')
    ax[1].set(xlim=(0, 15), ylim=(-0.7, -0.2))
    ax[1].set_xlabel(r'$R\,[\,\text{kpc}\,]$')
    ax[1].set_ylabel(r'$[Z/H]$')
    
    fig.tight_layout()
    fig.savefig('dil.pdf')

if __name__ == '__main__':
    run()

