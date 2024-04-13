import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import arepo

import os

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

def make_data(subs):
    dat = {}
    dat['Time'] = []
    dat['Mass'] = []
    dat['Size'] = []
    
    i = 0
    while True:
        if i in subs.keys():
            sub = subs[i]
            dat['Time'].append(sub.Time)
            dat['Mass'].append(sub.SubhaloMassInRadType[0,4])
            dat['Size'].append(sub.SubhaloHalfmassRadType[0,4])
            i += 1
        else:
            break
    
    for k in dat.keys():
        dat[k] = np.array(dat[k])
        
    return dat

def make_fig(subsMW, subsGSE, reset):
    
    fnameMW = 'datMW.npy'
    fnameGSE = 'datGSE.npy'
    if reset:
        datMW = make_data(subsMW)
        datGSE = make_data(subsGSE)
        np.save(fnameMW, datMW)
        np.save(fnameGSE, datGSE)
    else:
        if not os.path.isfile(fnameMW):
            return 1
        datMW = np.load(fnameMW, allow_pickle=True).tolist()
        
        if not os.path.isfile(fnameGSE):
            return 1
        datGSE = np.load(fnameGSE, allow_pickle=True).tolist()
    
    fig, ax = plt.subplots(2, 1, figsize=(columnwidth, 1.25*columnwidth), sharex=True)
    
    ax[0].plot(datMW['Time'], datMW['Mass'], c=tb_c[0], label=r'$\textrm{MW}$')
    ax[0].plot(datGSE['Time'], datGSE['Mass'], c=tb_c[1], label=r'$\textrm{GSE}$')
    
    ax[0].set(xlim=(0, 8))
    ax[0].set(ylim=(0.01, 2), ylabel=r'$\textrm{stellar mass}\,[\,10^{10}\,M_{\odot}\,]$')
    ax[0].set(yscale='log')
    
    ax[0].axhline(0.6, c=tb_c[0], ls='dashed')
    ax[0].axhline(0.05, c=tb_c[1], ls='dashed')
    ax[0].axvline(3, c='k', ls='dashed', alpha=0.5)
    
    ax[0].legend(frameon=False)
    ax[0].xaxis.set_tick_params(labelbottom=True)
    
    ax[1].plot(datMW['Time'], datMW['Size'], c=tb_c[0])
    ax[1].plot(datGSE['Time'], datGSE['Size'], c=tb_c[1])
    
    # ax[1].axhline(2*1.68, c=tb_c[0], ls='dashed')
    # ax[1].axhline(2.5, c=tb_c[1], ls='dashed')
    # ax[1].axvline(3, c='k', ls='dashed', alpha=0.5)
    
    xticks = np.arange(0, 8+1, 1)
    ax[1].set(xlim=(0, 8), xticks=xticks, xlabel=r'$t\,[\,\textrm{Gyr}\,]$')
    ax[1].set(ylim=(0, 5), ylabel=r'$\textrm{galaxy size}\,[\,\textrm{kpc}\,]$')

    xminorticks = np.arange(0, 8+1, 0.25)
    ax[1].set_xticks(xminorticks, minor=True)
    
    yticks = np.arange(0, 5+1, 1)
    yminorticks = np.arange(0, 5+1, 0.25)
    ax[1].set_yticks(yticks)
    ax[1].set_yticks(yminorticks, minor=True)
    
    fig.tight_layout()
    fig.savefig('mass_size.pdf')
    
    return 0

if __name__ == '__main__':
    basepath = '/n/holylfs05/LABS/hernquist_lab/Users/abeane/GSEgas/'
    
    MWname = 'MW7iso'
    MWlvl = 'lvl4-Ngb64'
    
    GSEname = 'GSE4iso'
    GSElvl = 'lvl4'

    reset = True if len(sys.argv) == 1 else False
    
    if reset:
        i, subsMW = 0, {}
        while True:
            fname = basepath+'runs/'+MWname+'/'+MWlvl+'/output/fof_subhalo_tab_'+str(i).zfill(3)+'.hdf5'
            if os.path.isfile(fname):
                subsMW[i] = arepo.Subfind(fname)
            else:
                break
            i += 1
        
        i, subsGSE = 0, {}
        while True:
            fname = basepath+'runs/'+GSEname+'/'+GSElvl+'/output/fof_subhalo_tab_'+str(i).zfill(3)+'.hdf5'
            if os.path.isfile(fname):
                subsGSE[i] = arepo.Subfind(fname)
            else:
                break
            i += 1    
    else:
        subsMW = subsGSE = None
    
    if make_fig(subsMW, subsGSE, reset):
        print('failed to make fig')
    else:
        print('succeeded to make fig')
    