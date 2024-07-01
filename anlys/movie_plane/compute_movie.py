import numpy as np
import arepo
import sys
from tqdm import tqdm
import glob
import os
import pickle
import h5py as h5
from numba import njit
from movie import make_movie
import logging
import matplotlib as mpl
from scipy.stats import binned_statistic_2d

from joblib import Parallel, delayed

def get_logFeH_logMgFe(sn, ptype=4):
    GFM_SOLAR_ABUNDANCE_HYDROGEN = 0.7388
    GFM_SOLAR_ABUNDANCE_MAGNESIUM = 0.0007
    GFM_SOLAR_ABUNDANCE_IRON   =   0.0013
    
    part = getattr(sn, 'part'+str(ptype))
    
    FeH = part.GFM_Metals[:,8] / GFM_SOLAR_ABUNDANCE_IRON
    logFeH = np.log10(FeH)

    MgH = part.GFM_Metals[:,6] / GFM_SOLAR_ABUNDANCE_MAGNESIUM
    MgFe = MgH/FeH
    logMgH = np.log10(MgH)
    logMgFe = np.log10(MgFe)
    
    return logFeH, logMgH, logMgFe

def _runner(path, snap, rcut):
    
    logFeH_bins  = np.linspace(-0.9, 0.6, 50)
    logMgFe_bins = np.linspace(-0.05, 0.45, 50) + 0.2
    
    try:
        sn = arepo.Snapshot(path+'/output', snap, parttype=[0])
        sub = arepo.Subfind(path+'/output', snap)
    except:
        return None
    
    logFeH, logMgH, logMgFe = get_logFeH_logMgFe(sn, ptype=0)
    r = np.linalg.norm(sn.part0.pos.value - sub.SubhaloPos[0], axis=1)
    
    mass = sn.part0.mass.value
    sfr = sn.part0.sfr.value
    
    key = r < rcut
    massbin, x_edge, y_edge, _ = binned_statistic_2d(logFeH[key], logMgFe[key], mass[key], 
                                                     statistic='sum', bins=[logFeH_bins, logMgFe_bins])
    
    sfrbin, x_edge, y_edge, _ = binned_statistic_2d(logFeH[key], logMgFe[key], sfr[key], 
                                                     statistic='sum', bins=[logFeH_bins, logMgFe_bins])
    
    output = (massbin, sfrbin, sn.Time.value)
    
    return output

def run(snap, path, fout, rcut=10):
    try:
        os.makedirs('frames/'+fout)
    except FileExistsError:
        pass
    
    if snap >= 0:
    
        out = _runner(path, snap, rcut)
        if out is None:
            return None

        massmap = out[0]
        sfrmap = out[1]
        Time = out[2]
        
        f = h5.File('frames/'+fout+'/frame'+str(snap).zfill(3)+'.h5', mode='w')
        
        f.create_dataset('Mmap', data=massmap)
        f.create_dataset('SFRmap', data=sfrmap)
        f.create_dataset('Time', data=Time)
        
        f.close()

    elif snap == -1:
        nsnap = len(glob.glob('frames/'+fout+'/frame*.h5'))
        Massmap = []
        SFRmap = []
        Time = []
        
        for i in range(nsnap):
            f = h5.File('frames/'+fout+'/frame'+str(i).zfill(3)+'.h5', mode='r')
            
            Massmap.append(f['Mmap'][:])
            SFRmap.append(f['SFRmap'][:])
            Time.append(f['Time'][()])
            
            f.close()
        
        Massmap = np.array(Massmap)
        SFRmap = np.array(SFRmap)
        Time = np.array(Time)
        
        # Make movies
        H_list = [Massmap, SFRmap]
        norm_list = [mpl.colors.LogNorm(vmin=8.7E-7, vmax=8.7E-4),
                     mpl.colors.Normalize(vmin=4, vmax=5.5),
                     ]
        
        fname_list = ['movies/'+fout+'_rho_rvr.mp4',
                      'movies/'+fout+'_T_rvr.mp4',
                       ]
        
        cmap_list = ['viridis', 'bwr']
        extent = [rbins[0], rbins[-1], vrbins[0], vrbins[-1]]

        _ = Parallel(n_jobs=2) (delayed(make_movie)(H, Time, H.shape[1], H.shape[2], norm, 
                                                    fname, cmap, extent) 
                                for H, norm, fname, cmap in 
                                zip(H_list, norm_list, fname_list, cmap_list))
        
if __name__ == '__main__':
    basepath = '../../runs/'

    MW7_GSE4 = 'MW7_GSE4'
    
    pair_list = [(MW7_GSE4, 'lvl5-denscut-Ngb64-steep1'), # 0
                 ]

    i = int(sys.argv[1])
    snap = int(sys.argv[2])
    
    name_list = [           p[0] + '-' + p[1] for p in pair_list]
    fout_list = [           p[0] + '-' + p[1] for p in pair_list]
    path_list = [basepath + p[0] + '/' + p[1] for p in pair_list]

    path = path_list[i]
    name = name_list[i]
    fout = fout_list[i]
    
    out = run(snap, path, fout)
