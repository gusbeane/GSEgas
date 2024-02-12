import numpy as np
import arepo
import sys
from tqdm import tqdm
import glob
import os
import pickle
from projection import compute_projections
import h5py as h5
from numba import njit
from movie import make_movie
import logging
import matplotlib as mpl

from joblib import Parallel, delayed

def _runner(pathMC, snap, rbins = np.linspace(20, 180, 180),
                        vrbins = np.linspace(-200, 150, 180)):
    
    MC = h5.File(pathMC + '/MC_Prop_' + str(snap).zfill(3) + '.h5', mode='r')
    
    Hrvr, Hrvr_T = compute_projections(MC, rbins=rbins, vrbin=vrbins)
    
    time = MC['Header'].attrs['Time']
    
    MC.close()

    # Grab time
    
    output = (Hrvr, Hrvr_T, time)
    
    return output

def run(snap, nsnap, pathMC, fout, rbins = np.linspace(20, 180, 180),
                        vrbins = np.linspace(-200, 150, 180)):

    try:
        os.makedirs('frames/'+fout)
    except FileExistsError:
        pass
    
    if snap >= nsnap:
        return None

    if snap >= 0:
    
        out = _runner(pathMC, snap, rbins=rbins, vrbins=vrbins)

        Hrvr = out[0]
        Hrvr_T = out[1]
        time = out[2]
        
        f = h5.File('frames/'+fout+'/frame'+str(snap).zfill(3)+'.h5', mode='w')
        
        f.create_dataset('Hrvr', data=Hrvr)
        f.create_dataset('Hrvr_T', data=Hrvr_T)
        f.create_dataset('Time', data=time)
        
        f.close()

    elif snap == -1:
        nsnap = len(glob.glob('frames/'+fout+'/frame*.h5'))
        Hrvr = []
        Hrvr_T = []
        Time = []
        
        for i in range(nsnap):
            f = h5.File('frames/'+fout+'/frame'+str(i).zfill(3)+'.h5', mode='r')
            
            Hrvr.append(f['Hrvr'][:])
            Hrvr_T.append(f['Hrvr_T'][:])
            
            Time.append(f['Time'][()])
        
            f.close()
        
        Hrvr = np.array(Hrvr)
        Hrvr_T = np.array(Hrvr_T)
        
        # Make movies
        H_list = [Hrvr, Hrvr_T]
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
    basepathMC = '../../anlys/MC/'

    Nbody = 'Nbody'
    MW3iso_fg05 = 'MW3iso_fg0.5'
    GSE2iso_fg07 = 'GSE2iso_fg0.7'
    GSE3iso_fg07 = 'GSE3iso_fg0.7'
    MW3iso_corona1 = 'MW3iso_fg0.7_MHG0.1_RC30'
    MW3iso_corona2 = 'MW3iso_fg0.7_MHG0.15_RC9'
    MW3iso_corona3 = 'MW3iso_fg0.7_MHG0.25_RC9'
    MW3iso_corona4 = 'MW3iso_fg0.7_MHG0.35_RC9'
    GSE2iso_corona1 = 'GSE2iso_fg0.7_MHG0.18_RC6.5'
    MW3_GSE2_merge0 = 'MW3_MHG0.25_GSE2_MHG0.18'
    MW3_GSE2_merge1 = 'MW3_MHG0.25_GSE2_MHG0.18_Rcut30'
    MW3_GSE2_merge2 = 'MW3_MHG0.25_GSE2_MHG0.18_Rcut10'
    MW3_GSE2_merge2_pro = 'MW3_MHG0.25_GSE2_MHG0.18_Rcut10_pro'
    MW3_GSE2_merge3 = 'MW3_MHG0.25_GSE2'
    MW3_GSE2_merge4 = 'MW3_MHG0.35_GSE2'
    MW3_GSE6_merge0 = 'MW3_MHG0.25_GSE6'
    MW3_GSE6_merge1 = 'MW3_MHG0.25_GSE6_kick'

    rng0 = [[-80, 80], [-80, 80]]
    rng1 = [[-5, 5], [-5, 5]]
    rng2 = [[-8, 8], [-8, 8]]
    rng3 = [[-30, 30], [-30, 30]]
    rng4 = [[-15, 15], [-15, 15]]
    rng5 = [[-140, 140], [-140, 140]]

    pair_list = [(MW3iso_fg05, 'lvl3', rng2, 'BoxCenter'), # 0
                 (MW3iso_fg05, 'lvl2', rng2, 'BoxCenter'), # 1
                 (MW3iso_fg05, 'lvl2-limiter', rng2, 'BoxCenter'), # 2
                 (MW3iso_fg05, 'lvl2-limiter2', rng2, 'BoxCenter'), # 3
                 (GSE2iso_fg07, 'lvl4', rng1, 'BoxCenter'), # 4
                 (GSE2iso_fg07, 'lvl3', rng1, 'BoxCenter'), # 5
                 (MW3iso_corona1, 'lvl4', rng2, 'BoxCenter'), # 6
                 (MW3iso_corona1, 'lvl3', rng2, 'BoxCenter'), # 7
                 (MW3iso_corona2, 'lvl4', rng2, 'BoxCenter'), # 8
                 (MW3iso_corona2, 'lvl3', rng2, 'BoxCenter'), # 9
                 (MW3iso_corona3, 'lvl4', rng0, 'BoxCenter'), # 10
                 (MW3iso_corona4, 'lvl4', rng0, 'BoxCenter'), # 11
                 (GSE2iso_corona1, 'lvl4', rng1, 'BoxCenter'), # 12
                 (GSE2iso_corona1, 'lvl3', rng1, 'BoxCenter'), # 13
                 (MW3_GSE2_merge0, 'lvl4', rng0, 'Tot_COM'), # 14
                 (MW3_GSE2_merge0, 'lvl4', rng4, 'GSE_COM'), # 15
                 (MW3_GSE2_merge1, 'lvl4', rng0, 'Tot_COM'), # 16
                 (MW3_GSE2_merge2, 'lvl4', rng0, 'Tot_COM'), # 17
                 (MW3_GSE2_merge2_pro, 'lvl4', rng0, 'Tot_COM'), # 18
                 (MW3_GSE2_merge0, 'lvl4', rng0, 'MW_COM'), # 19
                 (MW3_GSE2_merge1, 'lvl4', rng0, 'MW_COM'), # 20
                 (MW3_GSE2_merge2, 'lvl4', rng0, 'MW_COM'), # 21
                 (MW3_GSE2_merge3, 'lvl4', rng0, 'Tot_COM'), # 22
                 (MW3_GSE2_merge3, 'lvl4', rng5, 'Tot_COM'), # 23
                 (MW3_GSE2_merge4, 'lvl4', rng0, 'Tot_COM'), # 24
                 (MW3iso_corona4, 'lvl4', rng0, 'BoxCenter'), # 25
                 (GSE3iso_fg07, 'lvl4', rng1, 'BoxCenter'), # 26
                 (MW3_GSE6_merge0, 'lvl4', rng5, 'Tot_COM'), # 27
                 (MW3_GSE6_merge1, 'lvl4', rng5, 'Tot_COM'), # 28
                 ]

    i = int(sys.argv[1])
    snap = int(sys.argv[2])
    
    rng_list     = [                        p[2] for p in pair_list]
    rng_str_list = [str(rng).replace('[','').replace(']','').replace(', ','_') 
                    for rng in rng_list]
    
    name_list = [           p[0] + '-' + p[1] for p in pair_list]
    fout_list = [           p[0] + '-' + p[1] + '-rng_' + rng_str + '_' + p[3]
                            for p, rng_str in zip(pair_list, rng_str_list)]
    path_list = [basepath + p[0] + '/' + p[1] for p in pair_list]
    pathMC_list = [basepathMC + p[0] + '-' + p[1] + '/' for p in pair_list]
    
    COM_key_list = [p[3] for p in pair_list]

    nsnap_list = [len(glob.glob(pathMC+'/MC_Prop_*.hdf5')) for pathMC in pathMC_list]

  
    path = path_list[i]
    pathMC = pathMC_list[i]
    name = name_list[i]
    fout = fout_list[i]
    nsnap = nsnap_list[i]
    rng = rng_list[i]
    COM_key = COM_key_list[i]
    
    out = run(snap, nsnap, pathMC, fout)
