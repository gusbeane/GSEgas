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

from joblib import Parallel, delayed

def get_pos_vel(sn, ptypes):
    pos = []
    vel = []
    for pt in ptypes:
        if sn.NumPart_Total[pt] > 0:
            part = getattr(sn, 'part'+str(pt))
            pos.append(part.pos.value)
            vel.append(part.vel.value)
    
    pos = np.concatenate(pos)
    vel = np.concatenate(vel)
    return pos, vel

def _runner(path, snap, COM, nres, ptypes=[0, 2, 3, 4], 
            rng=[[-15, 15], [-15, 15]]):
    sn = arepo.Snapshot(path + '/output/', snap, 
                        parttype=ptypes, 
                        fields=['Coordinates', 'Velocities', 'Masses', 'Density'],
                        combineFiles=True)
    
    # Compute projection
    Hxy_s, Hxz_s, Hxy_g, Hxz_g = compute_projections(sn, COM, nres, rng=rng)
    # Hxy_s, Hxz_s, Hxy_g, Hxz_g = None, None, None, None

    # Grab time
    time = sn.Time.value

    # Package it all together
    output = (COM, Hxy_s, Hxz_s, Hxy_g, Hxz_g, time)
    
    return output

def run(snap, path, name, fout, nres, nsnap, rng, COM_key):

    try:
        os.makedirs('frames/'+fout)
    except FileExistsError:
        pass
    
    if COM_key == 'BoxCenter':
        print('using BoxCenter')
        print('path to output: ', path+'/output/')
        sn = arepo.Snapshot(path + '/output/', 0, 
                        parttype=0, 
                        fields=['Coordinates'],
                        combineFiles=True)
        COM = np.array([sn.BoxSize, sn.BoxSize, sn.BoxSize])/2.
        COM_list = np.full((nsnap, 3), COM)    
    else:
        basepath_COM = '../COM/'
        COM_fpath = basepath_COM + 'COM_' + name + '.npy'
        COM_file = np.load(COM_fpath, allow_pickle=True).item()

        COM_list = COM_file[COM_key]
        nsnap = len(COM_list)
    
    # Make sure that the snap exists
    print('nsnap =', nsnap)
    if snap >= nsnap:
        return -1

    if snap >= 0:
    
        logging.debug('before _runner')
        out = _runner(path, snap, COM_list[snap], nres, rng=rng)

        COM = out[0]
        Hxy_s = out[1]
        Hxz_s = out[2]
        Hxy_g = out[3]
        Hxz_g = out[4]
        time = out[5]
        
        f = h5.File('frames/'+fout+'/frame'+str(snap).zfill(3)+'.h5', mode='w')
        
        f.create_dataset('COM', data=COM)
        f.create_dataset('Hxy_s', data=Hxy_s)
        f.create_dataset('Hxz_s', data=Hxz_s)
        f.create_dataset('Hxy_g', data=Hxy_g)
        f.create_dataset('Hxz_g', data=Hxz_g)
        f.create_dataset('Time', data=time)

    elif snap == -1:
        nsnap = len(glob.glob('frames/'+fout+'/frame*.h5'))
        Hxy_s = []
        Hxz_s = []
        Hxy_g = []
        Hxz_g = []
        Time = []
        
        for i in range(nsnap):
            if not os.path.isfile('frames/'+fout+'/frame'+str(i).zfill(3)+'.h5'):
                break

            f = h5.File('frames/'+fout+'/frame'+str(i).zfill(3)+'.h5', mode='r')
            
            Hxy_s.append(f['Hxy_s'][:])
            Hxz_s.append(f['Hxz_s'][:])
            Hxy_g.append(f['Hxy_g'][:])
            Hxz_g.append(f['Hxz_g'][:])
            
            Time.append(f['Time'][()])
        
            f.close()
        
        Hxy_s = np.array(Hxy_s)
        Hxz_s = np.array(Hxz_s)
        Hxy_g = np.array(Hxy_g)
        Hxz_g = np.array(Hxz_g)
        
        print(Hxy_s.shape)
        print(Time)
        
        # Make movies
        H_list = [Hxy_s, Hxz_s, Hxy_g, Hxz_g, Hxy_g, Hxz_g, Hxy_g, Hxz_g]
        vmin_list = [1E-3, 1E-3, 1E-4, 1E-4, 1E-5, 1E-5, 1E-6, 1E-6]
        vmax_list = [1E0, 1E0, 1E-1, 1E-1, 1E-2, 1E-2, 1E-3, 1E-3]
        fname_list = ['movies/'+fout+'_star_xy.mp4',
                      'movies/'+fout+'_star_xz.mp4',
                      'movies/'+fout+'_gas_xy.mp4',
                      'movies/'+fout+'_gas_xz.mp4',
                      'movies/'+fout+'_lowdens_gas_xy.mp4',
                      'movies/'+fout+'_lowdens_gas_xz.mp4',
                      'movies/'+fout+'_lowdens2_gas_xy.mp4',
                      'movies/'+fout+'_lowdens2_gas_xz.mp4',
                      ]

        _ = Parallel(n_jobs=4) (delayed(make_movie)(H, Time, H.shape[1], vmin, vmax, fname) 
                                for H, vmin, vmax, fname in zip(H_list, vmin_list, vmax_list, fname_list))
        
if __name__ == '__main__':
    basepath = '../../runs/'

    GSE2iso_corona2 = 'GSE2iso_fg0.7_MHG0.5_RC6.5'
    MW4_GSE2_MHG05 = 'MW4_MHG0.25_GSE2_MHG0.5'
    MW4_GSE2_MHG05_pro = 'MW4_MHG0.25_GSE2_MHG0.5_pro'
    MW4_GSE2N = 'MW4_MHG0.25_GSE2N'
    MW4iso_corona3 = 'MW4iso_fg0.2_MHG0.25_RC9'

    MW4iso_corona4_vphi01 = 'MW4iso_fg0.2_MHG0.25_RC9_vphi0.1'
    MW4_vphi01_GSE2_MHG05 = 'MW4_MHG0.25_vphi0.1_GSE2_MHG0.5'
    MW4_vphi01_GSE2N = 'MW4_MHG0.25_vphi0.1_GSE2N'
 
    MW4iso_corona4_vphi02 = 'MW4iso_fg0.2_MHG0.25_RC9_vphi0.2'
    MW4_vphi02_GSE2_MHG05 = 'MW4_MHG0.25_vphi0.2_GSE2_MHG0.5'
    MW4_vphi02_GSE2N = 'MW4_MHG0.25_vphi0.2_GSE2N'

    MW5iso = 'MW5iso'

    MW6iso = 'MW6iso'

    MW7iso = 'MW7iso'
    GSE4iso = 'GSE4iso'
    MW7_GSE4 = 'MW7_GSE4'

    rng0 = [[-80, 80], [-80, 80]]
    rng1 = [[-5, 5], [-5, 5]]
    rng2 = [[-8, 8], [-8, 8]]
    rng3 = [[-30, 30], [-30, 30]]
    rng4 = [[-15, 15], [-15, 15]]
    rng5 = [[-40, 40], [-40, 40]]
    rng6 = [[-140, 140], [-140, 140]]

    pair_list = [
                 (MW4iso_corona3, 'lvl4', rng0, 'Tot_COM'), # 0
                 (MW4iso_corona3, 'lvl4', rng5, 'Tot_COM'), # 1
                 (MW4iso_corona3, 'lvl4', rng6, 'Tot_COM'), # 2
                 (MW4_GSE2_MHG05, 'lvl4', rng0, 'Tot_COM'), # 3
                 (MW4_GSE2_MHG05, 'lvl4', rng6, 'Tot_COM'), # 4
                 
                 (MW4iso_corona3, 'lvl4-noB', rng6, 'Tot_COM'), # 5
                 (MW4_GSE2_MHG05, 'lvl4-noB', rng6, 'Tot_COM'), # 6
                 (MW4_GSE2N, 'lvl4-noB', rng6, 'Tot_COM'), # 7
                 (MW4_GSE2_MHG05_pro, 'lvl4-noB', rng6, 'Tot_COM'), # 8

                 (MW4iso_corona4_vphi01, 'lvl4-noB', rng6, 'Tot_COM'), # 9
                 (MW4_vphi01_GSE2_MHG05, 'lvl4-noB', rng6, 'Tot_COM'), # 10
                 (MW4_vphi01_GSE2N, 'lvl4-noB', rng6, 'Tot_COM'), # 11
                 
                 (MW4iso_corona4_vphi02, 'lvl4-noB', rng6, 'Tot_COM'), # 12
                 (MW4_vphi02_GSE2_MHG05, 'lvl4-noB', rng6, 'Tot_COM'), # 13
                 (MW4_vphi02_GSE2N, 'lvl4-noB', rng6, 'Tot_COM'), # 14
                 
                 (MW5iso, 'lvl5-beta05', rng0, 'Tot_COM'), # 15
                 (MW5iso, 'lvl5-beta06', rng0, 'Tot_COM'), # 16
                 (MW5iso, 'lvl5-beta067', rng0, 'Tot_COM'), # 17
                 (MW5iso, 'lvl5-beta07', rng0, 'Tot_COM'), # 18
                 (MW5iso, 'lvl5-beta08', rng0, 'Tot_COM'), # 19
                 
                 (MW6iso, 'lvl5-beta05-adi-dens', rng0, 'BoxCenter'), # 20
                 (MW6iso, 'lvl5-beta06-adi-dens', rng0, 'BoxCenter'), # 21
                 (MW6iso, 'lvl5-beta067-adi-dens', rng0, 'BoxCenter'), # 22
                 (MW6iso, 'lvl5-beta07-adi-dens', rng0, 'BoxCenter'), # 23
                 (MW6iso, 'lvl5-beta08-adi-dens', rng0, 'BoxCenter'), # 24
                 
                 (MW6iso, 'lvl5-beta05-dens', rng0, 'BoxCenter'), # 25
                 (MW6iso, 'lvl5-beta06-dens', rng0, 'BoxCenter'), # 26
                 (MW6iso, 'lvl5-beta067-dens', rng0, 'BoxCenter'), # 27
                 (MW6iso, 'lvl5-beta07-dens', rng0, 'BoxCenter'), # 28
                 (MW6iso, 'lvl5-beta08-dens', rng0, 'BoxCenter'), # 29
                 
                 (MW6iso, 'lvl5-beta08-fbar02-dens', rng0, 'BoxCenter'), # 30
                 (MW6iso, 'lvl5-beta08-fbar02-vphi03-dens', rng0, 'BoxCenter'), # 31
                 (MW6iso, 'lvl5-beta08-fbar02-vphi04-dens', rng0, 'BoxCenter'), # 32
                 
                 (MW6iso, 'lvl4-beta08-fbar02-vphi03-dens', rng0, 'BoxCenter'), # 33
                    
                 (MW6iso, 'lvl5-beta08-fbar008-dens-Z0', rng0, 'BoxCenter'), # 34
                 
                 (MW7iso, 'lvl5', rng0, 'BoxCenter'), # 35
                 (GSE4iso, 'lvl5', rng0, 'BoxCenter'), # 36
                 (MW7_GSE4, 'lvl5', rng6, 'BoxCenter'), # 37
                 ]

    i = int(sys.argv[1])
    nres = int(sys.argv[2])
    snap = int(sys.argv[3])
    
    rng_list     = [                        p[2] for p in pair_list]
    rng_str_list = [str(rng).replace('[','').replace(']','').replace(', ','_') 
                    for rng in rng_list]
    
    name_list = [           p[0] + '-' + p[1] for p in pair_list]
    fout_list = [           p[0] + '-' + p[1] + '-rng_' + rng_str + '_' + p[3]
                            + '_nres' + str(nres)
                            for p, rng_str in zip(pair_list, rng_str_list)]
    path_list = [basepath + p[0] + '/' + p[1] for p in pair_list]
    
    COM_key_list = [p[3] for p in pair_list]

    nsnap_list = []
    for path in path_list:
        nsnap1 = len(glob.glob(path+'/output/snapdir*/*.0.hdf5'))
        nsnap2 = len(glob.glob(path+'/output/snapshot_*.hdf5'))
        nsnap_list.append(max(nsnap1, nsnap2))
  
    path = path_list[i]
    name = name_list[i]
    fout = fout_list[i]
    nsnap = nsnap_list[i]
    rng = rng_list[i]
    COM_key = COM_key_list[i]
    
    out = run(snap, path, name, fout, nres, nsnap, rng, COM_key)
