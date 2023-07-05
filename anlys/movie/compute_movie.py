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

def _runner(path, snap, COM, ptypes=[0, 2, 3, 4], 
            rng=[[-15, 15], [-15, 15]]):
    sn = arepo.Snapshot(path + '/output/', snap, 
                        parttype=ptypes, 
                        fields=['Coordinates', 'Velocities', 'Masses'],
                        combineFiles=True)
    
    # Get pos and vel of ptypes
    pos, vel = get_pos_vel(sn, [2, 3, 4])
    pos_gas, vel_gas = get_pos_vel(sn, [0])
    
    # if 'Nbody' in name:
    #     COM = np.array([0., 0., 0.])
    # #elif 'iso' in name:
    # #    COM = np.array([50., 50., 50.])
    # else:
    #     COM = np.array([sn.BoxSize])/2.

    # Compute projection
    Hxy_s, Hxz_s, Hxy_g, Hxz_g = compute_projections(sn, COM, rng=rng)
    # Hxy_s, Hxz_s, Hxy_g, Hxz_g = None, None, None, None

    # Grab time
    time = sn.Time.value

    # Package it all together
    output = (COM, Hxy_s, Hxz_s, Hxy_g, Hxz_g, time)
    
    return output

def run(path, name, fout, nsnap, nproc, rng, COM_key):

    if COM_key == 'BoxCenter':
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

    out = Parallel(n_jobs=nproc) (delayed(_runner)(path, i, COM_list[i], rng=rng) for i in tqdm(range(nsnap)))

    COM = np.array([out[i][0] for i in range(len(out))])
    Hxy_s = np.array([out[i][1] for i in range(len(out))])
    Hxz_s = np.array([out[i][2] for i in range(len(out))])
    Hxy_g = np.array([out[i][3] for i in range(len(out))])
    Hxz_g = np.array([out[i][4] for i in range(len(out))])
    time = np.array([out[i][5] for i in range(len(out))])

    # Output non-movie stuff
    output = {'COM': COM,
              'time': time}

    pickle.dump(output, open(fout+'.pickle', 'wb'))

    # Make movies
    make_movie(Hxy_s, time, Hxy_s.shape[1], 1E-3, 1E0, 'movies/'+fout+'_star_xy.mp4')
    make_movie(Hxz_s, time, Hxz_s.shape[1], 1E-3, 1E0, 'movies/'+fout+'_star_xz.mp4')
    make_movie(Hxy_g, time, Hxy_g.shape[1], 1E-4, 1E-1, 'movies/'+fout+'_gas_xy.mp4')
    make_movie(Hxz_g, time, Hxz_g.shape[1], 1E-4, 1E-1, 'movies/'+fout+'_gas_xz.mp4')

    # Output movie stuff
    # f = h5.File('proj-'+name+'.h5', 'w')
    # f.create_dataset("Hxy_s", data=Hxy_s)
    # f.create_dataset("Hxz_s", data=Hxz_s)
    # f.create_dataset("Hxy_g", data=Hxy_g)
    # f.create_dataset("Hxz_g", data=Hxz_g)
    # f.close()

if __name__ == '__main__':
    nproc = int(sys.argv[1])

    basepath = '../../runs/'
    

    Nbody = 'Nbody'
    fgMW05_fgGSE05 = 'fgGSE0.5_fgMW0.5'
    GSEiso = 'GSEiso_fg0.5_Z-1.2'
    GSEiso_Z0 = 'GSEiso_fg0.5_Z0'
    MWiso_Z0 = 'MWiso_fg0.5_Z0'
    MWiso_Z0_newS = 'MWiso_fg0.5_Z0-newSMUGGLE'
    MWiso_Z0_corona5 = 'MWiso_fg0.5_Z0_coronaMHG0.5'
    MW2iso_fg05 = 'MW2iso_fg0.5'
    MW3iso_fg05 = 'MW3iso_fg0.5'
    GSE2iso_fg07 = 'GSE2iso_fg0.7'

    rng0 = [[-80, 80], [-80, 80]]
    rng1 = [[-5, 5], [-5, 5]]
    rng2 = [[-8, 8], [-8, 8]]
    rng3 = [[-30, 30], [-30, 30]]
    rng4 = [[-15, 15], [-15, 15]]

    pair_list = [(Nbody, 'lvl4', rng0, 'Tot_COM'), # 0
                 (Nbody, 'lvl3', rng0, 'Tot_COM'), # 1
                 (fgMW05_fgGSE05, 'lvl4', rng0, 'Tot_COM'), # 2
                 (fgMW05_fgGSE05, 'lvl3', rng0, 'Tot_COM'), # 3
                 (fgMW05_fgGSE05, 'lvl4', rng3, 'Tot_COM'), # 4
                 (fgMW05_fgGSE05, 'lvl4', rng4, 'MW_COM'), # 5
                 (GSEiso, 'lvl4', rng1, 'BoxCenter'), # 6
                 (GSEiso_Z0, 'lvl4', rng1, 'BoxCenter'), # 7
                 (MWiso_Z0, 'lvl4', rng2, 'BoxCenter'), # 8
                 (MWiso_Z0_corona5, 'lvl4', rng1, 'BoxCenter'), # 9
                 (MW2iso_fg05, 'lvl4', rng2, 'BoxCenter'), # 10
                 (MW2iso_fg05, 'lvl3', rng2, 'BoxCenter'), # 11
                 (MW2iso_fg05, 'lvl4-nov', rng2, 'BoxCenter'), # 12
                 (MW2iso_fg05, 'lvl4-SFE1', rng2, 'BoxCenter'), # 13
                 (MW2iso_fg05, 'lvl4-novSFE1', rng2, 'BoxCenter'), # 14
                 (MW2iso_fg05, 'lvl4-hydrosoft', rng2, 'BoxCenter'), # 15
                 (MWiso_Z0_newS, 'lvl2-soft0.04', rng2, 'BoxCenter'), # 16
                 (MWiso_Z0_newS, 'lvl3', rng2, 'BoxCenter'), # 17
                 (MW3iso_fg05, 'lvl3', rng2, 'BoxCenter'), # 18
                 (MW3iso_fg05, 'lvl2', rng2, 'BoxCenter'), # 19
                 (MW3iso_fg05, 'lvl2-limiter', rng2, 'BoxCenter'), # 20
                 (MW3iso_fg05, 'lvl2-limiter2', rng2, 'BoxCenter'), # 21
                 (GSE2iso_fg07, 'lvl3', rng1, 'BoxCenter'), # 22
                 ]

    rng_list     = [                        p[2] for p in pair_list]
    rng_str_list = [str(rng).replace('[','').replace(']','').replace(', ','_') 
                    for rng in rng_list]

    name_list = [           p[0] + '-' + p[1] for p in pair_list]
    fout_list = [           p[0] + '-' + p[1] + '-rng_' + rng_str 
                            for p, rng_str in zip(pair_list, rng_str_list)]
    path_list = [basepath + p[0] + '/' + p[1] for p in pair_list]
    
    COM_key_list = [p[3] for p in pair_list]

    nsnap_list = [len(glob.glob(path+'/output/snapdir*/*.0.hdf5')) for path in path_list]

  
    i = int(sys.argv[2])
    path = path_list[i]
    name = name_list[i]
    fout = fout_list[i]
    nsnap = nsnap_list[i]
    rng = rng_list[i]
    COM_key = COM_key_list[i]

    out = run(path, name, fout, nsnap, nproc, rng, COM_key)
