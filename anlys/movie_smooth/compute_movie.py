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

def run(path, name, fout, nres, nsnap, nproc, rng, COM_key):

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
        nsnap = len(COM_list)

    out = Parallel(n_jobs=nproc) (delayed(_runner)(path, i, COM_list[i], nres, rng=rng) for i in tqdm(range(nsnap)))

    COM = np.array([out[i][0] for i in range(len(out))])
    Hxy_s = np.array([out[i][1] for i in range(len(out))])
    Hxz_s = np.array([out[i][2] for i in range(len(out))])
    Hxy_g = np.array([out[i][3] for i in range(len(out))])
    Hxz_g = np.array([out[i][4] for i in range(len(out))])
    time = np.array([out[i][5] for i in range(len(out))])

    # Output non-movie stuff
    output = {'COM': COM,
              'time': time}

    # pickle.dump(output, open(fout+'.pickle', 'wb'))

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
    MW3iso_fg05 = 'MW3iso_fg0.5'
    GSE2iso_fg07 = 'GSE2iso_fg0.7'
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

    rng0 = [[-80, 80], [-80, 80]]
    rng1 = [[-5, 5], [-5, 5]]
    rng2 = [[-8, 8], [-8, 8]]
    rng3 = [[-30, 30], [-30, 30]]
    rng4 = [[-15, 15], [-15, 15]]

    pair_list = [(MW3iso_fg05, 'lvl3', rng2, 'BoxCenter'), # 0
                 (MW3iso_fg05, 'lvl2', rng2, 'BoxCenter'), # 1
                 (MW3iso_fg05, 'lvl2-limiter', rng2, 'BoxCenter'), # 2
                 (MW3iso_fg05, 'lvl2-limiter2', rng2, 'BoxCenter'), # 3
                 (GSE2iso_fg07, 'lvl3', rng1, 'BoxCenter'), # 4
                 (MW3iso_corona1, 'lvl4', rng2, 'BoxCenter'), # 5
                 (MW3iso_corona1, 'lvl3', rng2, 'BoxCenter'), # 6
                 (MW3iso_corona2, 'lvl4', rng2, 'BoxCenter'), # 7
                 (MW3iso_corona2, 'lvl3', rng2, 'BoxCenter'), # 8
                 (MW3iso_corona3, 'lvl4', rng0, 'BoxCenter'), # 9
                 (GSE2iso_corona1, 'lvl4', rng1, 'BoxCenter'), # 10
                 (GSE2iso_corona1, 'lvl3', rng1, 'BoxCenter'), # 11
                 (MW3_GSE2_merge0, 'lvl4', rng0, 'Tot_COM'), # 12
                 (MW3_GSE2_merge0, 'lvl4', rng4, 'GSE_COM'), # 13
                 (MW3_GSE2_merge1, 'lvl4', rng0, 'Tot_COM'), # 14
                 (MW3_GSE2_merge2, 'lvl4', rng0, 'Tot_COM'), # 15
                 (MW3_GSE2_merge2_pro, 'lvl4', rng0, 'Tot_COM'), # 16
                 (MW3_GSE2_merge0, 'lvl4', rng0, 'MW_COM'), # 17
                 (MW3_GSE2_merge1, 'lvl4', rng0, 'MW_COM'), # 18
                 (MW3_GSE2_merge2, 'lvl4', rng0, 'MW_COM'), # 19
                 (MW3_GSE2_merge3, 'lvl4', rng0, 'Tot_COM'), # 20
                 (MW3_GSE2_merge4, 'lvl4', rng0, 'Tot_COM'), # 21
                 (MW3iso_corona4, 'lvl4', rng0, 'BoxCenter'), # 22
                 ]

    rng_list     = [                        p[2] for p in pair_list]
    rng_str_list = [str(rng).replace('[','').replace(']','').replace(', ','_') 
                    for rng in rng_list]

    nres = int(sys.argv[3])
    
    name_list = [           p[0] + '-' + p[1] for p in pair_list]
    fout_list = [           p[0] + '-' + p[1] + '-rng_' + rng_str + '_' + p[3]
                            + '_nres' + str(nres)
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

    out = run(path, name, fout, nres, nsnap, nproc, rng, COM_key)
