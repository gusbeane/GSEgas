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
        H_list = [Hxy_s, Hxz_s, Hxy_g, Hxz_g, Hxy_g, Hxz_g]
        vmin_list = [1E-3, 1E-3, 1E-4, 1E-4, 1E-5, 1E-5]
        vmax_list = [1E0, 1E0, 1E-1, 1E-1, 1E-2, 1E-2]
        fname_list = ['movies/'+fout+'_star_xy.mp4',
                      'movies/'+fout+'_star_xz.mp4',
                      'movies/'+fout+'_gas_xy.mp4',
                      'movies/'+fout+'_gas_xz.mp4',
                      'movies/'+fout+'_lowdens_gas_xy.mp4',
                      'movies/'+fout+'_lowdens_gas_xz.mp4',
                      ]

        _ = Parallel(n_jobs=4) (delayed(make_movie)(H, Time, H.shape[1], vmin, vmax, fname) 
                                for H, vmin, vmax, fname in zip(H_list, vmin_list, vmax_list, fname_list))
        
if __name__ == '__main__':
    basepath = '../../runs/'
    

    Nbody = 'Nbody'
    MW3iso_fg05 = 'MW3iso_fg0.5'
    GSE2iso_fg07 = 'GSE2iso_fg0.7'
    GSE3iso_fg07 = 'GSE3iso_fg0.7'
    GSE4iso_fg07 = 'GSE4iso_fg0.7'
    GSE5iso_fg07 = 'GSE5iso_fg0.7'
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
    MW3_GSE2_merge5 = 'MW3_MHG0.35_GSE2_Vvir110'
    MW3_GSE2_merge6 = 'MW3_MHG0.35_GSE2_e0.25'
    MW3_GSE5_merge0 = 'MW3_MHG0.25_GSE5'
    MW3_GSE3_merge0 = 'MW3_MHG0.35_GSE3'
    MW3_GSE6_merge0 = 'MW3_MHG0.25_GSE6'
    Seq1iso_fg07 = 'Seq1iso_fg0.7'
    Seq2iso_fg07 = 'Seq2iso_fg0.7'
    GSE2_Seq1 = 'GSE2_Seq1'
    GSE2_Seq1_R16 = 'GSE2_Seq1_R16'
    GSE2_Seq1_R24 = 'GSE2_Seq1_R24'
    GSE2_Seq1_R48 = 'GSE2_Seq1_R48'
    GSE2_Seq1_R48_eta05 = 'GSE2_Seq1_R48_eta0.5'
    GSE2_Seq1_R48_eta02 = 'GSE2_Seq1_R48_eta0.2'
    GSE6_Seq2_R48_eta02 = 'GSE6_Seq2_R48_eta0.2'
    GSE6_Seq2_R48_eta02_vz = 'GSE6_Seq2_R48_eta0.2_vz'
    GSE6_Seq2_R48_eta05 = 'GSE6_Seq2_R48_eta0.5'
    GSE6_Seq2_R48_eta05_vz = 'GSE6_Seq2_R48_eta0.5_vz'
    GSE6_Seq2_R96_eta05 = 'GSE6_Seq2_R96_eta0.5'
    GSE6_Seq2_R96_eta05_vz = 'GSE6_Seq2_R96_eta0.5_vz'
    GSE6_Seq2_R96_eta02 = 'GSE6_Seq2_R96_eta0.2'
    GSE6_Seq2_R96_eta02_vz = 'GSE6_Seq2_R96_eta0.2_vz'
    GSE6_Seq2_Rp30_Ri60_e06 = 'GSE6_Seq2_Rp30_Ri60_e0.6'
    GSE6_Seq2_Rp30_Ri60_e07 = 'GSE6_Seq2_Rp30_Ri60_e0.7'
    GSE6iso_fg07 = 'GSE6iso_fg0.7'
    MW3_GSE6_merge1 = 'MW3_MHG0.25_GSE6_kick'
    MW4_GSE6_merge0 = 'MW4_MHG0.25_GSE6'
    MW4iso_corona3 = 'MW4iso_fg0.2_MHG0.25_RC9'
    MW3_GSE2N_merge0 = 'MW3_MHG0.25_GSE2N'


    rng0 = [[-80, 80], [-80, 80]]
    rng1 = [[-5, 5], [-5, 5]]
    rng2 = [[-8, 8], [-8, 8]]
    rng3 = [[-30, 30], [-30, 30]]
    rng4 = [[-15, 15], [-15, 15]]
    rng5 = [[-40, 40], [-40, 40]]
    rng6 = [[-140, 140], [-140, 140]]

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
                 (MW3_GSE2_merge4, 'lvl4', rng0, 'Tot_COM'), # 23
                 (MW3_GSE2_merge5, 'lvl4', rng0, 'Tot_COM'), # 24
                 (MW3_GSE2_merge6, 'lvl4', rng0, 'Tot_COM'), # 25
                 (MW3_GSE3_merge0, 'lvl4', rng0, 'Tot_COM'), # 26
                 (MW3iso_corona4, 'lvl4', rng0, 'BoxCenter'), # 27
                 (GSE3iso_fg07, 'lvl4', rng1, 'BoxCenter'), # 28
                 (GSE2iso_fg07, 'lvl4', rng4, 'BoxCenter'), # 29
                 (GSE4iso_fg07, 'lvl4', rng4, 'BoxCenter'), # 30
                 (GSE5iso_fg07, 'lvl4', rng4, 'BoxCenter'), # 31
                 (GSE2iso_fg07, 'lvl4-fSN', rng4, 'BoxCenter'), # 32
                 (MW3_GSE5_merge0, 'lvl4', rng0, 'Tot_COM'), # 33
                 (Seq1iso_fg07, 'lvl4', rng1, 'BoxCenter'), # 34
                 (GSE2_Seq1, 'lvl4', rng4, 'BoxCenter'), # 35
                 (GSE2_Seq1_R16, 'lvl4', rng3, 'BoxCenter'), # 36
                 (GSE2_Seq1_R24, 'lvl4', rng5, 'BoxCenter'), # 37
                 (GSE2_Seq1_R48, 'lvl4', rng0, 'BoxCenter'), # 38
                 (GSE2_Seq1_R48_eta05, 'lvl4', rng0, 'BoxCenter'), # 39
                 (GSE6iso_fg07, 'lvl4', rng3, 'BoxCenter'), # 40
                 (GSE6iso_fg07, 'lvl4-adiabatic', rng3, 'BoxCenter'), # 41
                 (MW3_GSE6_merge0, 'lvl4', rng6, 'Tot_COM'), # 42
                 (GSE2_Seq1_R48_eta02, 'lvl4', rng0, 'BoxCenter'), # 43
                 (Seq2iso_fg07, 'lvl4', rng2, 'BoxCenter'), # 44
                 (GSE6_Seq2_R48_eta02, 'lvl4', rng0, 'BoxCenter'), # 45
                 (GSE6_Seq2_R48_eta02_vz, 'lvl4', rng0, 'BoxCenter'), # 46
                 (GSE6_Seq2_R48_eta05, 'lvl4', rng0, 'BoxCenter'), # 47
                 (GSE6_Seq2_R48_eta05_vz, 'lvl4', rng0, 'BoxCenter'), # 48
                 (GSE6_Seq2_R96_eta05, 'lvl4', rng0, 'BoxCenter'), # 49
                 (GSE6_Seq2_R96_eta05_vz, 'lvl4', rng0, 'BoxCenter'), # 50
                 (GSE6_Seq2_R96_eta02, 'lvl4', rng0, 'BoxCenter'), # 51
                 (GSE6_Seq2_R96_eta02_vz, 'lvl4', rng0, 'BoxCenter'), # 52
                 (GSE6_Seq2_Rp30_Ri60_e06, 'lvl4', rng0, 'BoxCenter'), # 53
                 (GSE6_Seq2_Rp30_Ri60_e07, 'lvl4', rng0, 'BoxCenter'), # 54
                 (MW3_GSE2_merge3, 'lvl3', rng0, 'Tot_COM'), # 55
                 (MW3_GSE6_merge1, 'lvl4', rng6, 'Tot_COM'), # 56
                 (MW4_GSE6_merge0, 'lvl4', rng0, 'Tot_COM'), # 57
                 (MW4iso_corona3, 'lvl4', rng0, 'Tot_COM'), # 58
                 (MW3_GSE2N_merge0, 'lvl4', rng6, 'Tot_COM'), # 59
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

    nsnap_list = [len(glob.glob(path+'/output/snapdir*/*.0.hdf5')) for path in path_list]

  
    path = path_list[i]
    name = name_list[i]
    fout = fout_list[i]
    nsnap = nsnap_list[i]
    rng = rng_list[i]
    COM_key = COM_key_list[i]
    
    out = run(snap, path, name, fout, nres, nsnap, rng, COM_key)
