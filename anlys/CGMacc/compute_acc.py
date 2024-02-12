import numpy as np
import arepo
import sys
from tqdm import tqdm
import glob
import os
import pickle
import h5py as h5
from numba import njit

from joblib import Parallel, delayed

basepath = '/n/holylfs05/LABS/hernquist_lab/Users/abeane/GSEgas/'

r0_list = np.arange(20, 200)
dr = 2.0

def accretion_rate(name, snap):
    fname = basepath + 'anlys/MC/'+name+'/MC_Prop_'+str(snap).zfill(3)+'.h5'
    f = h5.File(fname, mode='r')
    
    pos = f['PartType5/RotatedCoordinates'][:]
    vel = f['PartType5/RotatedVelocities'][:]
    memb = f['PartType5/Membership'][:]
    mass = f['Header'].attrs['TracerMass']
    
    r = np.linalg.norm(pos, axis=1)
    R = np.linalg.norm(pos[:,:2], axis=1)
    cphi = pos[:,0]/R
    sphi = pos[:,1]/R
    ctheta = pos[:,2]/r
    stheta = np.sqrt(1 - ctheta**2)
    
    vr = stheta * cphi * vel[:,0] + stheta * sphi * vel[:,1] + ctheta * vel[:,2]
    
    Mdotin = []
    Mdotout = []
    Menc = []
    
    for r0 in r0_list:
        key = np.logical_and(r > r0 - dr/2., r < r0 + dr/2.)
        key = np.logical_and(key, memb==1)
        keyin = np.logical_and(key, vr < 0)
        keyout = np.logical_and(key, vr > 0)
    
        Mdotin_ = mass * np.sum(vr[keyin]) / dr
        Mdotout_ = mass * np.sum(vr[keyout]) / dr
        
        enc = np.logical_and(r < r0, memb==1)
        Menc_ = mass * len(np.where(enc)[0])
        
        Mdotin.append(Mdotin_)
        Mdotout.append(Mdotout_)
        Menc.append(Menc_)
    
    Time = f['Header'].attrs['Time']
    
    return Time, Mdotin, Mdotout, Menc

def _runner(path, name, snap):
    Time, Mdotin, Mdotout, Menc = accretion_rate(name, snap)

    # Package it all together
    output = (Time, Mdotin, Mdotout, Menc)
    
    return output

def run(path, name, nsnap, nproc):

    fname = basepath + 'anlys/MC/'+name+'/MC_Prop_*.h5'
    nsnap = len(glob.glob(fname))
    
    out = Parallel(n_jobs=nproc) (delayed(_runner)(path, name, i) 
                                  for i in tqdm(range(nsnap)))

    Time        = np.array([out[i][0] for i in range(len(out))])
    Mdotin      = np.array([out[i][1] for i in range(len(out))])
    Mdotout     = np.array([out[i][2] for i in range(len(out))])
    Menc        = np.array([out[i][3] for i in range(len(out))])

    out = {'Time'    : Time,
           'Mdotin'  : Mdotin,
           'Mdotout' : Mdotout,
           'r0_list' : r0_list,
           'dr'      : dr,
           'Menc'    : Menc,
          }

    np.save('acc_'+name+'.npy', out)

if __name__ == '__main__':
    nproc = int(sys.argv[1])

    basepath = '../../'

    MW3iso_corona0   = 'MW3iso_fg0.7_MHG0.25_RC9'
    MW3_GSE6_merge0  = 'MW3_MHG0.25_GSE6'
    MW3_GSE2N_merge0 = 'MW3_MHG0.25_GSE2N'
    MW3_GSE6_merge1  = 'MW3_MHG0.25_GSE6_kick'

    MW4iso_corona4 = 'MW4iso_fg0.2_MHG0.15_RC9'
    MW4_GSE6_merge2 = 'MW4_MHG0.15_GSE6_kick'

    MW4iso_corona5 = 'MW4iso_fg0.2_MHG0.25_RC9'
    MW4_GSE2_corona = 'MW4_MHG0.25_GSE2_MHG0.5'

    pair_list = [(MW3iso_corona0, 'lvl4'), # 0
                 (MW3_GSE6_merge0, 'lvl4'), # 1
                 (MW3_GSE6_merge1, 'lvl4'), # 2
                 (MW3_GSE2N_merge0, 'lvl4'), # 3
                 (MW4iso_corona4, 'lvl4'), # 4
                 (MW4_GSE6_merge2, 'lvl4'), # 5
                 (MW4iso_corona5, 'lvl4'), # 6
                 (MW4_GSE2_corona, 'lvl4'), # 7
                 ]


    name_list = [           p[0] + '-' + p[1] for p in pair_list]
    path_list = [basepath + 'runs/' + p[0] + '/' + p[1] for p in pair_list]
    ic_list   = [basepath + 'ics/' + p[0] + '/' + p[1] for p in pair_list]
    
    nsnap_list = [len(glob.glob(path+'/output/snapdir*/*.0.hdf5')) for path in path_list]
  
    i = int(sys.argv[2])
    name = name_list[i]
    nsnap = nsnap_list[i]
    path = path_list[i]

    out = run(path, name, nsnap, nproc)
