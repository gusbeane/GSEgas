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

def get_COM(pos, vel, rmin=1, rmax=20, rfac=0.9):
    COM = np.mean(pos, axis=0)
    r = np.linalg.norm(pos-COM, axis=1)
    
    rcut = rmax
    while rcut > rmin:
        COM = np.mean(pos[r<rcut], axis=0)
        r = np.linalg.norm(pos-COM, axis=1)
        rcut *= rfac
    
    COMV = np.mean(vel[r<rcut], axis=0)
    
    return COM, COMV

def get_AngMom(pos, vel, mass, COM, COMV, rcut=8):
    r = np.linalg.norm(pos-COM)
    key = r < rcut
    
    ang = np.cross(pos-COM, vel-COMV)
    ang = np.sum(ang, axis=0)
    ang *= mass
    
    return ang

def _runner(path, ic, name, snap, ptypes=[2]):
    sn = arepo.Snapshot(path + '/output/', snap, 
                        parttype=ptypes, 
                        fields=['Coordinates', 'Velocities', 'Masses', 'ParticleIDs'],
                        combineFiles=True)
    
    if 'iso' in name:
        in_MW = in_GSE = None
        MW_pos = sn.part2.pos.value
        MW_vel = sn.part2.vel.value
        
        MW_COM, MW_COMV = get_COM(MW_pos, MW_vel)
        MW_AngMom = get_AngMom(MW_pos, MW_vel, sn.MassTable[2], MW_COM, MW_COMV)
        
        GSE_COM = GSE_COMV = GSE_AngMom = None
        
    else:    
        IDs = pickle.load(open(ic + '/IDs.p', 'rb'))

        # Get where each particle is
        in_MW = np.logical_and(sn.part2.id >= IDs['MW'][2][0], sn.part2.id <= IDs['MW'][2][1])
        in_GSE = np.logical_and(sn.part2.id >= IDs['GSE'][2][0], sn.part2.id <= IDs['GSE'][2][1])
    
        MW_pos  = sn.part2.pos.value[in_MW]
        MW_vel  = sn.part2.vel.value[in_MW]
        GSE_pos = sn.part2.pos.value[in_GSE]
        GSE_vel = sn.part2.vel.value[in_GSE]
    
        MW_COM, MW_COMV = get_COM(MW_pos, MW_vel)
        MW_AngMom = get_AngMom(MW_pos, MW_vel, sn.MassTable[2], MW_COM, MW_COMV)
    
        GSE_COM, GSE_COMV = get_COM(GSE_pos, GSE_vel)
        GSE_AngMom = get_AngMom(GSE_pos, GSE_vel, sn.MassTable[2], GSE_COM, GSE_COMV)
    
    Tot_COM  = np.mean(sn.part2.pos.value, axis=0)
    Tot_COMV = np.mean(sn.part2.vel.value, axis=0)
    
    Time = sn.Time.value

    # Package it all together
    output = (Tot_COM, Tot_COMV, 
              MW_COM, MW_COMV, MW_AngMom, 
              GSE_COM, GSE_COMV, GSE_AngMom, 
              Time)
    
    return output

def run(path, ic, name, nsnap, nproc):

    out = Parallel(n_jobs=nproc) (delayed(_runner)(path, ic, name, i) for i in tqdm(range(nsnap)))

    Tot_COM    = np.array([out[i][0] for i in range(len(out))])
    Tot_COMV   = np.array([out[i][1] for i in range(len(out))])
    MW_COM     = np.array([out[i][2] for i in range(len(out))])
    MW_COMV    = np.array([out[i][3] for i in range(len(out))])
    MW_AngMom  = np.array([out[i][4] for i in range(len(out))])
    GSE_COM    = np.array([out[i][5] for i in range(len(out))])
    GSE_COMV   = np.array([out[i][6] for i in range(len(out))])
    GSE_AngMom = np.array([out[i][7] for i in range(len(out))])
    Time       = np.array([out[i][8] for i in range(len(out))])

    out = {'Tot_COM'    : Tot_COM,
           'Tot_COMV'   : Tot_COMV,
           'MW_COM'     : MW_COM,
           'MW_COMV'    : MW_COMV,
           'MW_AngMom'  : MW_AngMom,
           'GSE_COM'    : GSE_COM,
           'GSE_COMV'   : GSE_COMV,
           'GSE_AngMom' : GSE_AngMom,
           'Time'       : Time}
    
    np.save('COM_'+name+'.npy', out)

if __name__ == '__main__':
    nproc = int(sys.argv[1])

    basepath = '../../'

    Nbody = 'Nbody'
    MW3iso_corona3 = 'MW3iso_fg0.7_MHG0.25_RC9'
    MW3_GSE2_merge2 = 'MW3_MHG0.25_GSE2_MHG0.18_Rcut10'
    MW3_GSE2_merge2_pro = 'MW3_MHG0.25_GSE2_MHG0.18_Rcut10_pro'
    MW3_GSE2_merge3 = 'MW3_MHG0.25_GSE2'
    MW3_GSE2_merge4 = 'MW3_MHG0.35_GSE2'

    pair_list = [(MW3iso_corona3, 'lvl4'), # 0
                 (MW3_GSE2_merge2, 'lvl4'), # 1
                 (MW3_GSE2_merge2_pro, 'lvl4'), # 2
                 (MW3_GSE2_merge3, 'lvl4'), # 3
                 (MW3_GSE2_merge4, 'lvl4'), # 4
                 ]


    name_list = [           p[0] + '-' + p[1] for p in pair_list]
    path_list = [basepath + 'runs/' + p[0] + '/' + p[1] for p in pair_list]
    ic_list   = [basepath + 'ics/' + p[0] + '/' + p[1] for p in pair_list]
    
    nsnap_list = [len(glob.glob(path+'/output/snapdir*/*.0.hdf5')) for path in path_list]
  
    i = int(sys.argv[2])
    path = path_list[i]
    name = name_list[i]
    nsnap = nsnap_list[i]
    ic = ic_list[i]

    out = run(path, ic, name, nsnap, nproc)
