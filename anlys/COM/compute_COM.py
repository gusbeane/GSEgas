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

def _runner(path, ic, name, snap, ptypes=[2]):
    sn = arepo.Snapshot(path + '/output/', snap, 
                        parttype=ptypes, 
                        fields=['Coordinates', 'Velocities', 'Masses', 'ParticleIDs'],
                        combineFiles=True)
    
    IDs = pickle.load(open(ic + '/IDs.p', 'rb'))

    # Get where each particle is
    in_MW = np.logical_and(sn.part2.id >= IDs['MW'][2][0], sn.part2.id <= IDs['MW'][2][1])
    in_GSE = np.logical_and(sn.part2.id >= IDs['GSE'][2][0], sn.part2.id <= IDs['GSE'][2][1])
    
    MW_pos  = sn.part2.pos.value[in_MW]
    MW_vel  = sn.part2.vel.value[in_MW]
    GSE_pos = sn.part2.pos.value[in_GSE]
    GSE_vel = sn.part2.vel.value[in_GSE]
    
    MW_COM, MW_COMV = get_COM(MW_pos, MW_vel)
    GSE_COM, GSE_COMV = get_COM(GSE_pos, GSE_vel)
    
    Tot_COM  = np.mean(sn.part2.pos.value, axis=0)
    Tot_COMV = np.mean(sn.part2.vel.value, axis=0)
    
    Time = sn.Time.value

    # Package it all together
    output = (Tot_COM, Tot_COMV, MW_COM, MW_COMV, GSE_COM, GSE_COMV, Time)
    
    return output

def run(path, ic, name, nsnap, nproc):

    out = Parallel(n_jobs=nproc) (delayed(_runner)(path, ic, name, i) for i in tqdm(range(nsnap)))

    Tot_COM  = np.array([out[i][0] for i in range(len(out))])
    Tot_COMV = np.array([out[i][1] for i in range(len(out))])
    MW_COM   = np.array([out[i][2] for i in range(len(out))])
    MW_COMV  = np.array([out[i][3] for i in range(len(out))])
    GSE_COM  = np.array([out[i][4] for i in range(len(out))])
    GSE_COMV = np.array([out[i][5] for i in range(len(out))])
    Time     = np.array([out[i][6] for i in range(len(out))])

    out = {'Tot_COM' : Tot_COM,
           'Tot_COMV': Tot_COMV,
           'MW_COM'  : MW_COM,
           'MW_COMV' : MW_COMV,
           'GSE_COM' : GSE_COM,
           'GSE_COMV': GSE_COMV,
           'Time'    : Time}
    
    np.save('COM_'+name+'.npy', out)

if __name__ == '__main__':
    nproc = int(sys.argv[1])

    basepath = '../../'

    Nbody = 'Nbody'
    fgMW05_fgGSE05 = 'fgGSE0.5_fgMW0.5'
    fgMW05_fgGSE05_COM = 'fgGSE0.5_fgMW0.5-COM'

    pair_list = [(Nbody, 'lvl4'), # 0
                 (Nbody, 'lvl3'), # 1
                 (fgMW05_fgGSE05, 'lvl4'), # 2
                 (fgMW05_fgGSE05, 'lvl3'), # 3
                 (fgMW05_fgGSE05_COM, 'lvl4'), # 4
                 (fgMW05_fgGSE05_COM, 'lvl3'), # 5
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
