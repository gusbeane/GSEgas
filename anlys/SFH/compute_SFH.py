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

def _runner(path, COM, name, snap, ptypes=[0], rcut=15):
    sn = arepo.Snapshot(path + '/output/', snap, 
                        parttype=ptypes, 
                        fields=['Coordinates', 'Velocities', 'Masses', 'ParticleIDs', 'StarFormationRate'],
                        combineFiles=True)
    
    IDs = pickle.load(open(ic + '/IDs.p', 'rb'))
    
    MW_COM = COM['MW_COM'][snap]
    GSE_COM = COM['GSE_COM'][snap]

    in_MW = sn.part0.pos.value - MW_COM
    in_MW = np.linalg.norm(in_MW, axis=1)
    in_MW = in_MW < rcut
    
    in_GSE = sn.part0.pos.value - GSE_COM
    in_GSE = np.linalg.norm(in_GSE, axis=1)
    in_GSE = np.logical_and(in_GSE < rcut, np.logical_not(in_MW))
    
    SFR_MW = np.sum(sn.part0.sfr.value[in_MW])
    SFR_GSE = np.sum(sn.part0.sfr.value[in_GSE])
    
    Time = sn.Time.value
    
    output = (SFR_MW, SFR_GSE, Time)
    
    return output

def run(path, ic, name, nsnap, nproc):

    out = Parallel(n_jobs=nproc) (delayed(_runner)(path, ic, name, i) for i in tqdm(range(nsnap)))

    SFR_MW  = np.array([out[i][0] for i in range(len(out))])
    SFR_GSE = np.array([out[i][1] for i in range(len(out))])
    Time    = np.array([out[i][2] for i in range(len(out))])

    out = {'SFR_MW' : SFR_MW,
           'SFR_GSE': SFR_GSE,
           'Time'   : Time}
    
    np.save('SFH_'+name+'.npy', out)

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
    # nsnap = nsnap_list[i]
    ic = ic_list[i]
    
    COM = np.load(basepath+'/anlys/COM/COM_'+name+'.npy', allow_pickle=True).item()

    # pull nsnap from COM in case new snaps have been written since COM was computed
    nsnap = COM['Tot_COM'].shape[0]

    out = run(path, COM, name, nsnap, nproc)
