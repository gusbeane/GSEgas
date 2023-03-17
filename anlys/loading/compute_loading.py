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

def mass_loading(sn, z0):
    center = np.array([sn.BoxSize, sn.BoxSize, sn.BoxSize])/2.
    pos = sn.part0.pos.value - center
    vel = sn.part0.vel.value
    mass = sn.part0.mass.value

    dr = 0.1 * z0

    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]

    r = np.linalg.norm(pos, axis=1)
    phi = np.arctan2(y, x)
    theta = np.arccos(z/r)

    rhat = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
    vel_r = rhat[0]*vel[:,0] + rhat[1]*vel[:,1] + rhat[2]*vel[:,2]

    key = np.logical_and(np.abs(z) > z0-dr/2, np.abs(z) < z0+dr/2)
    key_out = np.logical_and(key, vel_r > 0)
    key_in = np.logical_and(key, vel_r < 0)

    Mdot_out = np.sum(mass[key_out]*vel_r[key_out])/dr
    Mdot_in = np.sum(mass[key_in]*vel_r[key_in])/dr
    
    return Mdot_out, Mdot_in

def _runner(path, name, snap, ptypes=[0]):
    sn = arepo.Snapshot(path + '/output/', snap, 
                        parttype=ptypes, 
                        fields=['Coordinates', 'Velocities', 'Masses', 'StarFormationRate'],
                        combineFiles=True)
    
    Mdot_out_1,  Mdot_in_1  = mass_loading(sn, 1)
    Mdot_out_3,  Mdot_in_3  = mass_loading(sn, 3)
    Mdot_out_5,  Mdot_in_5  = mass_loading(sn, 5)
    Mdot_out_10, Mdot_in_10 = mass_loading(sn, 10)
    
    SFR = np.sum(sn.part0.sfr)
    SFR_CodeUnits = SFR * (9.7779E8/1E10)
    Time = sn.Time.value

    # Package it all together
    output = (Mdot_out_1,  Mdot_in_1,
              Mdot_out_3,  Mdot_in_3,
              Mdot_out_5,  Mdot_in_5,
              Mdot_out_10, Mdot_in_10,
              SFR, SFR_CodeUnits, Time)
    
    return output

def run(path, name, nsnap, nproc):

    out = Parallel(n_jobs=nproc) (delayed(_runner)(path, name, i) for i in tqdm(range(nsnap)))

    Mdot_out_1    = np.array([out[i][0] for i in range(len(out))])
    Mdot_in_1     = np.array([out[i][1] for i in range(len(out))])
    Mdot_out_3    = np.array([out[i][2] for i in range(len(out))])
    Mdot_in_3     = np.array([out[i][3] for i in range(len(out))])
    Mdot_out_5    = np.array([out[i][4] for i in range(len(out))])
    Mdot_in_5     = np.array([out[i][5] for i in range(len(out))])
    Mdot_out_10   = np.array([out[i][6] for i in range(len(out))])
    Mdot_in_10    = np.array([out[i][7] for i in range(len(out))])
    SFR           = np.array([out[i][8] for i in range(len(out))])
    SFR_CodeUnits = np.array([out[i][9] for i in range(len(out))])
    Time          = np.array([out[i][10] for i in range(len(out))])

    out = {'Mdot_out_1'    : Mdot_out_1,
           'Mdot_in_1'     : Mdot_in_1,
           'Mdot_out_3'    : Mdot_out_3,
           'Mdot_in_3'     : Mdot_in_3,
           'Mdot_out_5'    : Mdot_out_5,
           'Mdot_in_5'     : Mdot_in_5,
           'Mdot_out_10'   : Mdot_out_10,
           'Mdot_in_10'    : Mdot_in_10,
           'SFR'           : SFR,
           'SFR_CodeUnits' : SFR_CodeUnits,
           'Time'          : Time}
    
    np.save('loading_'+name+'.npy', out)

if __name__ == '__main__':
    nproc = int(sys.argv[1])

    basepath = '../../'

    MWfg05Z0 = 'MWiso_fg0.5_Z0'

    pair_list = [(MWfg05Z0, 'lvl4'), # 0
                 ]


    name_list = [           p[0] + '-' + p[1] for p in pair_list]
    path_list = [basepath + 'runs/' + p[0] + '/' + p[1] for p in pair_list]
    # ic_list   = [basepath + 'ics/' + p[0] + '/' + p[1] for p in pair_list]
    
    nsnap_list = [len(glob.glob(path+'/output/snapdir*/*.0.hdf5')) for path in path_list]
  
    i = int(sys.argv[2])
    path = path_list[i]
    name = name_list[i]
    nsnap = nsnap_list[i]

    out = run(path, name, nsnap, nproc)
