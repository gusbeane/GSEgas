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
    MW3iso_corona4 = 'MW3iso_fg0.7_MHG0.35_RC9'
    MW3_GSE2_merge2 = 'MW3_MHG0.25_GSE2_MHG0.18_Rcut10'
    MW3_GSE2_merge2_pro = 'MW3_MHG0.25_GSE2_MHG0.18_Rcut10_pro'
    MW3_GSE2_merge3 = 'MW3_MHG0.25_GSE2'
    MW3_GSE2_merge4 = 'MW3_MHG0.35_GSE2'
    MW3_GSE2_merge5 = 'MW3_MHG0.35_GSE2_Vvir110'
    MW3_GSE2_merge6 = 'MW3_MHG0.35_GSE2_e0.25'
    MW3_GSE5_merge0 = 'MW3_MHG0.25_GSE5'
    MW3_GSE3_merge0 = 'MW3_MHG0.35_GSE3'
    MW3_GSE6_merge0 = 'MW3_MHG0.25_GSE6'
    MW3_GSE6_merge1 = 'MW3_MHG0.25_GSE6_kick'
    MW3_GSE2N_merge0 = 'MW3_MHG0.25_GSE2N'
    MW4_GSE6_merge0 = 'MW4_MHG0.25_GSE6'
    MW4iso_corona3 = 'MW4iso_fg0.2_MHG0.25_RC9'
    MW4_GSE6_merge2 = 'MW4_MHG0.15_GSE6_kick'
    MW4_GSE6_merge3 = 'MW4_MHG0.15_GSE6_kick120'

    MW4iso_corona4 = 'MW4iso_fg0.2_MHG0.15_RC9'
    MWN0_GSEN0_BH = 'MWN0_GSEN0-BH'
    MW4_GSE2_MHG05 = 'MW4_MHG0.25_GSE2_MHG0.5'

    MW4_GSE2N = 'MW4_MHG0.25_GSE2N'
    MW4_GSE2_MHG05_pro = 'MW4_MHG0.25_GSE2_MHG0.5_pro'

    MW4iso_corona4_vphi01 = 'MW4iso_fg0.2_MHG0.25_RC9_vphi0.1'
    MW4_vphi01_GSE2_MHG05 = 'MW4_MHG0.25_vphi0.1_GSE2_MHG0.5'
    MW4_vphi01_GSE2N = 'MW4_MHG0.25_vphi0.1_GSE2N'
    
    MW4iso_corona4_vphi02 = 'MW4iso_fg0.2_MHG0.25_RC9_vphi0.2'
    MW4_vphi02_GSE2_MHG05 = 'MW4_MHG0.25_vphi0.2_GSE2_MHG0.5'
    MW4_vphi02_GSE2N = 'MW4_MHG0.25_vphi0.2_GSE2N'

    MW5iso = 'MW5iso'
    
    pair_list = [
                 (MW4iso_corona3, 'lvl4'), # 0
                 (MW4_GSE2_MHG05, 'lvl4'), # 1
                 
                 (MW4iso_corona3, 'lvl4-noB'), # 2
                 (MW4_GSE2_MHG05, 'lvl4-noB'), # 3
                 (MW4_GSE2N, 'lvl4-noB'), # 4
                 (MW4_GSE2_MHG05_pro, 'lvl4-noB'), # 5
                 
                 (MW4iso_corona4_vphi01, 'lvl4-noB'), # 6
                 (MW4_vphi01_GSE2_MHG05, 'lvl4-noB'), # 7
                 (MW4_vphi01_GSE2N, 'lvl4-noB'), # 8
                 
                 (MW4iso_corona4_vphi02, 'lvl4-noB'), # 9
                 (MW4_vphi02_GSE2_MHG05, 'lvl4-noB'), # 10
                 (MW4_vphi02_GSE2N, 'lvl4-noB'), # 11
                 
                 (MW5iso, 'lvl5-beta05'), # 12
                 (MW5iso, 'lvl5-beta06'), # 13
                 (MW5iso, 'lvl5-beta067'), # 14
                 (MW5iso, 'lvl5-beta07'), # 15
                 (MW5iso, 'lvl5-beta08'), # 16
                ]


    name_list = [           p[0] + '-' + p[1] for p in pair_list]
    path_list = [basepath + 'runs/' + p[0] + '/' + p[1] for p in pair_list]
    ic_list   = [basepath + 'ics/' + p[0] + '/' + p[1][:4] for p in pair_list]
    # p[1][:4] in ics is to remove anything like -GFM from the lvl, since the ics are just for
    # lvl4 and lvl4-GFM is a runtime modification (ie.. lvl4-GFM uses the same ics as lvl4)

    nsnap_list = []
    for path in path_list:
        nsnap0 = len(glob.glob(path+'/output/snap*.hdf5'))
        nsnap1 = len(glob.glob(path+'/output/snapdir*/*.0.hdf5'))
        nsnap_list.append(max(nsnap0, nsnap1))

    i = int(sys.argv[2])
    path = path_list[i]
    name = name_list[i]
    nsnap = nsnap_list[i]
    ic = ic_list[i]

    out = run(path, ic, name, nsnap, nproc)
