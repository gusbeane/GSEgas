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

MW_Disk_Metal = -0.3
MW_CGM_Metal = -1.3
GSE_Disk_Metal = -1.2
GSE_CGM_Metal = -1.7


@njit
def _MC_inner_loop(ParentIDs, ParticleIDs):
    key = np.full(ParentIDs.shape, -1)
    
    ilist = np.arange(len(ParentIDs))
    jlist = np.arange(len(ParticleIDs))
    
    ilist = ilist[np.argsort(ParentIDs)]
    jlist = jlist[np.argsort(ParticleIDs)]
    
    i_ = j_ = 0
    
    while i_ < len(ilist):
        if ParentIDs[ilist[i_]] == ParticleIDs[jlist[j_]]:
            key[ilist[i_]] = jlist[j_]
            i_ += 1
        else:
            j_ += 1
    
    return key

def extract_MC_info(sn0_prop, sn):
    ParentIDs = sn.part5.ParentID
    TracerIDs = sn.part5.TracerID
    
    # Get standard quantities
    if sn.NumPart_Total[4] > 0:
        ParticleIDs = np.concatenate((sn.part0.ParticleIDs, sn.part4.ParticleIDs))
        Coordinates = np.concatenate((sn.part0.Coordinates, sn.part4.Coordinates))
        Velocities = np.concatenate((sn.part0.Velocities, sn.part4.Velocities))
        Masses = np.concatenate((sn.part0.Masses, sn.part4.Masses))
        PartType = np.concatenate((np.full(sn.NumPart_Total[0], 0), np.full(sn.NumPart_Total[4], 4)))
    else:
        ParticleIDs = sn.part0.ParticleIDs
        Coordinates = sn.part0.Coordinates
        Velocities = sn.part0.Velocities
        Masses = sn.part0.Masses
        PartType = np.full(sn.NumPart_Total[0], 0)
    
    # Get initial membership based on metallicity
    MW_Disk = 10.**(MW_Disk_Metal) * 0.0127
    MW_CGM = 10.**(MW_CGM_Metal) * 0.0127
    GSE_Disk = 10.**(GSE_Disk_Metal) * 0.0127
    GSE_CGM = 10.**(GSE_CGM_Metal) * 0.0127
    
    PartMembership0 = np.full(sn0_prop['NumPart_Total[0]'], -1)
    
    in_MW_Disk = np.abs(sn0_prop['part0.GFM_Metallicity'] - MW_Disk) < 1e-8
    in_MW_CGM = np.abs(sn0_prop['part0.GFM_Metallicity'] - MW_CGM) < 1e-8
    in_GSE_Disk = np.abs(sn0_prop['part0.GFM_Metallicity'] - GSE_Disk) < 1e-8
    in_GSE_CGM = np.abs(sn0_prop['part0.GFM_Metallicity'] - GSE_CGM) < 1e-8
    
    PartMembership0[in_MW_Disk] = 0
    PartMembership0[in_MW_CGM] = 1
    PartMembership0[in_GSE_Disk] = 2
    PartMembership0[in_GSE_CGM] = 3
    
    if np.any(PartMembership0 < 0):
        return -1
    
    key = _MC_inner_loop(ParentIDs, ParticleIDs)
    key0 = _MC_inner_loop(sn0_prop['part5.ParentID'], sn0_prop['part0.ParticleIDs'])
    if np.any(key < 0) or np.any(key0 < 0):
        return -1
    
    MCMembership0 = PartMembership0[key0]
    key0sn = _MC_inner_loop(TracerIDs, sn0_prop['part5.TracerID'])
    MCMembershipsn = MCMembership0[key0sn]
    
    MC_Prop = {}
    MC_Prop['Coordinates'] = Coordinates[key]
    MC_Prop['Velocities']  = Velocities[key]
    MC_Prop['Masses']      = Masses[key]
    MC_Prop['PartType']    = PartType[key]
    MC_Prop['Membership']  = MCMembershipsn
    
    return MC_Prop

@njit
def rodrigues_formula(k, v, theta):
    N = v.shape[0]
    v_rot = np.zeros(np.shape(v))
    
    ctheta = np.cos(theta)
    stheta = np.sin(theta)
    
    for i in range(N):
        v_rot[i] = v[i] * ctheta + np.cross(k, v[i]) * stheta + k * (np.dot(k, v[i])) * (1-ctheta)
    
    return v_rot

def add_rotate_pos(MC_Prop, MW_COM, MW_COMV, MW_AngMom):
    pos = MC_Prop['Coordinates'] - MW_COM
    vel = MC_Prop['Velocities'] - MW_COMV
    
    ang_mom = MW_AngMom

    angmom_dir = ang_mom/np.linalg.norm(ang_mom)
    theta = np.arccos(np.dot(angmom_dir, np.array([0, 0, 1])))
    k = np.cross(ang_mom, np.array([0, 0, 1.]))
    k /= np.linalg.norm(k)

    pos_rot = rodrigues_formula(k, pos.astype(np.float64), theta)
    vel_rot = rodrigues_formula(k, vel.astype(np.float64), theta)

    MC_Prop['RotatedCoordinates'] = pos_rot
    MC_Prop['RotatedVelocities']  = vel_rot
    
    return MC_Prop
    
def _runner(path, ic, name, COM_file, sn0_prop, snap):
    sn = arepo.Snapshot(path + '/output/', snap, 
                        combineFiles=True)
    
    MC_Prop = extract_MC_info(sn0_prop, sn)
    
    MW_COM    = COM_file['MW_COM'][snap]
    MW_COMV   = COM_file['MW_COMV'][snap]
    MW_AngMom = COM_file['MW_AngMom'][snap]
    
    MC_Prop = add_rotate_pos(MC_Prop, MW_COM, MW_COMV, MW_AngMom)
    
    MC_Prop['Time'] = sn.Time.value

    np.save(name+'/'+'MC_Prop_'+str(snap).zfill(3)+'.npy', MC_Prop)
    
    # Package it all together
    # output = (MC_Prop,)
    
    return None

def get_sn0_prop(path):
    sn0 = arepo.Snapshot(path + '/output/', 0, 
                        combineFiles=True)
    sn0_prop = {}
    
    sn0_prop['NumPart_Total[0]'] = sn0.NumPart_Total[0]
    sn0_prop['part0.GFM_Metallicity'] = sn0.part0.GFM_Metallicity
    sn0_prop['part5.ParentID'] = sn0.part5.ParentID
    sn0_prop['part0.ParticleIDs'] = sn0.part0.ParticleIDs
    sn0_prop['part5.TracerID'] = sn0.part5.TracerID
    
    return sn0_prop
    
def run(path, ic, name, nproc):
    
    if not os.path.exists(name):
        os.mkdir(name)
    
    basepath_COM = '../COM/'
    COM_fpath = basepath_COM + 'COM_' + name + '.npy'
    COM_file = np.load(COM_fpath, allow_pickle=True).item()
    nsnap = len(COM_file['MW_COM'])
    
    sn0_prop = get_sn0_prop(path)
    
    _ = Parallel(n_jobs=nproc) (delayed(_runner)(path, ic, name, COM_file, sn0_prop, i) for i in tqdm(range(nsnap)))

if __name__ == '__main__':
    nproc = int(sys.argv[1])

    basepath = '../../'

    MW3_GSE2_merge0 = 'MW3_MHG0.25_GSE2_MHG0.18'
    MW3_GSE2_merge1 = 'MW3_MHG0.25_GSE2_MHG0.18_Rcut30'
    MW3_GSE2_merge2 = 'MW3_MHG0.25_GSE2_MHG0.18_Rcut10'

    pair_list = [(MW3_GSE2_merge0, 'lvl4'), # 0
                 (MW3_GSE2_merge1, 'lvl4'), # 1
                 (MW3_GSE2_merge2, 'lvl4'), # 2
                 ]


    name_list = [           p[0] + '-' + p[1] for p in pair_list]
    path_list = [basepath + 'runs/' + p[0] + '/' + p[1] for p in pair_list]
    ic_list   = [basepath + 'ics/' + p[0] + '/' + p[1] for p in pair_list]
    
    i = int(sys.argv[2])
    path = path_list[i]
    name = name_list[i]
    ic = ic_list[i]

    out = run(path, ic, name, nproc)
