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

metal_dict = {}
metal_dict[0] = {}
metal_dict[1] = {}

metal_dict[0]['MW_Disk'] = -0.3
metal_dict[0]['MW_CGM'] = -1.3
metal_dict[0]['GSE_Disk'] = -1.2
metal_dict[0]['GSE_CGM'] = -1.7

metal_dict[1]['MW_Disk'] = -0.3
metal_dict[1]['MW_CGM'] = -1.3
metal_dict[1]['GSE_Disk'] = -1.
metal_dict[1]['GSE_CGM'] = -1.2

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

def get_n_T(sn):
    UnitLength = sn.parameters.UnitLength_in_cm
    UnitMass = sn.parameters.UnitMass_in_g
    UnitVelocity = sn.parameters.UnitVelocity_in_cm_per_s

    UnitTime = UnitLength / UnitVelocity
    UnitEnergy = UnitMass * UnitVelocity**2

    HYDROGEN_MASSFRAC = 0.76
    GAMMA = 5./3.
    PROTONMASS = 1.67262178e-24
    BOLTZMANN = 1.38065e-16

    InternalEnergy = sn.part0.InternalEnergy.value
    ElectronAbundance = sn.part0.ElectronAbundance
    Density = sn.part0.Density.value
    
    mu = 4 * PROTONMASS / (1 + 3 * HYDROGEN_MASSFRAC + 4 * HYDROGEN_MASSFRAC * ElectronAbundance)
    T = (GAMMA - 1.) * (InternalEnergy / BOLTZMANN) * (UnitEnergy / UnitMass) * mu

    n = Density / mu
    n *= UnitMass/UnitLength**3
    
    return n, T
    
def extract_MC_info(sn0_prop, sn, metal):
    ParentIDs = sn.part5.ParentID
    TracerIDs = sn.part5.TracerID
    
    NumberDensity_Part0, Temperature_Part0 = get_n_T(sn)
    
    # Get standard quantities
    if sn.NumPart_Total[4] > 0:
        Temperature_Part4 = np.full(sn.NumPart_Total[4], -1)
        NumberDensity_Part4 = np.full(sn.NumPart_Total[4], -1)
        
        ParticleIDs = np.concatenate((sn.part0.ParticleIDs, sn.part4.ParticleIDs))
        Coordinates = np.concatenate((sn.part0.Coordinates, sn.part4.Coordinates))
        Velocities = np.concatenate((sn.part0.Velocities, sn.part4.Velocities))
        GFM_Metallicity = np.concatenate((sn.part0.GFM_Metallicity, sn.part4.GFM_Metallicity))
        NumberDensity = np.concatenate((NumberDensity_Part0, NumberDensity_Part4))
        Temperature = np.concatenate((Temperature_Part0, Temperature_Part4))
        PartType = np.concatenate((np.full(sn.NumPart_Total[0], 0), np.full(sn.NumPart_Total[4], 4)))
    else:
        ParticleIDs = sn.part0.ParticleIDs
        Coordinates = sn.part0.Coordinates
        Velocities = sn.part0.Velocities
        GFM_Metallicity = sn.part0.GFM_Metallicity
        NumberDensity = NumberDensity_Part0
        Temperature = Temperature_Part0
        PartType = np.full(sn.NumPart_Total[0], 0)
    
    # Get initial membership based on metallicity
    MW_Disk = 10.**(metal['MW_Disk']) * 0.0127
    MW_CGM = 10.**(metal['MW_CGM']) * 0.0127
    GSE_Disk = 10.**(metal['GSE_Disk']) * 0.0127
    GSE_CGM = 10.**(metal['GSE_CGM']) * 0.0127
    
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
    MC_Prop['GFM_Metallicity'] = GFM_Metallicity[key]
    MC_Prop['Temperature'] = Temperature[key]
    MC_Prop['NumberDensity'] = NumberDensity[key]
    MC_Prop['PartType']    = PartType[key]
    MC_Prop['Membership']  = MCMembershipsn.astype(np.int8)
    MC_Prop['TracerID'] = TracerIDs
    MC_Prop['ParentID'] = ParentIDs
    
    # now reorder based on TracerID
    key = np.argsort(TracerIDs)
    for field in ['Coordinates', 'Velocities', 'GFM_Metallicity', 'Temperature', 'NumberDensity',
                  'PartType', 'Membership', 'TracerID', 'ParentID']:
        MC_Prop[field] = MC_Prop[field][key]
    
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
    
def _runner(path, ic, name, COM_file, sn0_prop, snap, metal):
    sn = arepo.Snapshot(path + '/output/', snap, 
                        combineFiles=True)
    
    MC_Prop = extract_MC_info(sn0_prop, sn, metal)
    
    MW_COM    = COM_file['MW_COM'][snap]
    MW_COMV   = COM_file['MW_COMV'][snap]
    MW_AngMom = COM_file['MW_AngMom'][snap]
    
    MC_Prop = add_rotate_pos(MC_Prop, MW_COM, MW_COMV, MW_AngMom)
    
    t = h5.File(name+'/MC_Prop_'+str(snap).zfill(3)+'.h5', mode='w')
    
    t.create_dataset('PartType5/Coordinates', data=MC_Prop['Coordinates'])
    t.create_dataset('PartType5/Velocities', data=MC_Prop['Velocities'])
    t.create_dataset('PartType5/GFM_Metallicity', data=MC_Prop['GFM_Metallicity'])
    t.create_dataset('PartType5/NumberDensity', data=MC_Prop['NumberDensity'])
    t.create_dataset('PartType5/Temperature', data=MC_Prop['Temperature'])
    t.create_dataset('PartType5/PartType', data=MC_Prop['PartType'])
    t.create_dataset('PartType5/Membership', data=MC_Prop['Membership'])
    t.create_dataset('PartType5/TracerID', data=MC_Prop['TracerID'])
    t.create_dataset('PartType5/ParentID', data=MC_Prop['ParentID'])
    
    t.create_dataset('PartType5/RotatedCoordinates', data=MC_Prop['RotatedCoordinates'])
    t.create_dataset('PartType5/RotatedVelocities', data=MC_Prop['RotatedVelocities'])
    
    t.create_group('Header')
    t['Header'].attrs.create('Time', sn.Time.value)
    t['Header'].attrs.create('TracerMass', sn0_prop['TracerMass'])
    
    t.close()
    
    
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
    
    TotGasMass = np.sum(sn0.part0.mass.value)
    sn0_prop['TracerMass'] = TotGasMass / sn0.NumPart_Total[5]
    
    return sn0_prop
    
def run(path, ic, name, nproc, metal):
    
    if not os.path.exists(name):
        os.mkdir(name)
    
    basepath_COM = '../COM/'
    COM_fpath = basepath_COM + 'COM_' + name + '.npy'
    COM_file = np.load(COM_fpath, allow_pickle=True).item()
    nsnap = len(COM_file['MW_COM'])
    
    sn0_prop = get_sn0_prop(path)
    
    _ = Parallel(n_jobs=nproc) (delayed(_runner)(path, ic, name, COM_file, sn0_prop, i, metal) for i in tqdm(range(nsnap)))

if __name__ == '__main__':
    nproc = int(sys.argv[1])

    basepath = '../../'

    MW3iso_corona3 = 'MW3iso_fg0.7_MHG0.25_RC9'
    MW3iso_corona4 = 'MW3iso_fg0.7_MHG0.35_RC9'
    MW3_GSE2_merge0 = 'MW3_MHG0.25_GSE2_MHG0.18'
    MW3_GSE2_merge1 = 'MW3_MHG0.25_GSE2_MHG0.18_Rcut30'
    MW3_GSE2_merge2 = 'MW3_MHG0.25_GSE2_MHG0.18_Rcut10'
    MW3_GSE2_merge3 = 'MW3_MHG0.25_GSE2'
    MW3_GSE2_merge4 = 'MW3_MHG0.35_GSE2'
    MW3_GSE2_merge5 = 'MW3_MHG0.35_GSE2_Vvir110'
    MW3_GSE2_merge6 = 'MW3_MHG0.35_GSE2_e0.25'
    MW3_GSE5_merge0 = 'MW3_MHG0.25_GSE5'
    MW3_GSE3_merge0 = 'MW3_MHG0.35_GSE3'
    MW3_GSE2N_merge0 = 'MW3_MHG0.25_GSE2N'

    MW3_GSE6_merge0 = 'MW3_MHG0.25_GSE6'
    MW3_GSE6_merge1 = 'MW3_MHG0.25_GSE6_kick'

    MW4_GSE6_merge0 = 'MW4_MHG0.25_GSE6'
    MW4iso_corona0 = 'MW4iso_fg0.2_MHG0.25_RC9'

    MW4iso_corona4 = 'MW4iso_fg0.2_MHG0.15_RC9'
    MW4_GSE6_merge2 = 'MW4_MHG0.15_GSE6_kick'
    MW4_GSE6_merge3 = 'MW4_MHG0.15_GSE6_kick120'

    MW4_GSE2_MHG05 = 'MW4_MHG0.25_GSE2_MHG0.5'

    pair_list = [(MW3iso_corona3, 'lvl4', 0), # 0
                 (MW3iso_corona4, 'lvl4', 0), # 1
                 (MW3_GSE2_merge2, 'lvl4', 0), # 2
                 (MW3_GSE2_merge3, 'lvl4', 0), # 3
                 (MW3_GSE2_merge4, 'lvl4', 0), # 4
                 (MW3_GSE2_merge5, 'lvl4', 0), # 5
                 (MW3_GSE2_merge6, 'lvl4', 0), # 6
                 (MW3_GSE5_merge0, 'lvl4', 0), # 7
                 (MW3_GSE3_merge0, 'lvl4', 0), # 8
                 (MW3_GSE6_merge0, 'lvl4', 0), # 9
                 (MW3_GSE6_merge1, 'lvl4', 0), # 10
                 (MW3_GSE2N_merge0, 'lvl4', 0), # 11
                 (MW4_GSE6_merge0, 'lvl4', 0), # 12
                 (MW4iso_corona0, 'lvl4', 0), # 13
                 (MW4iso_corona4, 'lvl4', 0), # 14
                 (MW4_GSE6_merge2, 'lvl4', 0), # 15
                 (MW4_GSE6_merge3, 'lvl4', 0), # 16
                 (MW4_GSE2_MHG05, 'lvl4', 1), # 17
                 ]


    name_list  = [           p[0] + '-' + p[1] for p in pair_list]
    path_list  = [basepath + 'runs/' + p[0] + '/' + p[1] for p in pair_list]
    ic_list    = [basepath + 'ics/' + p[0] + '/' + p[1] for p in pair_list]
    metal_list = [p[2] for p in pair_list]
    
    i = int(sys.argv[2])
    path = path_list[i]
    name = name_list[i]
    ic = ic_list[i]
    metal = metal_dict[metal_list[i]]

    out = run(path, ic, name, nproc, metal)
