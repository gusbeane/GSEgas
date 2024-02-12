import numpy as np
import arepo
from numba import njit
import h5py as h5

from joblib import Parallel, delayed
import illustris_python as il

def get_time(time, redshift=False, 
             Omega0=0.3089, 
             OmegaLambda=0.6911,
             HubbleParam=0.6774):
    HUBBLE = 3.2407789e-18
    SEC_PER_MEGAYEAR = 3.15576e13
    
    if redshift:
        a = 1./(1.+time)
    else:
        a = time
    
    fac = 2. / (3. * np.sqrt(OmegaLambda))
    ans = fac * np.arcsinh(np.sqrt(a**3 * OmegaLambda/Omega0))

    ans /= HUBBLE * HubbleParam
    ans /= SEC_PER_MEGAYEAR * 1000
    
    return ans

@njit
def rodrigues_formula(k, v, theta):
    N = v.shape[0]
    v_rot = np.zeros(np.shape(v))
    
    ctheta = np.cos(theta)
    stheta = np.sin(theta)
    
    for i in range(N):
        v_rot[i] = v[i] * ctheta + np.cross(k, v[i]) * stheta + k * (np.dot(k, v[i])) * (1-ctheta)
    
    return v_rot

def _load_single_numpart(outputDir, snapnum, i):
    snapnum_str = str(snapnum).zfill(3)
    fname = outputDir + 'snapdir_' + snapnum_str + '/snap_' + snapnum_str + '.'+str(i)+'.hdf5'
    
    with h5.File(fname, mode='r') as f:
        NumPart_ThisFile = f['Header'].attrs['NumPart_ThisFile']
    
    return NumPart_ThisFile

def _load_single_snap(outputDir, snapnum, i, numToLoad, types=[0, 4, 5], fields=None):
    out = {}
    if np.all(numToLoad == 0):
        return out
    
    snapnum_str = str(snapnum).zfill(3)
    fname = outputDir + 'snapdir_' + snapnum_str + '/snap_' + snapnum_str + '.'+str(i)+'.hdf5'
    
    with h5.File(fname, mode='r') as f:
        for pt in types:
            ptype = 'PartType' + str(pt)
            if numToLoad[pt] > 0:
                out[ptype] = {}
                for fld in f[ptype].keys():
                    
                    if fields is not None:
                        if fld not in fields:
                            continue
                
                    out[ptype][fld] = f[ptype][fld][:numToLoad[pt]]
    return out
    
def load_zoom_group(outputDir, snapnum, nproc=16, types=[0, 4, 5], fields=None):
    NTypes = 6
    snapnum_str = str(snapnum).zfill(3)
        
    fof0_fname = outputDir + 'groups_' + snapnum_str + '/fof_subhalo_tab_' + snapnum_str + '.0.hdf5'
    with h5.File(fof0_fname, mode='r') as fof0:
        GroupLenType = fof0['Group/GroupLenType'][0]
    
    snap0_fname = outputDir + 'snapdir_' + snapnum_str + '/snap_' + snapnum_str + '.0.hdf5'
    with h5.File(snap0_fname, mode='r') as snap0:
        NumFiles = snap0['Header'].attrs['NumFilesPerSnapshot']
    
    NumToLoad = np.zeros((NumFiles, NTypes), dtype=int)
    NumPartLeft = np.copy(GroupLenType)
    
    NumPart_ThisFile_all = Parallel(n_jobs=nproc) (delayed(_load_single_numpart)(outputDir, snapnum, i) 
                                  for i in range(NumFiles))
    
    for i in range(NumFiles):
        NumPart_ThisFile = NumPart_ThisFile_all[i]
        for pt in range(NTypes):
            NumToLoad[i][pt] = min(NumPartLeft[pt], NumPart_ThisFile[pt])
            NumPartLeft[pt] -= NumToLoad[i][pt]
                
        if np.all(NumPartLeft == 0):
            break
        
    if np.any(NumPartLeft < 0):
        print('error')
        return -1

    # get all snapnums that have particles we want to load
    snapnum_toload = np.where(np.any(NumToLoad, axis=1))[0]
    
    # now we load
    out_all = Parallel(n_jobs=nproc) (delayed(_load_single_snap)(outputDir, snapnum, i, NumToLoad[i], types=types, fields=fields) 
                                  for i in snapnum_toload)
    
    out = {}
    for o in out_all:
        for ptype in o.keys():
            if ptype not in out.keys():
                out[ptype] = {}
            
            for fld in o[ptype].keys():
                if fld not in out[ptype].keys():
                    out[ptype][fld] = []
                
                out[ptype][fld].append(o[ptype][fld])
    
    for pt in range(NTypes):
        ptype = 'PartType' + str(pt)
        if ptype not in out.keys():
            continue

        for fld in out[ptype].keys():
            out[ptype][fld] = np.concatenate(out[ptype][fld])
            assert GroupLenType[pt] == len(out[ptype][fld])
    
    return out

class EmptyClass:
    pass

def load_galaxy(snapnum, subID_99, rhalf_fac=2, phys=True, load_DM=False,
                TNGbase='/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/L35n2160TNG/output/'):
    # create output object
    galaxy = EmptyClass()
    
    # subID_99 is the subhalo ID at snapshot 99
    header = arepo.Snapshot(TNGbase, snapnum, onlyHeader=True)
    h = header.HubbleParam
    a = header.Time
    galaxy.header = header
    
    treeMPB = il.sublink.loadTree(TNGbase, 99, subID_99, onlyMPB=True)
    galaxy.treeMPB = treeMPB
    
    # retrieve subID and subhalo at given snapnum
    k_snap = int(np.where(treeMPB['SnapNum'] == snapnum)[0])
    subID = treeMPB['SubfindID'][k_snap]
    subhalo = il.groupcat.loadSingle(TNGbase, snapnum, subhaloID=subID)
    haloID = subhalo['SubhaloGrNr']
    halo = il.groupcat.loadSingle(TNGbase, snapnum, haloID=subhalo['SubhaloGrNr'])
    galaxy.subID = subID
    galaxy.subhalo = subhalo
    galaxy.haloID = haloID
    galaxy.halo = halo
    
    # load snapshot
    sn = EmptyClass()
    sn.part0 = il.snapshot.loadHalo(TNGbase, snapnum, haloID, 0)
    if load_DM:
        sn.part1 = il.snapshot.loadHalo(TNGbase, snapnum, haloID, 1)
    sn.part4 = il.snapshot.loadHalo(TNGbase, snapnum, haloID, 4)
    sn.part5 = il.snapshot.loadHalo(TNGbase, snapnum, haloID, 5)
    
    # get COM, COMV, and ang mom of stars
    COM = subhalo['SubhaloPos']
    pos = sn.part4['Coordinates'] - COM
    r = np.linalg.norm(pos, axis=1)
    
    # sort out stars within rhalf_fac*rhalf of origin
    rhalf = subhalo['SubhaloHalfmassRadType'][4]
    in_rhalf = r < rhalf_fac * rhalf
    is_star = sn.part4['GFM_StellarFormationTime'] > 0
    is_star_in_rhalf = np.logical_and(is_star, in_rhalf)
    galaxy.is_star_in_rhalf = is_star_in_rhalf
    
    # compute COMV as mass weighted average of stars in rhalf_fac*rhalf
    vel_in_rhalf = sn.part4['Velocities'][is_star_in_rhalf]
    mass_in_rhalf = sn.part4['Masses'][is_star_in_rhalf]
    COMV = np.average(vel_in_rhalf, axis=0, weights=mass_in_rhalf)
    
    vel = sn.part4['Velocities'] - COMV
    
    angmom = np.cross(pos[is_star_in_rhalf], vel[is_star_in_rhalf])
    angmom *= mass_in_rhalf.reshape(-1, 1)
    angmom = np.sum(angmom, axis=0)
    
    galaxy.CenterOfMass = EmptyClass()
    galaxy.CenterOfMass.Coordinate = COM
    galaxy.CenterOfMass.Velocity = COMV
    galaxy.CenterOfMass.AngularMomentum = angmom
    
    angmom_dir = angmom/np.linalg.norm(angmom)
    theta = np.arccos(np.dot(angmom_dir, np.array([0, 0, 1])))
    k = np.cross(angmom, np.array([0, 0, 1.]))
    k /= np.linalg.norm(k)

    galaxy.CenterOfMass.theta = theta
    galaxy.CenterOfMass.k = k
    
    for ptype in [0, 1, 4, 5]:
        if ptype==1 and load_DM is False:
            continue
        
        part = getattr(sn, 'part'+str(ptype))
        pos = part['Coordinates'] - COM
        vel = part['Velocities'] - COMV
    
        pos_rot = rodrigues_formula(k, pos.astype(np.float64), theta)
        vel_rot = rodrigues_formula(k, vel.astype(np.float64), theta)
        
        part['GalacticCoordinates'] = pos_rot
        part['GalacticVelocities'] = vel_rot
    
        part['GalacticCoordinatesPhys'] = pos_rot * a / h
        part['GalacticVelocitiesPhys'] = vel_rot * np.sqrt(a)
        
        if 'GFM_StellarFormationTime' in part.keys():
            part['FormationTimeGyr'] = get_time(part['GFM_StellarFormationTime'])
    
    galaxy.sn = sn
    
    return galaxy
